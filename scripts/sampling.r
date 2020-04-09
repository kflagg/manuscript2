######################
# Setup environment. #
######################

# Load packages.
library(parallel)
library(spatstat)
library(INLA)

# Create a cluster for applying in parallel.
cl <- makeCluster(ceiling(0.75) * detectCores())
invisible(clusterEvalQ(cl, {
  library(spatstat)
  library(maptools)
  library(INLA)
}))

# Define user-specified parameters.
XSECT_WIDTH <- 2 # Width of transects.
XSECT_LENGTH_MIN <- 50
XSECT_LENGTH_MAX <- 500
XSECT_SEP_MIN <- 10
XSECT_SEP_MAX <- 250
DIST_INITIAL <- 5000
DIST_MAX <- 50000
WAYPOINT_MARGIN <- XSECT_WIDTH / 2 # Minimum distance between waypoints and site boundary

# Mesh parameters to experiment with.
MAX_EDGE_LENGTH <- 100
MAX_EDGE_EXT <- 200
MARGIN <- 200

# Graphics parameters.
NPIX_X <- 500
NPIX_Y <- 500


# Define derived parameters
XSECT_RADIUS <- XSECT_WIDTH / 2 # Half the transect width
XSECT_NUM_INITIAL <- DIST_INITIAL / XSECT_LENGTH_MAX


###################################
# Objects pertaining to the site. #
###################################

# Site window.

sim_R <- owin(poly = cbind(
    x = c(   0, 1600, 1600, 2000, 2000,  500,    0),
    y = c(   0,    0,  250,  250, 2000, 2000, 1500)
  ), unitnames = c('meter', 'meters'))

# Mesh covering the site.
R_boundary <- inla.mesh.segment(loc = do.call(cbind, vertices.owin(sim_R)))
R_full_mesh <- inla.mesh.create(
  boundary = R_boundary,
  refine = list(max.edge = MAX_EDGE_LENGTH)
)

# Mesh including a margin outside the site.
R_margin_mesh <- inla.mesh.2d(
  loc = R_full_mesh$loc[,1:2], # Include nodes from site.
  offset = MARGIN,
  max.edge = MAX_EDGE_EXT # Fill in the rest with a coarser triangulation.
)

# Convert the mesh edges to a psp.
R_mesh_loc <- R_margin_mesh$loc[,-3]
R_mesh_adj <- R_margin_mesh$graph$vv
R_mesh_adj[lower.tri(R_mesh_adj)] <- 0
R_mesh_seg_idx0 <- do.call(c,
  parApply(cl, cbind(seq_len(ncol(R_mesh_adj)), apply(R_mesh_adj, 2, sum)), 1,
    function(x){return(rep(x[1], x[2]))})
  )
R_mesh_seg_idx1 <- do.call(c, parApply(cl, R_mesh_adj == 1, 2, which))
R_mesh_seg0 <- R_mesh_loc[R_mesh_seg_idx0,]
R_mesh_seg1 <- R_mesh_loc[R_mesh_seg_idx1,]
R_mesh_psp <- psp(
  x0 = R_mesh_seg0[,1],
  y0 = R_mesh_seg0[,2],
  x1 = R_mesh_seg1[,1],
  y1 = R_mesh_seg1[,2],
  window = owin(range(R_mesh_loc[,1]), range(R_mesh_loc[,2])))
R_mesh_win <- convexhull(R_mesh_psp)
Window(R_mesh_psp) <- R_mesh_win

# Convert mesh triangles to owins.
clusterExport(cl, c('R_mesh_loc', 'sim_R'))
R_mesh_tris <- parApply(cl, R_margin_mesh$graph$tv, 1, function(x){
  return(owin(poly = R_mesh_loc[x,]))
})
R_tri_areas <- parSapply(cl, parLapply(cl, R_mesh_tris, '[', sim_R), area)
R_tri_margin_areas <- parSapply(cl, R_mesh_tris, area)

# Calculate the area represented by each node.
R_nodes_margin_area <- diag(inla.mesh.fem(R_margin_mesh)$c0)

# Calculate the interior area represented by each node.
# WRONG LENGTH! Need to include 0s for nodes in the margin.
R_nodes_area <- diag(inla.mesh.fem(R_full_mesh)$c0)

# Set up the spatial numerical integration.
margin_nV <- margin_mesh$n
margin_IntegrationMatrix <- sparseMatrix(i = 1:margin_nV, j = 1:margin_nV, x = rep(1, margin_nV))
margin_IntegrationWeights <- diag(inla.mesh.fem(margin_mesh)$c0)

# Projector to interpolate over site.
margin_proj <- inla.mesh.projector(margin_mesh, dims = c(NPIX_X, NPIX_Y))


###################################
# Objects pertaining to the data. #
###################################

# Function to simulate one realized point pattern. This includes an LGCP
# background process with Matern (alpha = 2) covariance and a foreground
# cluster process. The true intensity function is stored in an attribute
# called Lambda.
# Dots should include dimyx = c(NPIX_Y, NPIX_X).
sim_data <- function(
  R = sim_R,
  matern_mu = log(1000 / area(R)), matern_sd = 1, matern_range = 200,
  thomas_kappa = 2 / area(R), thomas_scale = 100, thomas_mu = 250,
  ...){

  # Simulate a background LGCP.
  maternlgcp <- rLGCP(
    model = 'matern',
    mu = matern_mu,
    param = list(nu = 1, var = matern_sd^2, scale = matern_range),
    win = R, saveLambda = TRUE, nsim = 1, ...
  )

  # Simulate foreground hotspots as a Thomas cluster process.
  thomas <- rThomas(
    kappa = thomas_kappa,
    scale = thomas_scale,
    mu = thomas_mu,
    win = R, saveLambda = TRUE, nsim = 1, expand = 0, ...
  )

  # Superimpose the two processes.
  result <- superimpose(maternlgcp, thomas)

  # Add the realized intensity functions.
  attr(result, 'Lambda') <- attr(maternlgcp, 'Lambda') + attr(thomas, 'Lambda')

  return(result)
}


####################################################################
# Objects pertaining to the model but independent of survey plans. #
####################################################################

# Define an SPDE representation of the spatial GP using INLA's default priors.
margin_spde <- inla.spde2.matern(mesh = margin_mesh)

# Model formula with no covariates.
bei_formula <- y ~ -1 + intercept + f(idx, model = margin_spde)


###################
# Fully surveyed. #
###################

# Weights are propotion of area represented by nodes inside the site.
full_mesh_weights <- mesh_area / mesh_margin_area

# Organize variables for inla().
full_pts <- as.matrix(as.data.frame(bei))

# Contruct the SPDE A matrix for nodes and points.
full_nData <- nrow(full_pts)
full_LocationMatrix <- inla.mesh.project(margin_mesh, full_pts)$A
full_ObservationMatrix <- rbind(margin_IntegrationMatrix, full_LocationMatrix)

# Get the integration weights.
full_E <- c(full_mesh_weights * margin_IntegrationWeights, rep(0, full_nData))

# Create the psuedodata.
full_psuedodata <- c(rep(0, margin_nV), rep(1, full_nData))
full_data_list <- list(
  y = full_psuedodata,
  idx = 1:margin_nV,
  intercept = rep(1, margin_nV)
)

# Fit the model.
system.time(
result_full <- inla(
  formula = bei_formula,
  data = full_data_list,
  family = 'poisson',
  control.predictor = list(A = full_ObservationMatrix),
  control.compute = list(config = TRUE),
  E = full_E,
  verbose = TRUE
)
)
print(result_full$summary.fixed)
print(result_full$summary.hyperpar)

# Plot surface.
dev.new()
par(mar = c(0.5, 0, 2, 2))
plot(im(t(inla.mesh.project(margin_proj, result_full$summary.fixed$mean + result_full$summary.random$idx$mean)),
        xrange = Frame(bei)$x + c(-MARGIN, MARGIN),
        yrange = Frame(bei)$y + c(-MARGIN, MARGIN),
        unitname = c('meter', 'meters')),
        main = 'Posterior Mean of Log-Intensity')
plot(bei_win, border = '#808080ff', add = TRUE)
points(bei, pch = '.', col = 'white')

# Plot uncertainty.
dev.new()
par(mar = c(0.5, 0, 2, 2))
plot(im(t(sqrt(result_full$summary.fixed$sd^2 + inla.mesh.project(margin_proj, result_full$summary.random$idx$sd^2))),
        xrange = Frame(bei)$x + c(-MARGIN, MARGIN),
        yrange = Frame(bei)$y + c(-MARGIN, MARGIN),
        unitname = c('meter', 'meters')),
        main = 'Posterior SD of Log-Intensity')
plot(bei_win, border = '#808080ff', add = TRUE)
points(bei, pch = '.', col = 'white')
# SD or variance is also approximately linear within each triangle.
# A strange assumption?
# Note the interpolated SD or variance is not the SD/variance of the
# interpolated mean surface. We would need the covariances between the nodes
# to get that.

# Spatially-averaged random effect variance.
full_avgvar <- (result_full$summary.random$idx$sd^2) %*% mesh_area / sum(mesh_area)
print(full_avgvar)

# Max random effect variance at any node. (Interpolated variance cannot exceed this.)
full_maxvar <- max(result_full$summary.random$idx$sd[mesh_area > 0]^2)
print(full_maxvar)


#######################################################
# Sampling parameters applicable to all survey plans. #
#######################################################

# Define the region where waypoints are allowed
sampleable <- erosion(bei_win, WAYPOINT_MARGIN)
min_x <- min(vertices(sampleable)$x)
max_x <- max(vertices(sampleable)$x)
min_y <- min(vertices(sampleable)$y)
max_y <- max(vertices(sampleable)$y)


###################################
# Initial SRS and sequential SRS. #
###################################

# Start with SRS.
set.seed(7352)
srs_x <- sort(runif(XSECT_NUM_INITIAL, min_x, max_x))
waypoints <- cbind(x = rep(srs_x, each = 2), y = c(min_y, max_y, max_y, min_y))
n_waypoints <- nrow(waypoints)

cog_psp <- psp(
    x0 = waypoints[-n_waypoints, 'x'],
    y0 = waypoints[-n_waypoints, 'y'],
    x1 = waypoints[-1, 'x'],
    y1 = waypoints[-1, 'y'],
    window = Window(bei)
  )
observed_window <- dilation(cog_psp, XSECT_RADIUS)
observed_ppp <- bei[observed_window]


# Get numerical integration weights.
#  - Get center and length of segment in each triangle.
#  - Get barycentric coordinate of each center.
#  - Multiply barycentric coordinates by length and width in each triangle?
#  - Sum products for each node.

# Take the intersection of cog_psp with each triangle.
clusterExport(cl, c('cog_psp'))
cog_subsegs <- parLapply(cl, mesh_tris, function(x){return(cog_psp[x])})

# Track which triangle each segment came from.
clusterExport(cl, 'cog_subsegs')
seg_tri_idx <- unlist(parLapply(cl, seq_along(cog_subsegs), function(i){return(rep(i, cog_subsegs[[i]]$n))}))

# Get the midpoints
cog_midpoints <- as.matrix(do.call(rbind, parLapply(cl, parLapply(cl, cog_subsegs, midpoints.psp), as.data.frame)))

# Get barycentric coordinates of the midpoints.
cog_bary <- inla.mesh.projector(margin_mesh, loc = cog_midpoints)$proj$A

# Get the lengths within each triangle.
cog_lengths <- unlist(parSapply(cl, cog_subsegs, lengths.psp))

# Calculate the integration weights.
# For each row of the barycentric coordinates matrix (which is number of
# segments by number of nodes), divide each entry by the portion of the
# triangle's area represented by that node, which is 1/3rd of the area
# of the triangle the segment is in.
clusterExport(cl, c('cog_bary', 'tri_margin_areas', 'seg_tri_idx'))
cog_bary_prop <- as(t(parSapply(cl, seq_along(seg_tri_idx),
    function(i){return(cog_bary[i,] / tri_margin_areas[seg_tri_idx[i]] * 3)})),
  'sparseMatrix')
# Multiply by the observed area represented in each triangle.
srs_mesh_weights <- as.vector(XSECT_WIDTH * cog_lengths %*% cog_bary_prop)


# Organize variables for inla().
srs_pts <- as.matrix(as.data.frame(observed_ppp))

# Contruct the SPDE A matrix for nodes and points.
srs_nData <- nrow(srs_pts)
srs_LocationMatrix <- inla.mesh.project(margin_mesh, srs_pts)$A
srs_ObservationMatrix <- rbind(margin_IntegrationMatrix, srs_LocationMatrix)

# Get the integration weights.
srs_E <- c(srs_mesh_weights * margin_IntegrationWeights, rep(0, srs_nData))

# Create the psuedodata.
srs_psuedodata <- c(rep(0, margin_nV), rep(1, srs_nData))
srs_data_list <- list(
  y = srs_psuedodata,
  idx = 1:margin_nV,
  intercept = rep(1, margin_nV)
)

system.time(
result_srs <- inla(
  formula = bei_formula,
  data = srs_data_list,
  family = 'poisson',
  control.predictor = list(A = srs_ObservationMatrix),
  control.compute = list(config = TRUE),
  E = srs_E,
  verbose = TRUE
)
)
print(result_srs$summary.fixed)
print(result_srs$summary.hyperpar)

# Plot surface.
dev.set(dev.prev())
par(mar = c(0.5, 0, 2, 2))
plot(im(t(inla.mesh.project(margin_proj, result_srs$summary.fixed$mean + result_srs$summary.random$idx$mean)),
        xrange = Frame(bei)$x + c(-MARGIN, MARGIN),
        yrange = Frame(bei)$y + c(-MARGIN, MARGIN),
        unitname = c('meter', 'meters')),
        main = 'Posterior Mean of Log-Intensity')
plot(bei_win, border = '#808080ff', add = TRUE)
plot(Window(observed_ppp), border = '#80808080', add = TRUE)
points(observed_ppp, pch = '.', col = 'white')

# Plot uncertainty.
dev.set(dev.next())
par(mar = c(0.5, 0, 2, 2))
plot(im(t(sqrt(result_srs$summary.fixed$sd^2 + inla.mesh.project(margin_proj, result_srs$summary.random$idx$sd^2))),
        xrange = Frame(bei)$x + c(-MARGIN, MARGIN),
        yrange = Frame(bei)$y + c(-MARGIN, MARGIN),
        unitname = c('meter', 'meters')),
        main = 'Posterior SD of Log-Intensity')
plot(bei_win, border = '#808080ff', add = TRUE)
plot(Window(observed_ppp), border = '#80808080', add = TRUE)
points(observed_ppp, pch = '.', col = 'white')

# Spatially-averaged random effect variance.
srs_avgvar <- (result_srs$summary.random$idx$sd^2) %*% mesh_area / sum(mesh_area)
print(srs_avgvar)
# Max random effect variance at any node. (Interpolated variance cannot exceed this.)
srs_maxvar <- max(result_srs$summary.random$idx$sd[mesh_area > 0]^2)
print(srs_maxvar)
# Distance traveled.
srs_pathdist <- sum(cog_lengths)
print(srs_pathdist)


# Open two more graphics devices to compare new plots to old plots.
dev.new()
dev.new()


# Sample sequentially.
while(tail(srs_pathdist, n = 1) < DIST_MAX){


# Choose a new segment.
srs_x_new <- runif(1, min_x, max_x)
waypoints_new <- cbind(
  x = rep(srs_x_new, 2),
  y = rev(tail(waypoints[,'y'], n = 2))
)

# Create a path from the last waypoint through the new segment.
cog_new <- psp(
    x0 = c(waypoints[n_waypoints, 'x'], waypoints_new[1, 'x']),
    y0 = c(waypoints[n_waypoints, 'y'], waypoints_new[1, 'y']),
    x1 = waypoints_new[,'x'],
    y1 = waypoints_new[,'y'],
    window = Window(bei)
  )
observed_window_new <- dilation(cog_new, XSECT_RADIUS)
observed_ppp_new <- bei[observed_window_new]

# Triangulate the new segments.
clusterExport(cl, 'cog_new')
cog_subsegs_new <- parLapply(cl, mesh_tris, function(x){return(cog_new[x])})

# Track which triangle each segment came from.
clusterExport(cl, 'cog_subsegs_new')
seg_tri_idx_new <- unlist(parLapply(cl, seq_along(cog_subsegs_new), function(i){return(rep(i, cog_subsegs_new[[i]]$n))}))

# Get the midpoints
cog_midpoints_new <- as.matrix(do.call(rbind, parLapply(cl, parLapply(cl, cog_subsegs_new, midpoints.psp), as.data.frame)))

# Get barycentric coordinates of the midpoints.
cog_bary_new <- inla.mesh.projector(margin_mesh, loc = cog_midpoints_new)$proj$A

# Get the lengths within each triangle.
cog_lengths_new <- unlist(parSapply(cl, cog_subsegs_new, lengths.psp))

# Calculate the integration weights.
clusterExport(cl, c('cog_bary_new', 'seg_tri_idx_new'))
cog_bary_prop_new <- as(t(parSapply(cl, seq_along(seg_tri_idx_new),
    function(i){return(cog_bary_new[i,] / tri_margin_areas[seg_tri_idx_new[i]] * 3)})),
  'sparseMatrix')
# Multiply by the observed area represented in each triangle.
srs_mesh_weights <- srs_mesh_weights + as.vector(XSECT_WIDTH * cog_lengths_new %*% cog_bary_prop_new)


# Combine old and new segments.
srs_x <- c(srs_x, srs_x_new)
waypoints <- rbind(waypoints, waypoints_new)
n_waypoints <- nrow(waypoints)
cog_psp <- psp(
    x0 = waypoints[-n_waypoints, 'x'],
    y0 = waypoints[-n_waypoints, 'y'],
    x1 = waypoints[-1, 'x'],
    y1 = waypoints[-1, 'y'],
    window = Window(bei)
  )
observed_window <- dilation(cog_psp, XSECT_RADIUS)
observed_ppp <- bei[observed_window]


# Organize variables for inla().
srs_pts <- rbind(srs_pts, as.matrix(as.data.frame(observed_ppp_new)))

# Contruct the SPDE A matrix for nodes and points.
srs_nData <- nrow(srs_pts)
srs_LocationMatrix <- inla.mesh.project(margin_mesh, srs_pts)$A
srs_ObservationMatrix <- rbind(margin_IntegrationMatrix, srs_LocationMatrix)

# Get the integration weights.
srs_E <- c(srs_mesh_weights * margin_IntegrationWeights, rep(0, srs_nData))

# Create the psuedodata.
srs_psuedodata <- c(rep(0, margin_nV), rep(1, srs_nData))

# Fit model to initial site.
srs_data_list <- list(
  y = srs_psuedodata,
  idx = 1:margin_nV,
  intercept = rep(1, margin_nV)
)

system.time(
result_srs <- inla(
  formula = bei_formula,
  data = srs_data_list,
  family = 'poisson',
  control.predictor = list(A = srs_ObservationMatrix),
  control.compute = list(config = TRUE),
  E = srs_E,
  verbose = TRUE
)
)
print(result_srs$summary.fixed)
print(result_srs$summary.hyperpar)

# Plot surface.
dev.set(dev.prev())
par(mar = c(0.5, 0, 2, 2))
plot(im(t(inla.mesh.project(margin_proj, result_srs$summary.fixed$mean + result_srs$summary.random$idx$mean)),
        xrange = Frame(bei)$x + c(-MARGIN, MARGIN),
        yrange = Frame(bei)$y + c(-MARGIN, MARGIN),
        unitname = c('meter', 'meters')),
        main = 'Posterior Mean of Log-Intensity')
plot(bei_win, border = '#808080ff', add = TRUE)
plot(Window(observed_ppp), border = '#80808080', add = TRUE)
points(observed_ppp, pch = '.', col = 'white')

# Plot uncertainty.
dev.set(dev.next())
par(mar = c(0.5, 0, 2, 2))
plot(im(t(sqrt(result_srs$summary.fixed$sd^2 + inla.mesh.project(margin_proj, result_srs$summary.random$idx$sd^2))),
        xrange = Frame(bei)$x + c(-MARGIN, MARGIN),
        yrange = Frame(bei)$y + c(-MARGIN, MARGIN),
        unitname = c('meter', 'meters')),
        main = 'Posterior SD of Log-Intensity')
plot(bei_win, border = '#808080ff', add = TRUE)
plot(Window(observed_ppp), border = '#80808080', add = TRUE)
points(observed_ppp, pch = '.', col = 'white')

# Spatially-averaged random effect variance.
srs_avgvar <- c(srs_avgvar, (result_srs$summary.random$idx$sd^2) %*% mesh_area / sum(mesh_area))
print(tail(srs_avgvar, n = 1))
# Max random effect variance at any node. (Interpolated variance cannot exceed this.)
srs_maxvar <- c(srs_maxvar, max(result_srs$summary.random$idx$sd[mesh_area > 0]^2))
print(tail(srs_maxvar, n = 1))
# Distance traveled.
srs_pathdist <- c(srs_pathdist, tail(srs_pathdist, n = 1) + sum(cog_lengths_new))
print(tail(srs_pathdist, n = 1))


} # End SRS SRS loop.


########################################
# Initial SRS and sequential perp SRS. #
########################################

# Start with SRS.
set.seed(352)
srs_perp_x <- sort(runif(XSECT_NUM_INITIAL, min_x, max_x))
waypoints <- cbind(x = rep(srs_perp_x, each = 2), y = c(min_y, max_y, max_y, min_y))
n_waypoints <- nrow(waypoints)

cog_psp <- psp(
    x0 = waypoints[-n_waypoints, 'x'],
    y0 = waypoints[-n_waypoints, 'y'],
    x1 = waypoints[-1, 'x'],
    y1 = waypoints[-1, 'y'],
    window = Window(bei)
  )
observed_window <- dilation(cog_psp, XSECT_RADIUS)
observed_ppp <- bei[observed_window]


# Get numerical integration weights.
#  - Get center and length of segment in each triangle.
#  - Get barycentric coordinate of each center.
#  - Multiply barycentric coordinates by length and width in each triangle?
#  - Sum products for each node.

# Take the intersection of cog_psp with each triangle.
clusterExport(cl, c('cog_psp'))
cog_subsegs <- parLapply(cl, mesh_tris, function(x){return(cog_psp[x])})

# Track which triangle each segment came from.
clusterExport(cl, 'cog_subsegs')
seg_tri_idx <- unlist(parLapply(cl, seq_along(cog_subsegs), function(i){return(rep(i, cog_subsegs[[i]]$n))}))

# Get the midpoints
cog_midpoints <- as.matrix(do.call(rbind, parLapply(cl, parLapply(cl, cog_subsegs, midpoints.psp), as.data.frame)))

# Get barycentric coordinates of the midpoints.
cog_bary <- inla.mesh.projector(margin_mesh, loc = cog_midpoints)$proj$A

# Get the lengths within each triangle.
cog_lengths <- unlist(parSapply(cl, cog_subsegs, lengths.psp))

# Calculate the integration weights.
# For each row of the barycentric coordinates matrix (which is number of
# segments by number of nodes), divide each entry by the portion of the
# triangle's area represented by that node, which is 1/3rd of the area
# of the triangle the segment is in.
clusterExport(cl, c('cog_bary', 'tri_margin_areas', 'seg_tri_idx'))
cog_bary_prop <- as(t(parSapply(cl, seq_along(seg_tri_idx),
    function(i){return(cog_bary[i,] / tri_margin_areas[seg_tri_idx[i]] * 3)})),
  'sparseMatrix')
# Multiply by the observed area represented in each triangle.
srs_perp_mesh_weights <- as.vector(XSECT_WIDTH * cog_lengths %*% cog_bary_prop)


# Organize variables for inla().
srs_perp_pts <- as.matrix(as.data.frame(observed_ppp))

# Contruct the SPDE A matrix for nodes and points.
srs_perp_nData <- nrow(srs_perp_pts)
srs_perp_LocationMatrix <- inla.mesh.project(margin_mesh, srs_perp_pts)$A
srs_perp_ObservationMatrix <- rbind(margin_IntegrationMatrix, srs_perp_LocationMatrix)

# Get the integration weights.
srs_perp_E <- c(srs_perp_mesh_weights * margin_IntegrationWeights, rep(0, srs_perp_nData))

# Create the psuedodata.
srs_perp_psuedodata <- c(rep(0, margin_nV), rep(1, srs_perp_nData))
srs_perp_data_list <- list(
  y = srs_perp_psuedodata,
  idx = 1:margin_nV,
  intercept = rep(1, margin_nV)
)

system.time(
result_srs_perp <- inla(
  formula = bei_formula,
  data = srs_perp_data_list,
  family = 'poisson',
  control.predictor = list(A = srs_perp_ObservationMatrix),
  control.compute = list(config = TRUE),
  E = srs_perp_E,
  verbose = TRUE
)
)
print(result_srs_perp$summary.fixed)
print(result_srs_perp$summary.hyperpar)

# Plot surface.
dev.set(dev.prev())
par(mar = c(0.5, 0, 2, 2))
plot(im(t(inla.mesh.project(margin_proj, result_srs_perp$summary.fixed$mean + result_srs_perp$summary.random$idx$mean)),
        xrange = Frame(bei)$x + c(-MARGIN, MARGIN),
        yrange = Frame(bei)$y + c(-MARGIN, MARGIN),
        unitname = c('meter', 'meters')),
        main = 'Posterior Mean of Log-Intensity')
plot(bei_win, border = '#808080ff', add = TRUE)
plot(Window(observed_ppp), border = '#80808080', add = TRUE)
points(observed_ppp, pch = '.', col = 'white')

# Plot uncertainty.
dev.set(dev.next())
par(mar = c(0.5, 0, 2, 2))
plot(im(t(sqrt(result_srs_perp$summary.fixed$sd^2 + inla.mesh.project(margin_proj, result_srs_perp$summary.random$idx$sd^2))),
        xrange = Frame(bei)$x + c(-MARGIN, MARGIN),
        yrange = Frame(bei)$y + c(-MARGIN, MARGIN),
        unitname = c('meter', 'meters')),
        main = 'Posterior SD of Log-Intensity')
plot(bei_win, border = '#808080ff', add = TRUE)
plot(Window(observed_ppp), border = '#80808080', add = TRUE)
points(observed_ppp, pch = '.', col = 'white')

# Spatially-averaged random effect variance.
srs_perp_avgvar <- (result_srs_perp$summary.random$idx$sd^2) %*% mesh_area / sum(mesh_area)
print(srs_perp_avgvar)
# Max random effect variance at any node. (Interpolated variance cannot exceed this.)
srs_perp_maxvar <- max(result_srs_perp$summary.random$idx$sd[mesh_area > 0]^2)
print(srs_perp_maxvar)
# Distance traveled.
srs_perp_pathdist <- sum(cog_lengths)
print(srs_perp_pathdist)


# Open two more graphics devices to compare new plots to old plots.
dev.new()
dev.new()


# Sample sequentially.
while(tail(srs_perp_pathdist, n = 1) < DIST_MAX){


# Choose a new segment.
srs_perp_y_new <- runif(1, min_y, max_y)
waypoints_new <- cbind(
  x = c(min_x, max_x)[order(abs(tail(waypoints[,'x'], n = 1) - c(min_x, max_x)))],
  y = rep(srs_perp_y_new, 2)
)

# Create a path from the last waypoint through the new segment.
cog_new <- psp(
    x0 = c(waypoints[n_waypoints, 'x'], waypoints_new[1, 'x']),
    y0 = c(waypoints[n_waypoints, 'y'], waypoints_new[1, 'y']),
    x1 = waypoints_new[,'x'],
    y1 = waypoints_new[,'y'],
    window = Window(bei)
  )
observed_window_new <- dilation(cog_new, XSECT_RADIUS)
observed_ppp_new <- bei[observed_window_new]

# Triangulate the new segments.
clusterExport(cl, 'cog_new')
cog_subsegs_new <- parLapply(cl, mesh_tris, function(x){return(cog_new[x])})

# Track which triangle each segment came from.
clusterExport(cl, 'cog_subsegs_new')
seg_tri_idx_new <- unlist(parLapply(cl, seq_along(cog_subsegs_new), function(i){return(rep(i, cog_subsegs_new[[i]]$n))}))

# Get the midpoints
cog_midpoints_new <- as.matrix(do.call(rbind, parLapply(cl, parLapply(cl, cog_subsegs_new, midpoints.psp), as.data.frame)))

# Get barycentric coordinates of the midpoints.
cog_bary_new <- inla.mesh.projector(margin_mesh, loc = cog_midpoints_new)$proj$A

# Get the lengths within each triangle.
cog_lengths_new <- unlist(parSapply(cl, cog_subsegs_new, lengths.psp))

# Calculate the integration weights.
clusterExport(cl, c('cog_bary_new', 'seg_tri_idx_new'))
cog_bary_prop_new <- as(t(parSapply(cl, seq_along(seg_tri_idx_new),
    function(i){return(cog_bary_new[i,] / tri_margin_areas[seg_tri_idx_new[i]] * 3)})),
  'sparseMatrix')
# Multiply by the observed area represented in each triangle.
srs_perp_mesh_weights <- srs_perp_mesh_weights + as.vector(XSECT_WIDTH * cog_lengths_new %*% cog_bary_prop_new)


# Combine old and new segments.
waypoints <- rbind(waypoints, waypoints_new)
n_waypoints <- nrow(waypoints)
cog_psp <- psp(
    x0 = waypoints[-n_waypoints, 'x'],
    y0 = waypoints[-n_waypoints, 'y'],
    x1 = waypoints[-1, 'x'],
    y1 = waypoints[-1, 'y'],
    window = Window(bei)
  )
observed_window <- dilation(cog_psp, XSECT_RADIUS)
observed_ppp <- bei[observed_window]


# Organize variables for inla().
srs_perp_pts <- rbind(srs_perp_pts, as.matrix(as.data.frame(observed_ppp_new)))

# Contruct the SPDE A matrix for nodes and points.
srs_perp_nData <- nrow(srs_perp_pts)
srs_perp_LocationMatrix <- inla.mesh.project(margin_mesh, srs_perp_pts)$A
srs_perp_ObservationMatrix <- rbind(margin_IntegrationMatrix, srs_perp_LocationMatrix)

# Get the integration weights.
srs_perp_E <- c(srs_perp_mesh_weights * margin_IntegrationWeights, rep(0, srs_perp_nData))

# Create the psuedodata.
srs_perp_psuedodata <- c(rep(0, margin_nV), rep(1, srs_perp_nData))

# Fit model to initial site.
srs_perp_data_list <- list(
  y = srs_perp_psuedodata,
  idx = 1:margin_nV,
  intercept = rep(1, margin_nV)
)

system.time(
result_srs_perp <- inla(
  formula = bei_formula,
  data = srs_perp_data_list,
  family = 'poisson',
  control.predictor = list(A = srs_perp_ObservationMatrix),
  control.compute = list(config = TRUE),
  E = srs_perp_E,
  verbose = TRUE
)
)
print(result_srs_perp$summary.fixed)
print(result_srs_perp$summary.hyperpar)

# Plot surface.
dev.set(dev.prev())
par(mar = c(0.5, 0, 2, 2))
plot(im(t(inla.mesh.project(margin_proj, result_srs_perp$summary.fixed$mean + result_srs_perp$summary.random$idx$mean)),
        xrange = Frame(bei)$x + c(-MARGIN, MARGIN),
        yrange = Frame(bei)$y + c(-MARGIN, MARGIN),
        unitname = c('meter', 'meters')),
        main = 'Posterior Mean of Log-Intensity')
plot(bei_win, border = '#808080ff', add = TRUE)
plot(Window(observed_ppp), border = '#80808080', add = TRUE)
points(observed_ppp, pch = '.', col = 'white')

# Plot uncertainty.
dev.set(dev.next())
par(mar = c(0.5, 0, 2, 2))
plot(im(t(sqrt(result_srs_perp$summary.fixed$sd^2 + inla.mesh.project(margin_proj, result_srs_perp$summary.random$idx$sd^2))),
        xrange = Frame(bei)$x + c(-MARGIN, MARGIN),
        yrange = Frame(bei)$y + c(-MARGIN, MARGIN),
        unitname = c('meter', 'meters')),
        main = 'Posterior SD of Log-Intensity')
plot(bei_win, border = '#808080ff', add = TRUE)
plot(Window(observed_ppp), border = '#80808080', add = TRUE)
points(observed_ppp, pch = '.', col = 'white')

# Spatially-averaged random effect variance.
srs_perp_avgvar <- c(srs_perp_avgvar, (result_srs_perp$summary.random$idx$sd^2) %*% mesh_area / sum(mesh_area))
print(tail(srs_perp_avgvar, n = 1))
# Max random effect variance at any node. (Interpolated variance cannot exceed this.)
srs_perp_maxvar <- c(srs_perp_maxvar, max(result_srs_perp$summary.random$idx$sd[mesh_area > 0]^2))
print(tail(srs_perp_maxvar, n = 1))
# Distance traveled.
srs_perp_pathdist <- c(srs_perp_pathdist, tail(srs_perp_pathdist, n = 1) + sum(cog_lengths_new))
print(tail(srs_perp_pathdist, n = 1))


} # End SRS perp SRS loop.


# Also do regular starting points and adapt by going through highest variance.
# And transpose plans and use shorter sequential xsects.


stopCluster(cl)
