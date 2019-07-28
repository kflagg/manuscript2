# Load packages.
library(spatstat)
library(maptools)
library(INLA)

# Define user-specified parameters.
XSECT_WIDTH <- 2 # Width of transects.
XSECT_LENGTH_MIN <- 50
XSECT_LENGTH_MAX <- 500
XSECT_SEP_MIN <- 10
XSECT_SEP_MAX <- 250
DIST_INITIAL <- 5000
DIST_MAX <- 20000
WAYPOINT_MARGIN <- XSECT_WIDTH # Minimum distance between waypoints and site boundary
NPIX_X <- 400
NPIX_Y <- 200

# Mesh parameters to experiment with.
MAX_EDGE_LENGTH <- 25
MAX_EDGE_EXT <- 50
MARGIN <- 100


# Define derived parameters
XSECT_RADIUS <- XSECT_WIDTH / 2 # Half the transect width
XSECT_NUM_INITIAL <- DIST_INITIAL / XSECT_LENGTH_MAX


# Mesh covering the site.
bei_boundary <- inla.mesh.segment(loc = do.call(cbind, vertices.owin(Window(bei))))
bei_full_mesh <- inla.mesh.create(
  boundary = bei_boundary,
  refine = list(max.edge = MAX_EDGE_LENGTH)
)

# Mesh including a margin outside the site.
margin_mesh <- inla.mesh.2d(
  loc = bei_full_mesh$loc[,1:2], # Include nodes from site.
  offset = MARGIN,
  max.edge = MAX_EDGE_EXT # Fill in the rest with a coarser triangulation.
)
margin_spde <- inla.spde2.matern(mesh = margin_mesh)


# Define the region where waypoints are allowed
sampleable <- erosion(Window(bei), WAYPOINT_MARGIN)
min_x <- min(vertices(sampleable)$x)
max_x <- max(vertices(sampleable)$x)
min_y <- min(vertices(sampleable)$y)
max_y <- max(vertices(sampleable)$y)

# Set up initial transect plan.
#set.seed(3829)
#seg_x <- runif(1, min_x, min_x + XSECT_SEP_MAX)
#while(max(seg_x) < max_x){
#  seg_x <- c(seg_x, max(seg_x) + runif(1, min_x, min_x + XSECT_SEP_MAX))
#}
#seg_x <- seg_x[seg_x < max_x]
#waypoints <- cbind(x = rep(seg_x, each = 2), y = c(min_y, max_y, max_y, min_y))
#n_waypoints <- nrow(waypoints)

# Start with SRS.
set.seed(7352)
seg_x <- sort(runif(XSECT_NUM_INITIAL, min_x, max_x))
waypoints <- cbind(x = rep(seg_x, each = 2), y = c(min_y, max_y, max_y, min_y))
n_waypoints <- nrow(waypoints)

#plot(bei)
#lines(waypoints)


cog_psp <- psp(
    x0 = waypoints[-n_waypoints, 'x'],
    y0 = waypoints[-n_waypoints, 'y'],
    x1 = waypoints[-1, 'x'],
    y1 = waypoints[-1, 'y'],
    window = Window(bei)
  )
observed_window <- dilation(cog_psp, XSECT_RADIUS)
observed_ppp <- bei[observed_window]

#plot(margin_mesh, asp = 1)
#plot(Window(observed_ppp), col = '#80000080', border = NA, add = TRUE)
#points(observed_ppp, pch = 4, cex = 0.25, col = 'red')


# TODO: Get numerical integration weights.
#  - Get center and length of segment in each triangle.
#  - Get barycentric coordinate of each center.
#  - Multiply barycentric coordinates by length and width in each triangle?
#  - Sum products for each node.

# Convert the mesh edges to a psp.
meshloc <- margin_mesh$loc[,-3]
meshadj <- margin_mesh$graph$vv
meshadj[lower.tri(meshadj)] <- 0
meshsegidx0 <- do.call(c,
  apply(cbind(seq_len(ncol(meshadj)), apply(meshadj, 2, sum)), 1,
    function(x){return(rep(x[1], x[2]))})
  )
meshsegidx1 <- do.call(c, apply(meshadj == 1, 2, which))
meshseg0 <- meshloc[meshsegidx0,]
meshseg1 <- meshloc[meshsegidx1,]
mesh_psp <- psp(
  x0 = meshseg0[,1],
  y0 = meshseg0[,2],
  x1 = meshseg1[,1],
  y1 = meshseg1[,2],
  window = owin(range(meshloc[,1]), range(meshloc[,2])))
mesh_win <- convexhull(mesh_psp)
Window(mesh_psp) <- mesh_win

# Take the intersection of cog_psp with each triangle.
mesh_tris <- apply(margin_mesh$graph$tv, 1, function(x){
  return(owin(poly = meshloc[x,]))
})
tri_areas <- sapply(mesh_tris, area)
cog_subsegs <- lapply(mesh_tris, `[.psp`, x = cog_psp)

# Join the segments and track which triangle each came from.
cog_split <- `Window<-`(cog_subsegs[[1]], value = mesh_win)
seg_tri_idx <- rep(1L, cog_subsegs[[1]]$n)
for(i in 2:length(cog_subsegs)){
  cog_split <- append.psp(cog_split, `Window<-`(cog_subsegs[[i]], value = mesh_win))
  seg_tri_idx <- c(seg_tri_idx, rep(i, cog_subsegs[[i]]$n))
}

# Get the midpoints
cog_midpoints <- as.matrix(as.data.frame(midpoints.psp(cog_split)))

# Get barycentric coordinates of the midpoints.
cog_bary <- inla.mesh.projector(margin_mesh, loc = cog_midpoints)$proj$A

# Get the lengths within each triangle.
cog_lengths <- lengths(cog_split)

# Calculate the integration weights.
# For each row of the barycentric coordinates matrix (which is number of
# segments by number of nodes), divide each entry by the portion of the
# triangle's area represented by that node, which is 1/3rd of the area
# of the triangle the segment is in.
cog_bary_prop <- array(NA_real_, dim = dim(cog_bary))
for(i in seq_len(nrow(cog_bary_prop))){
  cog_bary_prop[i,] <- cog_bary[i,] / tri_areas[seg_tri_idx[i]] * 3
}
cog_bary_prop <- as(cog_bary_prop, 'sparseMatrix')
# Multiply by the observed area represented in each triangle.
mesh_weights <- as.vector(XSECT_WIDTH * cog_lengths %*% cog_bary_prop)

#plot(margin_mesh, asp = 1)
#for(i in seq_along(cog_subsegs)){
#  plot(Window(cog_subsegs[[i]]), col = 'grey', border = 'darkgrey', add = TRUE)
#  plot(cog_subsegs[[i]], col = 'red', add = TRUE)
#}
#points(cog_midpoints, col = 'red', pch = 16, cex = 0.5)


# Organize variables for inla().
samp_pts <- cbind(observed_ppp$x, observed_ppp$y)
proj_margin_samp <- inla.mesh.projector(margin_mesh, dims = c(NPIX_X, NPIX_Y))

# Contruct the SPDE A matrix for nodes and points.
samp_nV <- margin_mesh$n
samp_nData <- dim(samp_pts)[1]
samp_LocationMatrix <- inla.mesh.project(margin_mesh, samp_pts)$A
samp_IntegrationMatrix <- sparseMatrix(i = 1:samp_nV, j = 1:samp_nV, x = rep(1, samp_nV))
samp_ObservationMatrix <- rbind(samp_IntegrationMatrix, samp_LocationMatrix)

# Get the integration weights.
samp_IntegrationWeights <- diag(inla.mesh.fem(margin_mesh)$c0)
samp_E_point_process <- c(mesh_weights * samp_IntegrationWeights, rep(0, samp_nData))

# Create the psuedodata.
samp_fake_data <- c(rep(0, samp_nV), rep(1, samp_nData))

# Fit model to initial site.
samp_formula <- y ~ -1 + intercept + f(idx, model = margin_spde) # No covariates.
samp_data <- list(y = samp_fake_data, idx = 1:samp_nV, intercept = rep(1, samp_nV))

system.time(
result_samp <- inla(
  formula = samp_formula,
  data = samp_data,
  family = 'poisson',
  control.predictor = list(A = samp_ObservationMatrix),
  E = samp_E_point_process,
  verbose = TRUE
)
)
result_samp$summary.fixed
result_samp$summary.hyperpar

# Plot surface.
#par(mar = c(0.5, 0, 2, 2))
#plot(im(t(inla.mesh.project(proj_margin_samp, result_samp$summary.fixed$mean + result_samp$summary.random$idx$mean)),
#        xrange = Frame(bei)$x + c(-MARGIN, MARGIN),
#        yrange = Frame(bei)$y + c(-MARGIN, MARGIN),
#        unitname = c('meter', 'meters')),
#        main = 'Posterior Mean of Log-Intensity')
#plot(Window(bei), border = 'white', add = TRUE)
#plot(Window(observed_ppp), border = 'white', add = TRUE)
#points(observed_ppp, pch = '.', col = 'white')
#par(mar = c(0.5, 0, 2, 2))
#plot(im(t(inla.mesh.project(proj_margin_samp, sqrt(result_samp$summary.fixed$sd^2 + result_samp$summary.random$idx$sd^2))),
#        xrange = Frame(bei)$x + c(-MARGIN, MARGIN),
#        yrange = Frame(bei)$y + c(-MARGIN, MARGIN),
#        unitname = c('meter', 'meters')),
#        main = 'Posterior SD of Log-Intensity')
#plot(Window(bei), border = 'white', add = TRUE)
#plot(Window(observed_ppp), border = 'white', add = TRUE)
#points(observed_ppp, pch = '.', col = 'white')
