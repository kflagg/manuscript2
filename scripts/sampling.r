#########################
#{{{ Setup environment. #
#########################

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


#}}}##################################
#{{{ Objects pertaining to the site. #
######################################

# Site window.

sim_R <- owin(poly = cbind(
    x = c(   0, 1600, 1600, 2000, 2000,  500,    0),
    y = c(   0,    0,  250,  250, 2000, 2000, 1500)
  ), unitnames = c('meter', 'meters'))

# Mesh covering the site.
R_boundary <- inla.mesh.segment(loc = do.call(cbind, vertices.owin(sim_R)))
R_mesh <- inla.mesh.create(
  boundary = R_boundary,
  refine = list(max.edge = MAX_EDGE_LENGTH)
)

# Mesh including a margin outside the site.
#R_margin_mesh <- inla.mesh.2d(
#  loc = R_mesh$loc[,1:2], # Include nodes from site.
#  offset = MARGIN,
#  max.edge = MAX_EDGE_EXT # Fill in the rest with a coarser triangulation.
#)

# Get the mesh nodes.
R_mesh_loc <- R_mesh$loc[,1:2]

# Convert mesh triangles to owins and create a tesselation.
clusterExport(cl, c('R_mesh_loc', 'sim_R'))
R_mesh_tess <- as.tess(parApply(cl, R_mesh$graph$tv, 1, function(x){
  return(owin(poly = R_mesh_loc[x,]))
}))

# Get the area of each triangle.
#R_tri_margin_areas <- tile.areas(R_margin_mesh_tess)

# Get the interior area of each triangle.
R_tri_areas <- parSapply(cl, parLapply(cl, tiles(R_mesh_tess), '[', sim_R), area)

# Construct the dual of the mesh.
# Start by getting a list of triangles for each node.
clusterExport(cl, 'R_mesh')
node_tris <- parLapply(cl, parLapply(cl, parLapply(cl, parLapply(cl,
  # Indicate which triangles have this node as a vertex.
  seq_len(R_mesh$n), `==`, R_mesh$graph$tv),
  # Flatten the indicator matrix.
  rowSums),
  # Make the indicator vector logical instead of integer 0/1.
  as.logical),
  # Find which rows of the matrix (triangles) had this node.
  which)

# Get the centroids of each triangle.
tri_centroids <- t(parSapply(cl, tiles(R_mesh_tess), centroid.owin))

# Now, construct the dual mesh. Get the quadrilaterals made by cutting each
# triangle from its edge midpoints to its centroid, indexed by the mesh nodes.
# Then union the quadrilaterals for each node.
clusterExport(cl, c('node_tris', 'tri_centroids'))
dual_tess <- as.tess(parLapply(cl, seq_len(R_mesh$n), function(i){return(
  do.call(union.owin, lapply(node_tris[[i]], function(j){return(
    convexhull.xy(rbind(
      # jth triangle centroid.
      tri_centroids[j,],
      # Average of the ith node and all vertices of the jth triangle.
      # This results in the ith node and the midpoints of its edges.
      t((R_mesh_loc[i,] + t(R_mesh_loc[
        R_mesh$graph$tv[j,],
      1:2])) / 2)
    ))
  )}))
)}))

# Get the area represented by each node.
#R_nodes_margin_area <- diag(inla.mesh.fem(R_margin_mesh)$c0)

# Calculate the interior area represented by each node.
R_nodes_area <- parSapply(cl, parLapply(cl, tiles(dual_tess), `[`, sim_R), area)

# Will repeat for each sampling plan later, subseting using D.

# Neat plot.
plot(dual_tess, border = '#80808020', do.col = TRUE, values = R_nodes_area,
     main = 'Mesh with dual colored by node weight')
plot(R_mesh_tess, add = TRUE, border = '#00000080')
plot(sim_R, border = 'white', add = TRUE)
points(R_mesh_loc[,], pch = 20)

# Set up the spatial numerical integration.
R_mesh_size <- R_mesh$n
R_node_bary <- inla.mesh.project(R_mesh, R_mesh_loc)$A

# Projector to interpolate over site.
R_proj <- inla.mesh.projector(R_mesh, dims = c(NPIX_X, NPIX_Y))


#}}}##################################
#{{{ Objects pertaining to the data. #
######################################

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

sim_dataset <- list(sim_data(dimyx = c(NPIX_Y, NPIX_X)))
idx <- 1


#}}}###################################################################
#{{{ Objects pertaining to the model but independent of survey plans. #
#######################################################################

# Define an SPDE representation of the spatial GP using PC priors.
# TODO: update this with priors.
R_spde <- inla.spde2.matern(R_mesh)

# Model formula with no covariates.
R_formula <- y ~ -1 + intercept + f(idx, model = R_spde)


#}}}###############
# Fully surveyed. #
###################

# Observed event locations.
full_pts <- cbind(sim_dataset[[idx]]$x, sim_dataset[[idx]]$y)
full_n_events <- nrow(full_pts)

# Create the psuedodata. This is a vector giving the count of events at each
# pseudopoint location, that is 0 at each mesh node and 1 at each event.
full_pseudodata <- c(rep(0, R_mesh_size), rep(1, full_n_events))

# Concatenate the weight vector with a vector of zeros for the observed events.
# This is the vector of Poisson exposure parameters for the pseudodata.
full_pseudodata_exp <- c(R_nodes_area, rep(0, full_n_events))

# Compute the barycentric coordinates of the observed events
# (i.e. project into the space spanned by the basis of mesh nodes).
full_bary <- inla.mesh.project(R_mesh, full_pts)$A

# Bind the node and event coordinates into a single matrix of pseudodata
# locations in barycentric coordinates.
full_pseudopoints <- rbind(R_node_bary, full_bary)

# Create the data list to pass to inla().
# Indices and intercepts are only needed for the nodes.
full_inla_data <- list(
  y = full_pseudodata, # The whole pseudodata vector.
  idx = seq_len(R_mesh_size), # Indices of the nodes.
  intercept = rep(1, R_mesh_size) # Intercept column.
)

# Fit the model as a Poisson GLM with exposures specified.
full_result <- inla(
  formula = R_formula,
  data = full_inla_data,
  family = 'poisson',
  control.fixed = list(
    # Prior means and precisions for coefficients.
    mean.intercept = 0,
    prec.intercept = 0,
    mean = 0,
    prec = 0
  ),
  control.predictor = list(A = full_pseudopoints),
  E = full_pseudodata_exp
)

# Summarize the posterior marginals of the parameters.
print(full_result$summary.fixed)
print(full_result$summary.hyperpar)

# Plot surface.
plot(im(t(inla.mesh.project(R_proj, full_result$summary.fixed$mean + full_result$summary.random$idx$mean)),
        xrange = Frame(sim_R)$x,
        yrange = Frame(sim_R)$y),
        main = 'Posterior Mean of Log-Intensity')
plot(sim_R, border = '#808080ff', add = TRUE)
points(sim_dataset[[idx]], pch = '.', col = 'white')

# Spatially-averaged random effect variance.
#full_avgvar <- (result_full$summary.random$idx$sd^2) %*% mesh_area / sum(mesh_area)
#print(full_avgvar)

# Max random effect variance at any node. (Interpolated variance cannot exceed this.)
#full_maxvar <- max(result_full$summary.random$idx$sd[mesh_area > 0]^2)
#print(full_maxvar)


#}}}######################################################
#{{{ Sampling parameters applicable to all survey plans. #
##########################################################

# Define the region where waypoints are allowed
sampleable <- erosion(sim_R, WAYPOINT_MARGIN)
min_x <- min(vertices(sampleable)$x)
max_x <- max(vertices(sampleable)$x)
#min_y <- min(vertices(sampleable)$y)
#max_y <- max(vertices(sampleable)$y)
min_y <- min(vertices(sim_R)$y) - XSECT_RADIUS
max_y <- max(vertices(sim_R)$y) + XSECT_RADIUS


#}}}#######
#{{{ SRS. #
###########

srs <- function(full_ppp, num_xsects = XSECT_NUM_INITIAL){
  srs_x <- sort(runif(num_xsects, min_x, max_x))
  waypoints <- cbind(x = rep(srs_x, each = 2), y = c(min_y, max_y, max_y, min_y))
  n_waypoints <- nrow(waypoints)

  cog_psp <- psp(
      x0 = waypoints[-n_waypoints, 'x'],
      y0 = waypoints[-n_waypoints, 'y'],
      x1 = waypoints[-1, 'x'],
      y1 = waypoints[-1, 'y'],
      window = dilation(Frame(sim_R), XSECT_RADIUS)
    )[sim_R]
    srs_D <- dilation(cog_psp, XSECT_RADIUS)
    srs_ppp <- full_ppp[srs_D]
    return(structure(srs_ppp, cog = cog_psp))
}


#}}}############################
#{{{ Systematic random sample. #
################################

sys <- function(full_ppp, num_xsects = XSECT_NUM_INITIAL){
  spacing <- (max_x - min_x) / num_xsects
  edge2edge <- spacing - XSECT_WIDTH
  sys_x <- runif(1, min_x, min_x + edge2edge) + (0:(num_xsects - 1)) * spacing
  waypoints <- cbind(x = rep(sys_x, each = 2), y = c(min_y, max_y, max_y, min_y))
  n_waypoints <- nrow(waypoints)

  cog_psp <- psp(
      x0 = waypoints[-n_waypoints, 'x'],
      y0 = waypoints[-n_waypoints, 'y'],
      x1 = waypoints[-1, 'x'],
      y1 = waypoints[-1, 'y'],
      window = dilation(Frame(sim_R), XSECT_RADIUS)
    )[sim_R]
    sys_D <- dilation(cog_psp, XSECT_RADIUS)
    sys_ppp <- full_ppp[sys_D]
    return(structure(sys_ppp, cog = cog_psp))
}


#}}}############################


stopCluster(cl)

# vim: foldmethod=marker:
