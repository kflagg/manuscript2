#########################
#{{{ Setup environment. #
#########################

# Load packages.
library(parallel)
library(spatstat)
library(lhs)
library(TSP)
library(HilbertVis)
library(INLA)
library(tibble)
library(dplyr)
library(tidyr)

# Create a cluster for applying in parallel.
cl <- makeCluster(ceiling(0.5) * detectCores(), outfile = '')
invisible(clusterEvalQ(cl, {
  library(spatstat)
  library(INLA)
}))

# Define user-specified parameters.
N_SIMS <- 100
N_REALIZED <- 5
XSECT_WIDTH <- 4 # Width of transects.
XSECT_LENGTH_MIN <- 50
XSECT_LENGTH_MAX <- 500
XSECT_LENGTH_CORR <- -0.8
XSECT_SEP_MIN <- 10
XSECT_SEP_MAX <- 250
DIST_INITIAL <- 10000
DIST_MAX <- 50000
WP_NUM_INITIAL <- 400
SERPS_INITIAL <- 5
HILBERT_ORDER_INITIAL <- 3
NUM_PAIRS <- 3
PAIR_RADIUS <- 20 * XSECT_WIDTH
WP_MARGIN <- XSECT_WIDTH / 2 # Minimum distance between waypoints and site boundary.

# Mesh parameters to experiment with.
MAX_EDGE_LENGTH <- 100
MAX_EDGE_EXT <- 200

# Graphics parameters.
NPIX_X <- 250 # Should be 500 but rLGCP to has trouble when NPIX_X != NPIX_Y.
NPIX_Y <- 250


# Define derived parameters
XSECT_RADIUS <- XSECT_WIDTH / 2 # Half the transect width
XSECT_NUM_INITIAL <- floor(DIST_INITIAL / XSECT_LENGTH_MAX)


#}}}##################################
#{{{ Objects pertaining to the site. #
######################################

# Rectangular window.

rect_R <- owin(c(0, 1500), c(0, 700))

# Mesh covering the site.
rect_R_boundary <- inla.mesh.segment(loc = do.call(cbind, vertices.owin(rect_R)))
rect_R_mesh <- inla.mesh.create(
  boundary = rect_R_boundary,
  refine = list(max.edge = MAX_EDGE_LENGTH)
)

# Get the mesh nodes.
rect_R_mesh_loc <- rect_R_mesh$loc[,1:2]

# Convert mesh triangles to owins and create a tesselation.
clusterExport(cl, c('rect_R_mesh_loc', 'rect_R'))
rect_R_mesh_tess <- as.tess(parApply(cl, rect_R_mesh$graph$tv, 1, function(x){
  return(owin(poly = rect_R_mesh_loc[x,]))
}))

# Convert mesh to linear network for plotting.
rect_R_mesh_net <- linnet(as.ppp(rect_R_mesh_loc, rect_R), as.matrix(rect_R_mesh$graph$vv == 1))

# Get the interior area of each triangle.
rect_R_tri_areas <- parSapply(cl, parLapply(cl, tiles(rect_R_mesh_tess), '[', rect_R), area)

# Construct the dual of the mesh.
# Start by getting a list of triangles for each node.
clusterExport(cl, 'rect_R_mesh')
rect_node_tris <- parLapply(cl, parLapply(cl, parLapply(cl, parLapply(cl,
  # Indicate which triangles have this node as a vertex.
  seq_len(rect_R_mesh$n), `==`, rect_R_mesh$graph$tv),
  # Flatten the indicator matrix.
  rowSums),
  # Make the indicator vector logical instead of integer 0/1.
  as.logical),
  # Find which rows of the matrix (triangles) had this node.
  which)

# Get the centroids of each triangle.
rect_tri_centroids <- t(parSapply(cl, tiles(rect_R_mesh_tess), centroid.owin))

# Now, construct the dual mesh. Get the quadrilaterals made by cutting each
# triangle from its edge midpoints to its centroid, indexed by the mesh nodes.
# Then union the quadrilaterals for each node.
clusterExport(cl, c('rect_node_tris', 'rect_tri_centroids'))
rect_dual_tess <- as.tess(parLapply(cl, seq_len(rect_R_mesh$n), function(i){return(
  do.call(union.owin, lapply(rect_node_tris[[i]], function(j){return(
    convexhull.xy(rbind(
      # jth triangle centroid.
      rect_tri_centroids[j,],
      # Average of the ith node and all vertices of the jth triangle.
      # This results in the ith node and the midpoints of its edges.
      t((rect_R_mesh_loc[i,] + t(rect_R_mesh_loc[
        rect_R_mesh$graph$tv[j,],
      1:2])) / 2)
    ))
  )}))
)}))

# Calculate the interior area represented by each node.
rect_R_nodes_area <- parSapply(cl, parLapply(cl, tiles(rect_dual_tess), `[`, rect_R), area)

# Neat plot.
#plot(rect_dual_tess, border = '#80808020', do.col = TRUE, values = rect_R_nodes_area,
#     main = 'Mesh with dual colored by node weight')
#plot(rect_R_mesh_net, add = TRUE, col = '#00000080')
#plot(rect_R, border = 'white', add = TRUE)
#points(rect_R_mesh_loc[,], pch = 20)

# Prediction projector.
rect_R_proj <- inla.mesh.projector(rect_R_mesh, dims = c(NPIX_X, NPIX_Y))


# Complicated window.

#sim_R <- owin(poly = cbind(
#    x = c(   0, 1400, 1400, 2000, 2000,  500,    0),
#    y = c(   0,    0,  600,  600, 2000, 2000, 1500)
#  ), unitnames = c('meter', 'meters'))

# Mesh covering the site.
#sim_R_boundary <- inla.mesh.segment(loc = do.call(cbind, vertices.owin(sim_R)))
#sim_R_mesh <- inla.mesh.create(
#  boundary = sim_R_boundary,
#  refine = list(max.edge = MAX_EDGE_LENGTH)
#)

# Get the mesh nodes.
#sim_R_mesh_loc <- sim_R_mesh$loc[,1:2]

# Convert mesh triangles to owins and create a tesselation.
#clusterExport(cl, c('sim_R_mesh_loc', 'sim_R'))
#sim_R_mesh_tess <- as.tess(parApply(cl, sim_R_mesh$graph$tv, 1, function(x){
#  return(owin(poly = sim_R_mesh_loc[x,]))
#}))

# Convert mesh to linear network for plotting.
#sim_R_mesh_net <- linnet(as.ppp(sim_R_mesh_loc, sim_R), as.matrix(sim_R_mesh$graph$vv == 1))

# Get the interior area of each triangle.
#sim_R_tri_areas <- parSapply(cl, parLapply(cl, tiles(sim_R_mesh_tess), '[', sim_R), area)

# Construct the dual of the mesh.
# Start by getting a list of triangles for each node.
#clusterExport(cl, 'sim_R_mesh')
#sim_node_tris <- parLapply(cl, parLapply(cl, parLapply(cl, parLapply(cl,
#  # Indicate which triangles have this node as a vertex.
#  seq_len(sim_R_mesh$n), `==`, sim_R_mesh$graph$tv),
#  # Flatten the indicator matrix.
#  rowSums),
#  # Make the indicator vector logical instead of integer 0/1.
#  as.logical),
#  # Find which rows of the matrix (triangles) had this node.
#  which)

# Get the centroids of each triangle.
#sim_tri_centroids <- t(parSapply(cl, tiles(sim_R_mesh_tess), centroid.owin))

# Now, construct the dual mesh. Get the quadrilaterals made by cutting each
# triangle from its edge midpoints to its centroid, indexed by the mesh nodes.
# Then union the quadrilaterals for each node.
#clusterExport(cl, c('sim_node_tris', 'sim_tri_centroids'))
#sim_dual_tess <- as.tess(parLapply(cl, seq_len(sim_R_mesh$n), function(i){return(
#  do.call(union.owin, lapply(sim_node_tris[[i]], function(j){return(
#    convexhull.xy(rbind(
#      # jth triangle centroid.
#      sim_tri_centroids[j,],
#      # Average of the ith node and all vertices of the jth triangle.
#      # This results in the ith node and the midpoints of its edges.
#      t((sim_R_mesh_loc[i,] + t(sim_R_mesh_loc[
#        sim_R_mesh$graph$tv[j,],
#      1:2])) / 2)
#    ))
#  )}))
#)}))

# Calculate the interior area represented by each node.
#sim_R_nodes_area <- parSapply(cl, parLapply(cl, tiles(sim_dual_tess), `[`, sim_R), area)

# Neat plot.
#plot(sim_dual_tess, border = '#80808020', do.col = TRUE, values = sim_R_nodes_area,
#     main = 'Mesh with dual colored by node weight')
#plot(sim_R_mesh_net, add = TRUE, col = '#00000080')
#plot(sim_R, border = 'white', add = TRUE)
#points(sim_R_mesh_loc[,], pch = 20)

# Prediction projector.
#sim_R_proj <- inla.mesh.projector(sim_R_mesh, dims = c(NPIX_X, NPIX_Y))


#}}}##################################
#{{{ Objects pertaining to the data. #
######################################

# Function to simulate one realized point pattern. This includes an LGCP
# background process with Matern (alpha = 2) covariance and a foreground
# cluster process. The true intensity function is stored in an attribute
# called Lambda.
# Dots should include dimyx = c(NPIX_Y, NPIX_X).
simulate_data <- function(
  R = rect_R,
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

lgcp_design <- tibble(
  Model = 'LGCP',
  matern_mu = log(250 / area(rect_R)),
  matern_sd = 2,
  matern_range = 200,
  thomas_kappa = 0,
  thomas_scale = 50,
  thomas_mu = 0,
  DataID = sprintf('LGCP%06d', seq_len(N_REALIZED))
)

clust_design <- tibble(
  Model = 'Cluster',
  matern_mu = log(250 / area(rect_R)),
  matern_sd = 1,
  matern_range = 200,
  thomas_kappa = 3 / area(rect_R),
  thomas_scale = 50,
  thomas_mu = 200,
  DataID = sprintf('Cluster%06d', seq_len(N_REALIZED))
)

rect_design <- bind_rows(lgcp_design, clust_design)

invisible(clusterEvalQ(cl, {library(tibble);library(dplyr)}))
clusterExport(cl, c('simulate_data', 'rect_design', 'NPIX_X', 'NPIX_Y'))
rect_datasets <- bind_rows(
  parLapply(cl, seq_len(nrow(rect_design)), function(r){
    message(rect_design$DataID[r])
    dataset <- simulate_data(
      rect_R, rect_design$matern_mu[r], rect_design$matern_sd[r],
      rect_design$matern_range[r], rect_design$thomas_kappa[r],
      rect_design$thomas_scale[r], rect_design$thomas_mu[r],
      dimyx = c(NPIX_Y, NPIX_X))
    return(tibble_row(
      DataID = rect_design$DataID[r],
      Data = list(dataset)
    ))
  })
)


#}}}###################################################################
#{{{ Objects pertaining to the model but independent of survey plans. #
#######################################################################

# Define an SPDE representation of the spatial GP using PC priors.
rect_R_spde <- inla.spde2.pcmatern(
  mesh = rect_R_mesh,
  alpha = 2,
  prior.range = c(100, 0.1),
  prior.sigma = c(3, 0.1)
)

# Prior means and precisions for coefficients.
rect_prior_fixed = list(
    mean.intercept = 0,
    prec.intercept = 0,
    mean = 0,
    prec = 0
)

# Model formula with no covariates.
rect_R_formula <- y ~ -1 + intercept + f(idx, model = rect_R_spde)


#}}}######################################################
#{{{ Sampling parameters applicable to all survey plans. #
##########################################################

# Function to subset a ppp to a region along a psp.
sample_ppp <- function(full_ppp, path, xsect_radius = XSECT_RADIUS){
  if(!is.psp(path)){
    path <- as.psp(path)
  }
  obs_D <- dilation(path, xsect_radius)
  obs_ppp <- full_ppp[obs_D]
  return(structure(obs_ppp, path = path, Lambda = attr(full_ppp, 'Lambda')))
}

# Function to compute the total length of a path.
totallength <- function(x){
  if(inherits(x, 'linnet')){
    return(volume(x))
  }
  return(sum(lengths.psp(as.psp(x))))
}

# Lengths of segments in a linear network.
lengths.linnet <- function(x){
  return(lengths.psp(as.psp(x)))
}

# Number of corners in a linear network.
cornercount.linnet <- function(x){
  return(sum(x$m) - x$vertices$n)
}

# Angles of corners in a linear network.
angles.linnet <- function(x){
  adj <- x$m + diag(x$vertices$n)
  links <- apply(adj > 0, 1, which)
  links <- links[lengths(links) == 3]
  verts <- cbind(x$vertices$x, x$vertices$y)
  return(sapply(links, function(idx){
    vs <- verts[idx,]
    v1 <- vs[2,] - vs[1,]
    v3 <- vs[3,] - vs[2,]
    a <- atan2(v3[1], v3[2]) - atan2(v1[1], v1[2])
    return(if(a > pi){a - 2 * pi}else if(a < -pi){a + 2 * pi}else{a})
  }))
}

# Nearest neighbor distance on a linear network.
nndist.linnet <- function(X, k = 1, agg = 'avg', ...){
  numsegs <- X$lines$n

  if(k < 1){
    return(0)
  }

  # Find which other segments each segment shares vertices with.
  segshare <- lapply(seq_len(numsegs), function(s){
    return(sort(unique(c(
      which(X$from == X$from[s]),
      which(X$to == X$from[s]),
      which(X$from == X$to[s]),
      which(X$to == X$to[s])
    ))))
  })
  segfrom <- unlist(lapply(seq_len(numsegs), function(s){return(rep(s, length(segshare[[s]])))}))
  segto <- unlist(segshare)

  # Construct the adjacency matrix for segments.
  adj <- sparseMatrix(i = segfrom, j = segto, symmetric = TRUE)

  # Construct the k-1 order adjacency matrix.
  adj_k <- sparseMatrix(i = seq_len(numsegs), j = seq_len(numsegs), symmetric = TRUE)
  if(k > 1){
    for(power_counter in 2:k){
      adj_k <- adj_k %*% adj
    }
  }

  # Get the kth order inclusion matrix.
  incl <- !adj_k

  rawdist <- sapply(seq_len(numsegs), function(s){
    return(min(crossdist.psp(X$lines[s], X$lines[incl[s,]], type = 'separation')))
  })

  # Return the average if that was requested.
  if(agg == 'avg'){
    seglengths <- lengths.linnet(X)
    return(sum(rawdist * seglengths) / sum(seglengths))
  }

  # Otherwise return the minimum.
  return(min(rawdist))
}

clusterExport(cl, c('totallength', 'lengths.linnet', 'angles.linnet', 'cornercount.linnet', 'nndist.linnet'))


#}}}#######
#{{{ SRS. #
###########

srs <- function(full_win, num_xsects = XSECT_NUM_INITIAL, xsect_radius = XSECT_RADIUS){
  full_frame <- Frame(full_win)
  min_x <- min(full_frame$x) + xsect_radius
  max_x <- max(full_frame$x) - xsect_radius
  min_y <- min(full_frame$y)
  max_y <- max(full_frame$y)

  srs_x <- sort(runif(num_xsects, min_x, max_x))
  waypoints <- cbind(
    x = rep(srs_x, each = 2),
    y = rep(c(min_y, max_y, max_y, min_y), ceiling(num_xsects / 2))[
      1:(2 * num_xsects)
    ]
  )
  n_waypoints <- nrow(waypoints)
  n_segs <- n_waypoints / 2

  path_linnet <- linnet(
    vertices = as.ppp(waypoints, W = full_frame),
    edges = cbind(seq_len(n_segs) * 2 - 1, seq_len(n_segs) * 2),
    sparse = TRUE,
    warn = FALSE
  )
  return(path_linnet)[full_win]
}

srs_design <- tibble(
  num_xsects = rep(c(10, 25, 50, 70), each = N_SIMS),
  xsect_radius = XSECT_RADIUS
) %>% mutate(PlanID = sprintf('SRS%06d', 1:n()))
clusterExport(cl, c('srs', 'srs_design'))
srs_plans <- bind_rows(
  parLapply(cl, seq_len(nrow(srs_design)), function(r){
    plan <- srs(rect_R, srs_design$num_xsects[r], srs_design$xsect_radius[r])
    message(srs_design$PlanID[r])
    return(tibble_row(
      PlanID = srs_design$PlanID[r],
      Plan = list(plan),
      xsect_radius = srs_design$xsect_radius[r],
      Distance = totallength(plan),
      Lengths = list(lengths.linnet(plan)),
      Angles = list(angles.linnet(plan)),
      Corners = cornercount.linnet(plan),
      MinNND_2 = nndist(plan, 2, 'min'),
      AvgNND_2 = nndist(plan, 2, 'avg'),
      AvgNND_1 = nndist(plan, 1, 'avg')
   ))})
)


#}}}############################
#{{{ Systematic random sample. #
################################

sys <- function(full_win, num_xsects = XSECT_NUM_INITIAL, xsect_radius = XSECT_RADIUS){
  full_frame <- Frame(full_win)
  min_x <- min(full_frame$x) + xsect_radius
  max_x <- max(full_frame$x) - xsect_radius
  min_y <- min(full_frame$y)
  max_y <- max(full_frame$y)

  spacing <- (max_x - min_x) / num_xsects
  edge2edge <- max(spacing - xsect_radius * 2, 0)
  sys_x <- runif(1, min_x, min_x + edge2edge) + (0:(num_xsects - 1)) * spacing
  waypoints <- cbind(
    x = rep(sys_x, each = 2),
    y = rep(c(min_y, max_y, max_y, min_y), ceiling(num_xsects / 2))[
      1:(2 * num_xsects)
    ]
  )
  n_waypoints <- nrow(waypoints)
  n_segs <- n_waypoints / 2

  path_linnet <- linnet(
    vertices = as.ppp(waypoints, W = full_frame),
    edges = cbind(seq_len(n_segs) * 2 - 1, seq_len(n_segs) * 2),
    sparse = TRUE,
    warn = FALSE
  )
  return(path_linnet)[full_win]
}

sys_design <- tibble(
  num_xsects = rep(c(10, 25, 50, 70), each = N_SIMS),
  xsect_radius = XSECT_RADIUS
) %>% mutate(PlanID = sprintf('Sys%06d', 1:n()))
clusterExport(cl, c('sys', 'sys_design'))
sys_plans <- bind_rows(
  parLapply(cl, seq_len(nrow(sys_design)), function(r){
    plan <- sys(rect_R, sys_design$num_xsects[r], sys_design$xsect_radius[r])
    message(sys_design$PlanID[r])
    return(tibble_row(
      PlanID = sys_design$PlanID[r],
      Plan = list(plan),
      xsect_radius = sys_design$xsect_radius[r],
      Distance = totallength(plan),
      Lengths = list(lengths.linnet(plan)),
      Angles = list(angles.linnet(plan)),
      Corners = cornercount.linnet(plan),
      MinNND_2 = nndist(plan, 2, 'min'),
      AvgNND_2 = nndist(plan, 2, 'avg'),
      AvgNND_1 = nndist(plan, 1, 'avg')
   ))})
)


#}}}########################
#{{{ Serpentine transects. #
############################

serp <- function(full_win, num_xsects = XSECT_NUM_INITIAL,
                 serp_dist = XSECT_SEP_MIN, serp_num = SERPS_INITIAL,
                 xsect_radius = XSECT_RADIUS){
  full_frame <- Frame(full_win)
  min_x <- min(full_frame$x) + xsect_radius - ifelse(serp_dist < 0, serp_dist, 0)
  max_x <- max(full_frame$x) - xsect_radius - ifelse(serp_dist > 0, serp_dist, 0)
  min_y <- min(full_frame$y)
  max_y <- max(full_frame$y)

  spacing <- (max_x - min_x) / num_xsects
  edge2edge <- max(spacing - xsect_radius * 2, 0)
  waypoints <- cbind(
    x <- runif(1, min_x, min_x + edge2edge) +
      rep(0:(num_xsects - 1), each = serp_num * 2) * spacing +
      rep(
        rep(c(0, 0, 1, 1), serp_num)[c(1:(2 * serp_num), (2 * serp_num):1)],
        ceiling(num_xsects / 2)
      )[seq_len(num_xsects * serp_num * 2)] * serp_dist,
    y = rep(rep(seq(min_y, max_y, length.out = serp_num + 1), each = 2)[
      c(2:(2 * serp_num + 1), (2 * serp_num + 1):2)
    ], ceiling(num_xsects / 2))[seq_len(num_xsects * serp_num * 2)]
  )
  n_waypoints <- nrow(waypoints)
  n_segs <- n_waypoints - 1

  path_linnet <- linnet(
    vertices = as.ppp(waypoints, W = full_frame, xsect_radius),
    edges = cbind(1:n_segs, 2:n_waypoints)[-seq(0, n_waypoints, 2 * serp_num),],
    sparse = TRUE,
    warn = FALSE
  )
  return(path_linnet)[full_win]
}

serp_design <- tibble(
    reps = N_SIMS,
    expand.grid(
      num_xsects = c(7, 22, 47, 67),
      serp_num = c(5, 8)
    ),
    xsect_radius = XSECT_RADIUS
  ) %>%
  mutate(serp_dist = 2100 / num_xsects / (serp_num - 1)) %>%
  uncount(reps) %>%
  mutate(PlanID = sprintf('Serp%06d', 1:n()))
clusterExport(cl, c('serp', 'serp_design'))
serp_plans <- bind_rows(
  parLapply(cl, seq_len(nrow(serp_design)), function(r){
    plan <- serp(rect_R, serp_design$num_xsects[r], serp_design$serp_dist[r],
                 serp_num = serp_design$serp_num[r], serp_design$xsect_radius[r])
    message(serp_design$PlanID[r])
    return(tibble_row(
      PlanID = serp_design$PlanID[r],
      Plan = list(plan),
      xsect_radius = serp_design$xsect_radius[r],
      Distance = totallength(plan),
      Lengths = list(lengths.linnet(plan)),
      Angles = list(angles.linnet(plan)),
      Corners = cornercount.linnet(plan),
      MinNND_2 = nndist(plan, 2, 'min'),
      AvgNND_2 = nndist(plan, 2, 'avg'),
      AvgNND_1 = nndist(plan, 1, 'avg')
   ))})
)


#}}}###############################
#{{{ Inhibitory plus close pairs. #
###################################

inhib <- function(full_win, num_primary = XSECT_NUM_INITIAL - num_paired,
                  num_paired = NUM_PAIRS, pair_radius = PAIR_RADIUS,
                  antirepulsion = 0.05, xsect_radius = XSECT_RADIUS,
                  mc_iter = 1000, mc_radius = pair_radius * 2){
  full_frame <- Frame(full_win)
  min_x <- min(full_frame$x) + xsect_radius
  max_x <- max(full_frame$x) - xsect_radius
  min_y <- min(full_frame$y)
  max_y <- max(full_frame$y)
  num_xsects <- num_primary + num_paired

  initial_x <- runif(num_primary, min_x, max_x)

  # A simple Metropolis-Hastings algorithm.
  for(iter in seq_len(mc_iter)){
    # Count the number of current close pairs.
    current_pairs <- 0L
    for(i in 1:(num_primary - 1)){
      current_pairs <- current_pairs +
        sum(abs(initial_x[i] - initial_x[-(1:i)]) < pair_radius)
    }

    # Choose one to perturb.
    perturb_idx <- sample.int(num_primary, 1)

    # Perturb it.
    prop_x <- initial_x
    prop_x[perturb_idx] <- runif(1,
       max(prop_x[perturb_idx] - mc_radius, min_x),
       min(prop_x[perturb_idx] + mc_radius, max_x)
    )

    # Count the number of proposed close pairs.
    prop_pairs <- 0L
    for(i in 1:(num_primary - 1)){
      prop_pairs <- prop_pairs +
        sum(abs(prop_x[i] - prop_x[-(1:i)]) < pair_radius)
    }

    # Accept or reject the proposal.
    if(runif(1) < antirepulsion^(prop_pairs - current_pairs)){
      initial_x <- prop_x
    }
  }

  paired_x <- sample(initial_x, num_paired)
  pairs_x <- runif(num_paired,
                   pmax(paired_x - pair_radius, min_x),
                   pmin(paired_x + pair_radius, max_x))
  inhib_x <- sort(c(initial_x, pairs_x))
  waypoints <- cbind(
    x = rep(inhib_x, each = 2),
    y = rep(c(min_y, max_y, max_y, min_y), ceiling(num_xsects / 2))[
      1:(2 * num_xsects)
    ]
  )
  n_waypoints <- nrow(waypoints)
  n_segs <- n_waypoints / 2

  path_linnet <- linnet(
    vertices = as.ppp(waypoints, W = full_frame),
    edges = cbind(seq_len(n_segs) * 2 - 1, seq_len(n_segs) * 2),
    sparse = TRUE,
    warn = FALSE
  )
  return(path_linnet)[full_win]
}

inhib_design <- tibble(
    reps = N_SIMS,
    expand.grid(
      num_xsects = c(10, 25, 50, 70),
      prop_pairs = c(0.1, 0.2)
    ),
    xsect_radius = XSECT_RADIUS
  ) %>%
  mutate(
    num_pair = round(num_xsects * prop_pairs),
    num_prim = num_xsects - num_pair,
    pair_radius = 1500 / num_xsects
  ) %>%
  uncount(reps) %>%
  mutate(PlanID = sprintf('Inhib%06d', 1:n()))
clusterExport(cl, c('inhib', 'inhib_design'))
inhib_plans <- bind_rows(
  parLapply(cl, seq_len(nrow(inhib_design)), function(r){
    plan <- inhib(rect_R, inhib_design$num_prim[r], inhib_design$num_pair[r], inhib_design$pair_radius[r], xsect_radius = inhib_design$xsect_radius[r])
    message(inhib_design$PlanID[r])
    return(tibble_row(
      PlanID = inhib_design$PlanID[r],
      Plan = list(plan),
      xsect_radius = inhib_design$xsect_radius[r],
      Distance = totallength(plan),
      Lengths = list(lengths.linnet(plan)),
      Angles = list(angles.linnet(plan)),
      Corners = cornercount.linnet(plan),
      MinNND_2 = nndist(plan, 2, 'min'),
      AvgNND_2 = nndist(plan, 2, 'avg'),
      AvgNND_1 = nndist(plan, 1, 'avg')
   ))})
)


#}}}#########################
#{{{ Latin traveling sales. #
#############################

lhstsp <- function(full_win, num_bins = round(sqrt(WP_NUM_INITIAL)), margin = WP_MARGIN, ...){
  full_frame <- Frame(full_win)
  min_x <- min(full_frame$x) + margin
  max_x <- max(full_frame$x) - margin
  min_y <- min(full_frame$y) + margin
  max_y <- max(full_frame$y) - margin
  sampleable <- erosion(full_win, margin)

  wp_unscaled <- maximinLHS(n = num_bins, k = 2, ...)
  wp_x <- wp_unscaled[,1] * (max_x - min_x) + min_x
  wp_y <- wp_unscaled[,2] * (max_y - min_y) + min_y
  wp_ppp <- ppp(wp_x, wp_y, window = full_frame)[sampleable]
  waypoints <- as.data.frame(wp_ppp)
  n_waypoints <- nrow(waypoints)
  n_segs <- n_waypoints - 1

  wp_tsp <- ETSP(waypoints)
  wp_tour <- waypoints[solve_TSP(wp_tsp, ...),]

  path_linnet <- linnet(
    vertices = as.ppp(wp_tour, W = full_frame),
    edges = cbind(seq_len(n_waypoints), c(n_waypoints, seq_len(n_segs))),
    sparse = TRUE,
    warn = FALSE
  )
  return(path_linnet)[full_win]
}

lhs_design <- tibble(
    bins = rep(c(50, 300, 1200, 2400), each = N_SIMS),
    xsect_radius = XSECT_RADIUS
  ) %>%
  mutate(PlanID = sprintf('LHS-TSP%06d', 1:n()))
invisible(clusterEvalQ(cl, {library(lhs);library(TSP)}))
clusterExport(cl, c('lhstsp', 'lhs_design'))
lhs_plans <- bind_rows(
  parLapply(cl, seq_len(nrow(lhs_design)), function(r){
    plan <- lhstsp(rect_R, lhs_design$bins[r], lhs_design$xsect_radius[r])
    message(lhs_design$PlanID[r])
    return(tibble_row(
      PlanID = lhs_design$PlanID[r],
      Plan = list(plan),
      xsect_radius = lhs_design$xsect_radius[r],
      Distance = totallength(plan),
      Lengths = list(lengths.linnet(plan)),
      Angles = list(angles.linnet(plan)),
      Corners = cornercount.linnet(plan),
      MinNND_2 = nndist(plan, 2, 'min'),
      AvgNND_2 = nndist(plan, 2, 'avg'),
      AvgNND_1 = nndist(plan, 1, 'avg')
   ))})
)


#}}}#################
#{{{ Hilbert curve. #
#####################

hilbert <- function(full_win, h_order = HILBERT_ORDER_INITIAL, margin = WP_MARGIN, ...){
  full_frame <- Frame(full_win)
  min_x <- min(full_frame$x) + margin
  max_x <- max(full_frame$x) - margin
  min_y <- min(full_frame$y) + margin
  max_y <- max(full_frame$y) - margin
  sampleable <- erosion(full_win, margin)

  wp_unscaled <- hilbertCurve(h_order)
  wp_x <- (wp_unscaled[,1] + runif(1)) / (2^h_order) * (max_x - min_x) + min_x
  wp_y <- (wp_unscaled[,2] + runif(1)) / (2^h_order) * (max_y - min_y) + min_y
  wp_ppp <- ppp(wp_x, wp_y, window = full_frame)#[sampleable]
  waypoints <- as.data.frame(wp_ppp)
  n_waypoints <- nrow(waypoints)
  n_segs <- n_waypoints - 1

  path_linnet <- linnet(
    vertices = as.ppp(waypoints, W = full_frame),
    edges = cbind(seq_len(n_segs), 2:n_waypoints),
    sparse = TRUE,
    warn = FALSE
  )
  return(path_linnet)[full_win]
}

hilb_design <- tibble(
  order = rep(c(3, 4, 5, 6), each = N_SIMS),
  xsect_radius = XSECT_RADIUS
) %>% mutate(PlanID = sprintf('Hilbert%06d', 1:n()))
invisible(clusterEvalQ(cl, library(HilbertVis)))
clusterExport(cl, c('hilbert', 'hilb_design'))
hilb_plans <- bind_rows(
  parLapply(cl, seq_len(nrow(hilb_design)), function(r){
    message(hilb_design$PlanID[r])
    plan <- hilbert(rect_R, hilb_design$order[r], hilb_design$xsect_radius[r])
    return(tibble_row(
      PlanID = hilb_design$PlanID[r],
      Plan = list(plan),
      xsect_radius = hilb_design$xsect_radius[r],
      Distance = totallength(plan),
      Lengths = list(lengths.linnet(plan)),
      Angles = list(angles.linnet(plan)),
      Corners = cornercount.linnet(plan),
      MinNND_2 = nndist(plan, 2, 'min'),
      AvgNND_2 = nndist(plan, 2, 'avg'),
      AvgNND_1 = nndist(plan, 1, 'avg')
   ))})
)


#}}}############################
#{{{ Random particle movement. #
################################

rpm <- function(full_win, dist_cutoff = DIST_MAX, corr = XSECT_LENGTH_CORR,
                seg_min = XSECT_LENGTH_MIN, seg_max = XSECT_LENGTH_MAX,
                angle_m = pi/3, angle_s = pi/6, angle_prob = c(0.5, 0.5),
                a = 1, b = 1, pair_radius = PAIR_RADIUS, antirepulsion = 0.8,
                margin = WP_MARGIN, animate = FALSE, ...){
  if(corr < 0 & (a != 1 | b != 1)){
    warning('Negative correlation only available for a = b = 1!')
  }
  full_frame <- Frame(full_win)
  min_x <- min(full_frame$x) + margin
  max_x <- max(full_frame$x) - margin
  min_y <- min(full_frame$y) + margin
  max_y <- max(full_frame$y) - margin
  sampleable <- erosion(full_win, margin)

  dist_range <- seg_max - seg_min

  pos_corr <- corr > 0
  if(!pos_corr){
    corr <- -corr
  }
  p <- corr * (a + b) / (corr + b)
  if(p >= a){
    p <- a
    warning('corr is to big for the given a!')
  }

  # Starting segment.
  wp <- as.data.frame(runifpoint(1, sampleable))

  reject <- TRUE
  while(reject){
    angle <- runif(1, 0, 2*pi)
    dist_unscaled <- rbeta(1, a, b)
    new_dist <- dist_unscaled * dist_range + seg_min

    new_x <- tail(wp$x, 1) + cos(angle) * new_dist
    new_y <- tail(wp$y, 1) + sin(angle) * new_dist

    if(inside.owin(new_x, new_y, sampleable)){
        reject <- FALSE
    }
  }
  path_psp <- psp(wp$x, wp$y, new_x, new_y, rect_R)
  cum_dist <- new_dist
  current_x <- new_x
  current_y <- new_y

  if(animate){
    plot(path_psp, ...)
  }

  # Loop until exceeding the distance cutoff.
  while(cum_dist < dist_cutoff){
    reject <- TRUE

    while(reject){
      new_angle <- (angle + sample(c(-1, 1), 1, prob = angle_prob) *
        rnorm(1, angle_m, angle_s)) %% (2*pi)
      new_dist_unscaled <- rbeta(1, b, a - p) *
        (1 - rbeta(1, p, a - p) * dist_unscaled)
      if(pos_corr){
        new_dist_unscaled <- 1 - new_dist_unscaled
      }

      dist_unscaled <- new_dist_unscaled
      new_dist <- dist_unscaled * dist_range + seg_min

      new_x <- current_x + cos(new_angle) * new_dist
      new_y <- current_y + sin(new_angle) * new_dist

      if(inside.owin(new_x, new_y, sampleable)){
        if(runif(1) < antirepulsion^(sum(crossdist(
            psp(current_x, current_y, new_x, new_y, rect_R),
            path_psp, type = 'separation') < pair_radius) - 1)){
          reject <- FALSE
        }
      }
    }
    angle <- new_angle
    cum_dist <- cum_dist + new_dist
    wp <- rbind(wp, data.frame(x = new_x, y = new_y))
    path_psp$ends <- rbind(path_psp$ends, data.frame(
      x0 = current_x,
      y0 = current_y,
      x1 = new_x,
      y1 = new_y
    ))
    path_psp$n <- path_psp$n + 1
    current_x <- new_x
    current_y <- new_y

    if(animate){
      plot(path_psp, add = TRUE)
    }
  }

  path_linnet <- linnet(
    vertices = ppp(c(path_psp$ends$x0, current_x),
                   c(path_psp$ends$y0, current_y),
                   window = full_frame),
    edges = cbind(seq_len(path_psp$n), 2:(path_psp$n + 1)),
    sparse = TRUE,
    warn = FALSE
  )

  return(path_linnet[full_win])
}

rpm_design <- tibble(
    reps = N_SIMS,
    expand.grid(
      max_dist = c(6700, 17200, 34700, 49700),
      length_corr = c(0, -0.8),
      seg_max = 500,
      seg_min_prop = 0.1,
      pair_radius = c(80, 300),
      angle_m = c(pi/3, pi/2)
    ),
    xsect_radius = XSECT_RADIUS
  ) %>%
  filter(pair_radius == 80 | max_dist < 10000) %>%
  mutate(
    seg_min = seg_max * seg_min_prop,
    angle_s = ifelse(angle_m > pi/3, pi/12, pi/6)
  ) %>%
  uncount(reps) %>%
  mutate(PlanID = sprintf('RPM%06d', 1:n()))
clusterExport(cl, c('rpm', 'rpm_design'))
rpm_plans <- bind_rows(
  parLapply(cl, seq_len(nrow(rpm_design)), function(r){
    plan <- rpm(rect_R, rpm_design$max_dist[r], rpm_design$length_corr[r], rpm_design$seg_min[r], rpm_design$seg_max[r], rpm_design$angle_m[r], rpm_design$angle_s[r], pair_radius = rpm_design$pair_radius[r], margin = rpm_design$xsect_radius[r])
    message(rpm_design$PlanID[r])
    return(tibble_row(
      PlanID = rpm_design$PlanID[r],
      Plan = list(plan),
      xsect_radius = rpm_design$xsect_radius[r],
      Distance = totallength(plan),
      Lengths = list(lengths.linnet(plan)),
      Angles = list(angles.linnet(plan)),
      Corners = cornercount.linnet(plan),
      MinNND_2 = nndist(plan, 2, 'min'),
      AvgNND_2 = nndist(plan, 2, 'avg'),
      AvgNND_1 = nndist(plan, 1, 'avg')
   ))})
)


#}}}##############################
#{{{ Fit model to observed data. #
##################################

model_fit <- function(model_formula, obs_ppp, rect_R_mesh, dual_tess, rect_R_proj,
  control.fixed = list(
    # Prior means and precisions for coefficients.
    mean.intercept = 0,
    prec.intercept = 0,
    mean = 0,
    prec = 0
  ),
  ...){

  # Observation window.
  S_win <- Window(obs_ppp)

  # Observed event locations.
  obs_pts <- cbind(obs_ppp$x, obs_ppp$y)

  # Get the numbers of mesh nodes and real events.
  # The sum of these will be the number of pseudodata points.
  mesh_size <- rect_R_mesh$n
  n_events <- nrow(obs_pts)

  # Create the psuedodata. This is a vector giving the count of events at each
  # pseudopoint location, that is 0 at each mesh node and 1 at each event.
  pseudodata <- c(rep(0, mesh_size), rep(1, n_events))

  # Calculate the interior area represented by each node.
  mesh_weights <- sapply(lapply(tiles(dual_tess), `[`, S_win), area)

  # Concatenate the weight vector with a vector of zeros for the observed events.
  # This is the vector of Poisson exposure parameters for the pseudodata.
  pseudodata_exp <- c(mesh_weights, rep(0, n_events))

  # Compute the barycentric coordinates of the observed events
  # (i.e. project into the space spanned by the basis of mesh nodes).
  bary <- inla.mesh.project(rect_R_mesh, obs_pts)$A

  # Compute the barycentric coordinates of the nodes. Because the
  # node coordinatess are the basis vectors, this is an identity matrix.
  int_matrix <- sparseMatrix(
    i = seq_len(mesh_size),
    j = seq_len(mesh_size),
    x = rep(1, mesh_size)
  )

  # Bind the node and event coordinates into a single matrix of pseudodata
  # locations in barycentric coordinates.
  pseudopoints <- rbind(int_matrix, bary)

  # Create the data list to pass to inla().
  # Indices and intercepts are only needed for the nodes.
  inla_data <- list(
    y = pseudodata, # The whole pseudodata vector.
    idx = seq_len(mesh_size), # Indices of the nodes.
    intercept = rep(1, mesh_size) # Intercept column.
  )

  # Fit the model as a Poisson GLM with exposures specified.
  result <- inla(
    formula = model_formula,
    data = inla_data,
    family = 'poisson',
    control.fixed = control.fixed,
    control.predictor = list(A = pseudopoints),
    E = pseudodata_exp,
    ...
  )

  gpmap <- im(
    t(inla.mesh.project(rect_R_proj, result$summary.random[[1]]$mean)) +
      result$summary.fixed['intercept', 'mean'],
    xrange = rect_R$x,
    yrange = rect_R$y
  )

  return(tibble_row(
    Fit = list(result),
    Prediction = list(gpmap),
    MSPE = mean((log(attr(obs_ppp, 'Lambda')) - gpmap)^2),
    APV = sum(mesh_weights * result$summary.random[[1]]$sd^2) / sum(mesh_weights),
    MaxPV = max(result$summary.random[[1]]$sd^2)
  ))
}


#}}}##############################

allplans <- bind_rows(
  tibble(
    Scheme = 'SRS',
    srs_plans
  ),
  tibble(
    Scheme = 'Sys',
    sys_plans
  ),
  tibble(
    Scheme = 'Inhib',
    inhib_plans
  ),
  tibble(
    Scheme = 'Serp',
    serp_plans
  ),
  tibble(
    Scheme = 'LHS-TSP',
    lhs_plans
  ),
  tibble(
    Scheme = 'Hilbert',
    hilb_plans
  ),
  tibble(
    Scheme = 'RPM',
    rpm_plans
  )
)

fit_design <- expand.grid(PlanID = allplans$PlanID, DataID = rect_datasets$DataID)

clusterExport(cl, c('sample_ppp', 'model_fit', 'rect_R_spde', 'rect_R_formula', 'rect_dual_tess', 'rect_R_proj', 'allplans', 'fit_design', 'rect_datasets', 'rect_prior_fixed'))
rect_results <- bind_rows(parLapply(cl, seq_len(nrow(fit_design)), function(r){
  thisplan <- allplans %>% filter(PlanID == fit_design$PlanID[r])
  thisdataset <- rect_datasets %>% filter(DataID == fit_design$DataID[r])
  obs_ppp <- sample_ppp(thisdataset$Data[[1]], thisplan$Plan[[1]], thisplan$xsect_radius)
  message(sprintf('Dataset %s, Plan %s, observed %d points',
                  thisdataset$DataID, thisplan$PlanID, obs_ppp$n))
  return(bind_rows(
    tibble_row(DataID = thisdataset$DataID, PlanID = thisplan$PlanID),
    Fit = model_fit(rect_R_formula, obs_ppp, rect_R_mesh, rect_dual_tess, rect_R_proj, rect_prior_fixed)
  ))}))

saveRDS(rect_datasets, 'rect_data.rds')
saveRDS(allplans, 'rect_plans.rds')
saveRDS(rect_results, 'rect_results.rds')

stopCluster(cl)

# vim: foldmethod=marker:
