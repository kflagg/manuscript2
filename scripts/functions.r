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
cl <- makeCluster(ceiling(0.75 * detectCores()), outfile = '')
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
FINE_EDGE_LENGTH <- 20

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
finemesh <- inla.mesh.create(
  boundary = rect_R_boundary,
  refine = list(max.edge = FINE_EDGE_LENGTH)
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

# Prediction projector.
rect_R_proj <- inla.mesh.projector(rect_R_mesh, dims = c(NPIX_X, NPIX_Y))
fineproj <- inla.mesh.projector(finemesh, dims = c(NPIX_X, NPIX_Y))


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

# Number of segments in a linear network.
segmentcount.linnet <- function(x){
  return(x$lines$n)
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

# Distances from nodes to the path.
pointdist <- function(path, mesh = finemesh, meshproj = fineproj){
  path <- as.psp(path)
  win <- Window(path)
  pts <- mesh$loc[,1:2]
  meshweights <- diag(inla.mesh.fem(mesh)$c0)
  node_dists <- apply(pts, 1, function(xy){
    return(min(crossdist.psp(psp(xy[1], xy[2], xy[1], xy[2], win),
                             path, type = 'separation')))
  })
  interp_dists <- inla.mesh.project(meshproj, node_dists)
  dmap <- im(
    t(interp_dists),
    xrange = win$x,
    yrange = win$y
  )
  return(structure(
    dmap,
    avg = sum(meshweights * node_dists) / sum(meshweights),
    max = max(interp_dists)
  ))
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

clusterExport(cl, c('totallength', 'lengths.linnet', 'segmentcount.linnet', 'angles.linnet', 'cornercount.linnet', 'pointdist', 'nndist.linnet'))


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
    Subscheme = paste0('SRS', rep(seq_len(4), each = N_SIMS)),
    num_xsects = rep(c(10, 25, 50, 70), each = N_SIMS),
    xsect_radius = XSECT_RADIUS
  ) %>% mutate(PlanID = sprintf('SRS%06d', 1:n()))
clusterExport(cl, c('srs', 'srs_design'))


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
    Subscheme = paste0('Sys', rep(seq_len(4), each = N_SIMS)),
    num_xsects = rep(c(10, 25, 50, 70), each = N_SIMS),
    xsect_radius = XSECT_RADIUS
  ) %>% mutate(PlanID = sprintf('Sys%06d', 1:n()))
clusterExport(cl, c('sys', 'sys_design'))


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
    Subscheme = paste0('Serp', seq_len(8)),
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
    Subscheme = paste0('Inhib', seq_len(8)),
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
    Subscheme = paste0('LHS-TSP', rep(seq_len(4), each = N_SIMS)),
    bins = rep(c(50, 300, 1200, 2400), each = N_SIMS),
    xsect_radius = XSECT_RADIUS
  ) %>%
  mutate(PlanID = sprintf('LHS-TSP%06d', 1:n()))
invisible(clusterEvalQ(cl, {library(lhs);library(TSP)}))
clusterExport(cl, c('lhstsp', 'lhs_design'))


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
    Subscheme = paste0('Hilbert', rep(seq_len(4), each = N_SIMS)),
    order = rep(c(3, 4, 5, 6), each = N_SIMS),
    xsect_radius = XSECT_RADIUS
  ) %>% mutate(PlanID = sprintf('Hilbert%06d', 1:n()))
invisible(clusterEvalQ(cl, library(HilbertVis)))
clusterExport(cl, c('hilbert', 'hilb_design'))


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
  save_fit = FALSE, save_pred = TRUE,
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
  result <- tryCatch(
    inla(
      formula = model_formula,
      data = inla_data,
      family = 'poisson',
      control.fixed = control.fixed,
      control.predictor = list(A = pseudopoints),
      E = pseudodata_exp,
      ...
    ),
    error = function(...){
      return(NULL)
    }
  )

  if(is.null(result)){
    return(tibble_row(
      Fit = list(NULL),
      IntMean = NA_real_,
      IntSD = NA_real_,
      Int025 = NA_real_,
      Int975 = NA_real_,
      RangeMean = NA_real_,
      RangeSD = NA_real_,
      Range025 = NA_real_,
      Range975 = NA_real_,
      SigMean = NA_real_,
      SigSD = NA_real_,
      Sig025 = NA_real_,
      Sig975 = NA_real_,
      Prediction = list(NULL),
      PredictionSD = list(NULL),
      MSPE = NA_real_,
      APV = NA_real_,
      MedPV = NA_real_,
      MaxPV = NA_real_,
      Area = NA_real_,
      SurveyProp = NA_real_
    ))
  }

  gpmap <- im(
    t(inla.mesh.project(rect_R_proj, result$summary.random[[1]]$mean)) +
      result$summary.fixed['intercept', 'mean'],
    xrange = rect_R$x,
    yrange = rect_R$y
  )

  return(tibble_row(
    Fit = list(if(save_fit) result),
    IntMean = result$summary.fixed['intercept', 'mean'],
    IntSD = result$summary.fixed['intercept', 'sd'],
    Int025 = result$summary.fixed['intercept', '0.025quant'],
    Int975 = result$summary.fixed['intercept', '0.975quant'],
    RangeMean = result$summary.hyperpar['Range for idx', 'mean'],
    RangeSD = result$summary.hyperpar['Range for idx', 'sd'],
    Range025 = result$summary.fixed['Range for idx', '0.025quant'],
    Range975 = result$summary.fixed['Range for idx', '0.975quant'],
    SigMean = result$summary.hyperpar['Stdev for idx', 'mean'],
    SigSD = result$summary.hyperpar['Stdev for idx', 'sd'],
    Sig025 = result$summary.fixed['Stdev for idx', '0.025quant'],
    Sig975 = result$summary.fixed['Stdev for idx', '0.975quant'],
    Prediction = list(if(save_pred) result$summary.random[[1]]$mean),
    PredictionSD = list(if(save_pred) result$summary.random[[1]]$sd),
    MSPE = mean((log(attr(obs_ppp, 'Lambda')) - gpmap)^2),
    APV = sum(mesh_weights * result$summary.random[[1]]$sd^2) / sum(mesh_weights),
    MedPV = if(any(is.na(result$summary.random[[1]]$sd))){
      NA_real_
    }else{
      weighted.median(result$summary.random[[1]]$sd^2, rect_R_nodes_area, na.rm = TRUE)
    },
    MaxPV = max(result$summary.random[[1]]$sd^2),
    Area = area(Window(obs_ppp)),
    SurveyProp = area(Window(obs_ppp)) / area(rect_R)
  ))
}


#}}}##############################

# vim: foldmethod=marker:
