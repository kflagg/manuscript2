# Load packages, create functions, create designs.
library(ggplot2)
theme_set(theme_classic())
theme_update(plot.title = element_text(hjust = 0.5))

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
WP_MARGIN <- XSECT_WIDTH / 2 # Minimum distance between waypoints and site boundary.

# Mesh parameters to experiment with.
MAX_EDGE_LENGTH <- 50
MAX_EDGE_EXT <- 100
FINE_EDGE_LENGTH <- 20

# Graphics parameters.
NPIX_X <- 250 # Should be 500 but rLGCP to has trouble when NPIX_X != NPIX_Y.
NPIX_Y <- 250


# Define derived parameters
XSECT_RADIUS <- XSECT_WIDTH / 2 # Half the transect width
XSECT_NUM_INITIAL <- floor(DIST_INITIAL / XSECT_LENGTH_MAX)


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


# Function to fit model to observed data.
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


stopCluster(cl)

options(scipen = 5)


# Read the data and plans.
rect_datasets <- readRDS('../data/rect_data.rds')
allplans <- readRDS('../data/rect_plans.rds')


# Parellel transects at one end.
ex_plan <- allplans %>%
  filter(PlanID == 'Serp000148') %>%
  `$`('Plan') %>%
  `[[`(1)

finer_results <- bind_rows(lapply(seq_len(nrow(rect_datasets)), function(r){
  obs_ppp <- sample_ppp(rect_datasets$Data[[r]], ex_plan, XSECT_RADIUS)
  return(bind_cols(
    tibble_row(DataID = rect_datasets$DataID[r], PlanID = 'Finer'),
    Fit = model_fit(rect_R_formula, obs_ppp, rect_R_mesh, rect_dual_tess, rect_R_proj, rect_prior_fixed)
  ))}))


for(thisdataset in rect_datasets$DataID){
  # Plot the mesh.
  pdf(paste0('../graphics/finermesh-', thisdataset, '.pdf'), width = 9, height = 4)
  par(mar = c(0, 1, 2, 1))
  rect_datasets %>%
    filter(DataID == thisdataset) %>%
    `$`('Data') %>%
    `[[`(1) %>%
    attr('Lambda') %>%
    log %>%
    plot(main = 'Finer Mesh over Realized Log-Intensity', ribbon = FALSE)
  plot(rect_dual_tess, add = TRUE, border = '#808080', do.col = TRUE,
       values = rep(1, rect_dual_tess$n), col = '#ffffff40')
  plot(rect_R_mesh_net, add = TRUE, col = '#00000080')
  points(rect_R_mesh_loc[,], pch = 20)
  dev.off()

  pdf(paste0('../graphics/lambda-Finer-', thisdataset, '.pdf'), width = 9, height = 4)
  par(mar = c(0, 0, 2, 2))
  thisresult <- finer_results %>%
    filter(DataID == thisdataset)
  (thisresult$IntMean + inla.mesh.project(rect_R_proj, thisresult$Prediction[[1]])) %>%
    t %>%
    im(xrange = rect_R$x, yrange = rect_R$y) %>%
    plot(main = sprintf('Prediction Surface\n(MSPE = %.2f)',
                        finer_results %>%
                        filter(DataID == thisdataset) %>%
                        `[`(1, 'MSPE')
                        ), ribsep = 0.05, ribargs = list(las = 1))
  plot(rect_R_mesh_tess, border = '#00000010', add = TRUE)
  plot(ex_plan, col = '#ffffff40', add = TRUE)
  sample_ppp(
    rect_datasets %>% filter(DataID == thisdataset) %>% `$`('Data') %>% `[[`(1),
    ex_plan
  ) %>%
    points(col = '#ffffff80', bg = '#ffffff40', pch = 21, cex = 0.5)
  dev.off()

  pdf(paste0('../graphics/lambdaSD-Finer-', thisdataset, '.pdf'), width = 9, height = 4)
  par(mar = c(0, 0, 2, 2))
  plot(rect_dual_tess, border = '#80808020', do.col = TRUE,
       values = thisresult$PredictionSD[[1]], ribargs = list(las = 1),
       main = sprintf('Prediction SD of the GP\n(APV = %.2f)',
                      finer_results %>%
                      filter(DataID == thisdataset) %>%
                      `[`(1, 'APV')
                      ), ribsep = 0.05)
  plot(rect_R_mesh_tess, border = '#00000010', add = TRUE)
  plot(ex_plan, col = '#ffffff40', add = TRUE)
  sample_ppp(
    rect_datasets %>% filter(DataID == thisdataset) %>% `$`('Data') %>% `[[`(1),
    ex_plan
  ) %>%
    points(col = '#ffffff80', bg = '#ffffff40', pch = 21, cex = 0.5)
  dev.off()
}
