# Load packages, create functions, create designs.
source('functions.r')

# Read the datasets.
rect_datasets <- readRDS('../data/rect_data.rds')

# Mesh nodes will be sorted by y.
lambda_grid <- as.data.frame(attr(rect_datasets$Data[[1]], 'Lambda')) %>%
  arrange(y) %>%
  as.matrix
lambda_mesh <- inla.mesh.create(
  lattice = inla.mesh.lattice(
    unique(lambda_grid[,'x']), unique(lambda_grid[,'y'])
  ), refine = FALSE, extend = FALSE
)
proj_gridtomesh <- inla.mesh.projector(lambda_mesh, rect_R_mesh$loc[,1:2])

invisible(clusterEvalQ(cl, library(dplyr)))
clusterExport(cl, c('rect_datasets', 'proj_gridtomesh'))
lambda_at_nodes <- parSapply(cl, seq_len(nrow(rect_datasets)), function(r){return(
  inla.mesh.project(proj_gridtomesh,
    rect_datasets$Data[[r]] %>%
      attr('Lambda') %>%
      as.data.frame %>%
      arrange(y) %>%
      `$`('value')
  )
)})
colnames(lambda_at_nodes) <- rect_datasets$DataID
saveRDS(lambda_at_nodes, '../data/lambda_at_nodes.rds')

stopCluster(cl)
