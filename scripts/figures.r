#########################
#{{{ Setup environment. #
#########################

# Load packages, create functions, create designs.
source('functions.r')


#}}}##################################
#{{{ Objects pertaining to the site. #
######################################

# Neat plot.
pdf('../writeup/mesh_full.pdf', width = 9, height = 4)
par(mar = c(0, 0, 2, 2))
plot(rect_dual_tess, border = '#80808020', do.col = TRUE,
     values = rect_R_nodes_area, ribargs = list(las = 1),
     main = 'Mesh with dual colored by node weight')
plot(rect_R_mesh_net, add = TRUE, col = '#00000080')
plot(rect_R, border = 'white', add = TRUE)
points(rect_R_mesh_loc[,], pch = 20)
dev.off()


#}}}##################################
#{{{ Objects pertaining to the data. #
######################################

rect_datasets <- readRDS('../data/rect_data.rds')


#}}}##################################

allplans <- readRDS('../data/rect_plans.rds')

plotplans <- c(
  'Hilbert000180',
  'Hilbert000237',
  'Inhib000171',
  'LHS-TSP000161',
  'Serp000124',
  'Serp000539',
  'SRS000176',
  'Sys000141'
)
plottitles <- c(
  'Hilbert design of order 4',
  'Hilbert design of order 5',
  'Inhibitory plus close pairs line transect design',
  'LHS-TSP design',
  'Serpentine transect design with 5 zigzags',
  'Serpentine transect design with 8 zigzags',
  'Simple random sample line transect design',
  'Systematic line transect design'
)

for(i in seq_along(plotplans)){
  thisplan <- allplans %>% filter(PlanID == plotplans[i])
  pdf(paste0('../writeup/', plotplans[i], '.pdf'), width = 6, height = 4)
  par(mar = c(3, 0, 2, 0))
  plot(thisplan$Plan[[1]], main = plottitles[i])
  plot(rect_R, border = 'grey', add = TRUE)
  dev.off()
}

# Create combinations of plans and data.
#fit_design <- expand.grid(PlanID = allplans$PlanID, DataID = rect_datasets$DataID)

rect_results <- readRDS('../data/rect_results.rds')

stopCluster(cl)

# vim: foldmethod=marker:
