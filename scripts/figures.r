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

rect_datasets <- readRDS('rect_data.rds')


#}}}##################################

allplans <- readRDS('rect_plans.rds')

plotplans <- c(
  'Hilbert000180',
  'Inhib000171',
  'LHS-TSP000259',
  'RPM001107',
  'Serp000277',
  'SRS000214',
  'Sys000013'
)

for(planid in plotplans){
  thisplan <- allplans %>% filter(PlanID == planid)
  pdf(paste0('../writeup/', planid, '.pdf'), width = 9, height = 4)
  par(mar = c(0, 0, 2, 0))
  plot(thisplan$Plan[[1]], main = paste(thisplan$Scheme, 'Example'))
  plot(rect_R, border = 'grey', add = TRUE)
  dev.off()
}

# Create combinations of plans and data.
#fit_design <- expand.grid(PlanID = allplans$PlanID, DataID = rect_datasets$DataID)

rect_results <- readRDS('rect_results.rds')

stopCluster(cl)

# vim: foldmethod=marker:
