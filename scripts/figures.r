# Load packages, create functions, create designs.
library(ggplot2)
theme_set(theme_classic())
theme_update(plot.title = element_text(hjust = 0.5))
source('functions.r')

# Create subscheme identifiers based on the combination of scheme options.
srs_design$Subscheme <- paste0('SRS', as.numeric(factor(apply(
    srs_design %>% select(-PlanID), 1, paste, collapse = ','
  ))))
sys_design$Subscheme <- paste0('Sys', as.numeric(factor(apply(
    sys_design %>% select(-PlanID), 1, paste, collapse = ','
  ))))
serp_design$Subscheme <- paste0('Serp', as.numeric(factor(apply(
    serp_design %>% select(-PlanID), 1, paste, collapse = ','
  ))))
inhib_design$Subscheme <- paste0('Inhib', as.numeric(factor(apply(
    inhib_design %>% select(-PlanID), 1, paste, collapse = ','
  ))))
lhs_design$Subscheme <- paste0('LHS-TSP', as.numeric(factor(apply(
    lhs_design %>% select(-PlanID), 1, paste, collapse = ','
  ))))
hilb_design$Subscheme <- paste0('Hilbert', as.numeric(factor(apply(
    hilb_design %>% select(-PlanID), 1, paste, collapse = ','
  ))))
subschemes <- bind_rows(
  srs_design %>% select(PlanID, Subscheme),
  sys_design %>% select(PlanID, Subscheme),
  serp_design %>% select(PlanID, Subscheme),
  inhib_design %>% select(PlanID, Subscheme),
  lhs_design %>% select(PlanID, Subscheme),
  hilb_design %>% select(PlanID, Subscheme)
)


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


# Read the data and plans.
rect_datasets <- readRDS('../data/rect_data.rds')
allplans <- readRDS('../data/rect_plans.rds') %>%
  left_join(subschemes)


# Plot a selection of plans.
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


# Read the model fitting results.
rect_results <- readRDS('../data/rect_results.rds')


# Examine which data/plan combinations could not be fit.
rect_results %>%
  filter(sapply(Prediction, is.null)) %>%
  print


# Focus on LGCP000001 for now.
thisdataset <- 'LGCP000001'

  rect_results %>%
    filter(DataID == thisdataset) %>%
    left_join(allplans %>% select(Scheme, Subscheme, PlanID, Distance, Segments)) %>%
    ggplot(aes(y = APV, x = Distance)) +
    # TODO: find a better way to add the averages to the plot.
    geom_smooth(se = FALSE) +
    geom_point(aes(col = Segments), alpha = 0.25) +
    facet_wrap(~Scheme) +
    ylim(0, 160) +
    ggtitle('Average Prediction Variance vs Distance Surveyed')

  rect_results %>%
    filter(DataID == thisdataset, APV > 160) %>%
    print


stopCluster(cl)
