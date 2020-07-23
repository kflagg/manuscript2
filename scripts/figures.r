# Load packages, create functions, create designs.
library(ggplot2)
theme_set(theme_classic())
theme_update(plot.title = element_text(hjust = 0.5))
source('functions.r')

stopCluster(cl)


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
allplans <- readRDS('../data/rect_plans.rds')


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


# Read the model fitting results and prepare for plotting.
rect_results <- readRDS('../data/rect_results.rds') %>%
  left_join(allplans %>% select(Scheme, Subscheme, PlanID, Distance, Segments)) %>%
  group_by(DataID, Subscheme) %>%
  mutate(AvgDistance = mean(Distance, na.rm = TRUE)) %>%
  ungroup %>%
  mutate(Variant = case_when(
    Subscheme %in% c('Inhib1', 'Inhib2', 'Inhib3', 'Inhib4') ~ '10% pairs',
    Subscheme %in% c('Inhib5', 'Inhib6', 'Inhib7', 'Inhib8') ~ '20% pairs',
    Subscheme %in% c('Serp1', 'Serp2', 'Serp3', 'Serp4') ~ '5 zigzags',
    Subscheme %in% c('Serp5', 'Serp6', 'Serp7', 'Serp8') ~ '8 zizags',
    TRUE ~ NA_character_
  ))

# Summarize the results.
rect_summary <- rect_results %>%
  left_join(allplans %>% select(Scheme, Subscheme, PlanID, Distance)) %>%
  group_by(DataID, Subscheme) %>%
  summarize(
    Scheme = unique(Scheme),
    AvgAPV = mean(APV, na.rm = TRUE),
    SDAPV = sd(APV, na.rm = TRUE),
    MinAPV = min(APV, na.rm = TRUE),
    Q1APV = quantile(APV, 0.25, na.rm = TRUE),
    MedAPV = median(APV, na.rm = TRUE),
    Q3APV = quantile(APV, 0.75, na.rm = TRUE),
    MaxAPV = max(APV, na.rm = TRUE),
    AvgMaxPV = mean(MaxPV, na.rm = TRUE),
    SDMaxPV = sd(MaxPV, na.rm = TRUE),
    MinMaxPV = min(MaxPV, na.rm = TRUE),
    Q1MaxPV = quantile(MaxPV, 0.25, na.rm = TRUE),
    MedMaxPV = median(MaxPV, na.rm = TRUE),
    Q3MaxPV = quantile(MaxPV, 0.75, na.rm = TRUE),
    MaxMaxPV = max(MaxPV, na.rm = TRUE),
    AvgMSPE = mean(MSPE, na.rm = TRUE),
    SDMSPE = sd(MSPE, na.rm = TRUE),
    MinMSPE = min(MSPE, na.rm = TRUE),
    Q1MSPE = quantile(MSPE, 0.25, na.rm = TRUE),
    MedMSPE = median(MSPE, na.rm = TRUE),
    Q3MSPE = quantile(MSPE, 0.75, na.rm = TRUE),
    MaxMSPE = max(MSPE, na.rm = TRUE),
    AvgDistance = mean(Distance),
    SDDistance = sd(Distance)
  ) %>%
  ungroup


# Examine which data/plan combinations could not be fit.
rect_results %>%
  filter(is.na(IntMean)) %>%
  print


# Plot APV and MSPE. Focus on Cluster000004 and LGCP000004 in the paper.
for(thisdataset in rect_datasets$DataID){
  png(paste0('../writeup/APV-', thisdataset, '.png'), width = 9, height = 4, units = 'in', res = 300)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = APV, x = Distance, col = Variant)) +
    geom_line(aes(x = AvgDistance), stat = 'summary', fun = median) +
    geom_point(alpha = 0.25) +
    facet_wrap(~Scheme) +
    scale_y_log10() +
    ggtitle('Average Prediction Variance vs Distance Surveyed')
  )
  dev.off()

  png(paste0('../writeup/MSPE-', thisdataset, '.png'), width = 9, height = 4, units = 'in', res = 300)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = MSPE, x = Distance, col = Variant)) +
    geom_line(aes(x = AvgDistance), stat = 'summary', fun = median) +
    geom_point(alpha = 0.25) +
    facet_wrap(~Scheme) +
    scale_y_log10() +
    ggtitle('MSPE vs Distance Surveyed')
  )
  dev.off()

  png(paste0('../writeup/APV-Inhib-', thisdataset, '.png'), width = 9, height = 4, units = 'in', res = 300)
  print(
    rect_results %>%
    filter(DataID == thisdataset, Scheme == 'Inhib') %>%
    left_join(inhib_design) %>%
    mutate(`Proportion Pairs` = paste0(100 * prop_pairs, '%')) %>%
    rename(`Total Transects` = num_xsects) %>%
    ggplot(aes(y = APV, x = `Total Transects`, col = `Proportion Pairs`)) +
    geom_line(stat = 'summary', fun = median) +
    geom_point(alpha = 0.5, position = position_jitterdodge(dodge.width = 5, jitter.width = 2, jitter.height = 0)) +
    scale_x_continuous(breaks = unique(inhib_design$num_xsects)) +
    scale_y_log10() +
    ggtitle('Average Prediction Variance for Inhibitory Plus Close Pairs Designs')
  )
  dev.off()

  png(paste0('../writeup/MSPE-Inhib-', thisdataset, '.png'), width = 9, height = 4, units = 'in', res = 300)
  print(
    rect_results %>%
    filter(DataID == thisdataset, Scheme == 'Inhib') %>%
    left_join(inhib_design) %>%
    mutate(`Proportion Pairs` = paste0(100 * prop_pairs, '%')) %>%
    rename(`Total Transects` = num_xsects) %>%
    ggplot(aes(y = MSPE, x = `Total Transects`, col = `Proportion Pairs`)) +
    geom_line(stat = 'summary', fun = median) +
    geom_point(alpha = 0.5, position = position_jitterdodge(dodge.width = 5, jitter.width = 2, jitter.height = 0)) +
    scale_x_continuous(breaks = unique(inhib_design$num_xsects)) +
    scale_y_log10() +
    ggtitle('MSPE for Inhibitory Plus Close Pairs Designs')
  )
  dev.off()

  png(paste0('../writeup/APV-Serp-', thisdataset, '.png'), width = 9, height = 4, units = 'in', res = 300)
  print(
    rect_results %>%
    filter(DataID == thisdataset, Scheme == 'Serp') %>%
    left_join(serp_design) %>%
    mutate(
      Zigzags = factor(serp_num),
    ) %>%
    rename(
      `Total Transects` = num_xsects
    ) %>%
    ggplot(aes(y = APV, x = `Total Transects`, col = Zigzags)) +
    geom_line(stat = 'summary', fun = median) +
    geom_point(alpha = 0.5, position = position_jitterdodge(dodge.width = 3, jitter.width = 2, jitter.height = 0)) +
    scale_x_continuous(breaks = unique(serp_design$num_xsects)) +
    scale_y_log10() +
    ggtitle('Average Prediction Variance for Serpentine Transect Designs')
  )
  dev.off()

  png(paste0('../writeup/MSPE-Serp-', thisdataset, '.png'), width = 9, height = 4, units = 'in', res = 300)
  print(
    rect_results %>%
    filter(DataID == thisdataset, Scheme == 'Serp') %>%
    left_join(serp_design) %>%
    mutate(
      Zigzags = factor(serp_num),
    ) %>%
    rename(
      `Total Transects` = num_xsects
    ) %>%
    ggplot(aes(y = MSPE, x = `Total Transects`, col = Zigzags)) +
    geom_line(stat = 'summary', fun = median) +
    geom_point(alpha = 0.5, position = position_jitterdodge(dodge.width = 3, jitter.width = 2, jitter.height = 0)) +
    scale_x_continuous(breaks = unique(serp_design$num_xsects)) +
    scale_y_log10() +
    ggtitle('MSPE for Serpentine Transect Designs')
  )
  dev.off()

#  rect_results %>%
#    filter(DataID == thisdataset, APV > 160) %>%
#    print
}

# Explore some of the high-MSPE vs low-MSPE fits.
mspe_focus <- c(
  'Hilbert000131',
  'Hilbert000130',
  'Inhib000184',
  'Inhib000138',
  'Inhib000534',
  'Inhib000514',
  'LHS-TSP000143',
  'LHS-TSP000131',
  'Serp000148',
  'Serp000150',
  'Serp000553',
  'SRS000179',
  'SRS000187',
  'Sys000174',
  'Sys000127'
)
mspe_results <- rect_results %>%
  filter(DataID == 'LGCP000004', PlanID %in% mspe_focus) %>%
  mutate(`MSPE Cluster` = ifelse(MSPE > 100, 'High', 'Low'))

rect_datasets %>%
  filter(DataID == 'LGCP000004') %>%
  `$`('Data') %>%
  `[[`(1) %>%
  attr('Lambda') %>%
  log %>%
  plot(main = 'Realized Log-Intensity of LGCP000004')

for(thisplan in mspe_focus){
  mspe_results %>%
  filter(DataID == 'LGCP000004', PlanID == thisplan) %>%
  `$`('Prediction') %>%
  `[[`(1) %>%
  inla.mesh.project.inla.mesh.projector(projector = rect_R_proj) %>% t %>%
  im(xrange = rect_R$x, yrange = rect_R$y) %>%
  plot(main = sprintf('Prediction Surface for LGCP000004, %s\n(MSPE = %.2f)',
                      thisplan,
                      mspe_results %>%
                      filter(DataID == 'LGCP000004', PlanID == thisplan) %>%
                      `[`(1, 'MSPE')
                      ), ribsep = 0.05)
  allplans %>%
  filter(PlanID == thisplan) %>%
  `$`('Plan') %>%
  `[[`(1) %>%
  plot(add = TRUE)
Sys.sleep(2)
}
