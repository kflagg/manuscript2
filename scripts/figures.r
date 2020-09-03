# Load packages, create functions, create designs.
library(ggplot2)
theme_set(theme_classic())
theme_update(plot.title = element_text(hjust = 0.5))
source('functions.r')

options(scipen = 5)


# Neat plot.
pdf('../graphics/mesh_full.pdf', width = 9, height = 4)
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
lambda_at_nodes <- readRDS('../data/lambda_at_nodes.rds')
log_lambda <- log(lambda_at_nodes)
intercept <- apply(log_lambda, 2, weighted.mean, rect_R_nodes_area, na.rm = TRUE)
allplans <- readRDS('../data/rect_plans.rds') %>%
  mutate(
    Segments = sapply(Lengths, length),
    Area = Distance * XSECT_WIDTH,
    SurveyProp = Area / area(rect_R),
    Effort = factor(case_when(
      Distance < 12000 ~ 'Low',
      Distance < 25000 ~ 'Medium',
      Distance < 45000 ~ 'High',
      TRUE ~ 'Very High'
    ), levels = c('Low', 'Medium', 'High',  'Very High'))
  )


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
  pdf(paste0('../graphics/', plotplans[i], '.pdf'), width = 6, height = 3)
  par(mar = c(0, 0, 1, 0))
  plot(thisplan$Plan[[1]], main = plottitles[i])
  plot(rect_R, border = 'grey', add = TRUE)
  dev.off()
}


# Read the model fitting results and prepare for plotting.
rect_results <- readRDS('../data/rect_results.rds') %>%
  left_join(allplans %>% select(Scheme, Subscheme, Effort, PlanID, Distance, Area, SurveyProp, Segments)) %>%
  group_by(DataID, Subscheme) %>%
  mutate(AvgDistance = mean(Distance, na.rm = TRUE)) %>%
  ungroup %>%
  mutate(
    IntMeanError = NA_real_,
    MedAPE = NA_real_,
    MedPV = NA_real_,
    Variant = case_when(
      Subscheme %in% c('Inhib1', 'Inhib2', 'Inhib3', 'Inhib4') ~ '10% pairs',
      Subscheme %in% c('Inhib5', 'Inhib6', 'Inhib7', 'Inhib8') ~ '20% pairs',
      Subscheme %in% c('Serp1', 'Serp2', 'Serp3', 'Serp4') ~ '5 zigzags',
      Subscheme %in% c('Serp5', 'Serp6', 'Serp7', 'Serp8') ~ '8 zizags',
      TRUE ~ NA_character_
    )
  )

clusterExport(cl, c('rect_results', 'log_lambda', 'intercept', 'rect_R_nodes_area'))
rect_results$IntMeanError <- parSapply(cl, seq_len(nrow(rect_results)), function(r){
  return(intercept[rect_results$DataID[r]] - rect_results$IntMean[r])
#  thesepreds <- rect_results$Prediction[[r]]
#  if(all(is.na(thesepreds))){return(NA_real_)}
#  return(weighted.median(abs(thesepreds - log_lambda[,rect_results$DataID[r]]),
#                         rect_R_nodes_area, na.rm = TRUE))
})

rect_results$MedAPE <- parSapply(cl, seq_len(nrow(rect_results)), function(r){
  thesepreds <- rect_results$Prediction[[r]]
  if(all(is.na(thesepreds))){return(NA_real_)}
  return(weighted.median(abs(thesepreds - log_lambda[,rect_results$DataID[r]]),
                         rect_R_nodes_area, na.rm = TRUE))
})

badpreds <- which(parSapply(cl, parLapply(cl, rect_results$Prediction, is.na), any))
rect_results$MedPV[-badpreds] <- parSapply(cl,
  rect_results$PredictionSD[-badpreds], weighted.median, rect_R_nodes_area
)^2

stopCluster(cl)


# Summarize the results.
rect_summary <- rect_results %>%
  left_join(allplans %>% select(Scheme, Subscheme, PlanID, Distance)) %>%
  group_by(DataID, Subscheme) %>%
  summarize(
    Scheme = unique(Scheme),
    Reps = sum(!is.na(IntMean)),
    AvgAPV = mean(APV, na.rm = TRUE),
    SDAPV = sd(APV, na.rm = TRUE),
    MinAPV = min(APV, na.rm = TRUE),
    Q1APV = quantile(APV, 0.25, na.rm = TRUE),
    MedAPV = median(APV, na.rm = TRUE),
    Q3APV = quantile(APV, 0.75, na.rm = TRUE),
    MaxAPV = max(APV, na.rm = TRUE),
    IQRAPV = IQR(APV, na.rm = TRUE),
    RangeAPV = max(APV, na.rm = TRUE) - min(APV, na.rm = TRUE),
    AvgMaxPV = mean(MaxPV, na.rm = TRUE),
    SDMaxPV = sd(MaxPV, na.rm = TRUE),
    MinMaxPV = min(MaxPV, na.rm = TRUE),
    Q1MaxPV = quantile(MaxPV, 0.25, na.rm = TRUE),
    MedMaxPV = median(MaxPV, na.rm = TRUE),
    Q3MaxPV = quantile(MaxPV, 0.75, na.rm = TRUE),
    MaxMaxPV = max(MaxPV, na.rm = TRUE),
    IQRMaxPV = IQR(MaxPV, na.rm = TRUE),
    RangeMaxPV = max(MaxPV, na.rm = TRUE) - min(MaxPV, na.rm = TRUE),
    AvgMSPE = mean(MSPE, na.rm = TRUE),
    SDMSPE = sd(MSPE, na.rm = TRUE),
    MinMSPE = min(MSPE, na.rm = TRUE),
    Q1MSPE = quantile(MSPE, 0.25, na.rm = TRUE),
    MedMSPE = median(MSPE, na.rm = TRUE),
    Q3MSPE = quantile(MSPE, 0.75, na.rm = TRUE),
    MaxMSPE = max(MSPE, na.rm = TRUE),
    IQRMSPE = IQR(MSPE, na.rm = TRUE),
    RangeMSPE = max(MSPE, na.rm = TRUE) - min(MSPE, na.rm = TRUE),
    AvgMedAPE = mean(MedAPE, na.rm = TRUE),
    SDMedAPE = sd(MedAPE, na.rm = TRUE),
    MinMedAPE = min(MedAPE, na.rm = TRUE),
    Q1MedAPE = quantile(MedAPE, 0.25, na.rm = TRUE),
    MedMedAPE = median(MedAPE, na.rm = TRUE),
    Q3MedAPE = quantile(MedAPE, 0.75, na.rm = TRUE),
    MaxMedAPE = max(MedAPE, na.rm = TRUE),
    IQRMedAPE = IQR(MedAPE, na.rm = TRUE),
    RangeMedAPE = max(MedAPE, na.rm = TRUE) - min(MedAPE, na.rm = TRUE),
    AvgMedPV = mean(MedPV, na.rm = TRUE),
    SDMedPV = sd(MedPV, na.rm = TRUE),
    MinMedPV = min(MedPV, na.rm = TRUE),
    Q1MedPV = quantile(MedPV, 0.25, na.rm = TRUE),
    MedMedPV = median(MedPV, na.rm = TRUE),
    Q3MedPV = quantile(MedPV, 0.75, na.rm = TRUE),
    MaxMedPV = max(MedPV, na.rm = TRUE),
    IQRMedPV = IQR(MedPV, na.rm = TRUE),
    RangeMedPV = max(MedPV, na.rm = TRUE) - min(MedPV, na.rm = TRUE),
    MinInt = min(IntMeanError, na.rm = TRUE),
    Q1Int = quantile(IntMeanError, 0.25, na.rm = TRUE),
    MedInt = median(IntMeanError, na.rm = TRUE),
    Q3Int = quantile(IntMeanError, 0.75, na.rm = TRUE),
    MaxInt = max(IntMeanError, na.rm = TRUE),
    IQRInt = IQR(IntMeanError, na.rm = TRUE),
    RangeInt = max(IntMeanError, na.rm = TRUE) - min(IntMeanError, na.rm = TRUE),
    AvgDistance = mean(Distance),
    SDDistance = sd(Distance),
    AvgArea = mean(Area),
    SDArea = sd(Area),
    AvgProp = mean(SurveyProp),
    SDProp = sd(SurveyProp),
    AvgSegments = mean(Segments),
    SDSegments = sd(Segments)
  ) %>%
  ungroup


# Examine which data/plan combinations could not be fit.
rect_results %>%
  filter(is.na(IntMean)) %>%
  print


# Plot APV and MSPE. Focus on Cluster000004 and LGCP000004 in the paper.
for(thisdataset in rect_datasets$DataID){
  png(paste0('../graphics/Int-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = IntMeanError, x = Distance, col = Variant)) +
    geom_line(aes(x = AvgDistance), stat = 'summary', fun = median) +
    geom_point(alpha = 0.25) +
    facet_wrap(~Scheme) +
    ggtitle('Error in Posterior Mean vs Distance Surveyed')
  )
  dev.off()

  png(paste0('../graphics/Int-notpaneled-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = IntMeanError, x = Distance, col = Variant)) +
    geom_line(aes(x = AvgDistance), stat = 'summary', fun = median) +
    geom_point(alpha = 0.25) +
    ggtitle('Error in Posterior Mean vs Distance Surveyed')
  )
  dev.off()

  png(paste0('../graphics/Int-effort-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = IntMeanError, x = Effort, col = Variant, group = Variant)) +
    geom_line(stat = 'summary', fun = median) +
    geom_point(alpha = 0.25) +
    facet_wrap(~Scheme) +
    ggtitle('Error in Posterior Mean vs Survey Effort')
  )
  dev.off()

  png(paste0('../graphics/Int-effort-notpaneled-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = IntMeanError, x = Effort, col = Variant, group = Variant)) +
    geom_line(stat = 'summary', fun = median) +
    geom_point(alpha = 0.25) +
    ggtitle('Error in Posterior Mean vs Survey Effort')
  )
  dev.off()

  png(paste0('../graphics/MaxPV-MSPE-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = MaxPV, x = MSPE, col = Scheme)) +
    geom_point(alpha = 0.25) +
    scale_x_log10() +
    scale_y_log10() +
    ggtitle('MaxPV vs MSPE')
  )
  dev.off()

  png(paste0('../graphics/APV-MSPE-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = APV, x = MSPE, col = Scheme)) +
    geom_point(alpha = 0.25) +
    scale_x_log10() +
    scale_y_log10() +
    ggtitle('APV vs MSPE')
  )
  dev.off()

  png(paste0('../graphics/MSPE-Segments-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    left_join(allplans %>% select(PlanID, Segments)) %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = MSPE, x = Segments, col = Scheme)) +
    geom_point(alpha = 0.25) +
    scale_x_log10() +
    scale_y_log10() +
    ggtitle('MSPE vs Number of Segments')
  )
  dev.off()

  png(paste0('../graphics/APV-Segments-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    left_join(allplans %>% select(PlanID, Segments)) %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = APV, x = Segments, col = Scheme)) +
    geom_point(alpha = 0.25) +
    scale_x_log10() +
    scale_y_log10() +
    ggtitle('APV vs Number of Segments')
  )
  dev.off()

  png(paste0('../graphics/MSPE-Coverage-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    left_join(allplans %>% select(PlanID, CoverageMaxDist)) %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = MSPE, x = CoverageMaxDist, col = Scheme)) +
    geom_point(alpha = 0.25) +
    scale_x_log10() +
    scale_y_log10() +
    ggtitle('MSPE vs Max Distance to Path')
  )
  dev.off()

  png(paste0('../graphics/APV-Coverage-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    left_join(allplans %>% select(PlanID, CoverageMaxDist)) %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = APV, x = CoverageMaxDist, col = Scheme)) +
    geom_point(alpha = 0.25) +
    scale_x_log10() +
    scale_y_log10() +
    ggtitle('APV vs Max Distance to Path')
  )
  dev.off()

  png(paste0('../graphics/APV-notpaneled-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = APV, x = Distance, col = Scheme)) +
    geom_line(aes(x = AvgDistance), stat = 'summary', fun = median) +
    geom_point(alpha = 0.25) +
    scale_y_log10() +
    ggtitle('APV vs Distance Surveyed')
  )
  dev.off()

  png(paste0('../graphics/APV-effort-notpaneled-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = APV, x = Effort, col = Scheme, group = Scheme)) +
    geom_line(stat = 'summary', fun = median) +
    geom_point(alpha = 0.25) +
    scale_y_log10() +
    ggtitle('APV vs Survey Effort')
  )
  dev.off()

  png(paste0('../graphics/MSPE-notpaneled-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = MSPE, x = Distance, col = Scheme)) +
    geom_line(aes(x = AvgDistance), stat = 'summary', fun = median) +
    geom_point(alpha = 0.25) +
    scale_y_log10() +
    ggtitle('MSPE vs Distance Surveyed')
  )
  dev.off()

  png(paste0('../graphics/MSPE-effort-notpaneled-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = MSPE, x = Effort, col = Scheme, group = Scheme)) +
    geom_line(stat = 'summary', fun = median) +
    geom_point(alpha = 0.25) +
    scale_y_log10() +
    ggtitle('MSPE vs Survey Effort')
  )
  dev.off()

  png(paste0('../graphics/APV-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = APV, x = Distance, col = Variant)) +
    geom_line(aes(x = AvgDistance), stat = 'summary', fun = median) +
    geom_point(alpha = 0.25) +
    facet_wrap(~Scheme) +
    scale_y_log10() +
    ggtitle('APV vs Distance Surveyed')
  )
  dev.off()

  png(paste0('../graphics/APV-effort-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = APV, x = Effort, col = Variant, group = Variant)) +
    geom_line(stat = 'summary', fun = median) +
    geom_point(alpha = 0.25) +
    facet_wrap(~Scheme) +
    scale_y_log10() +
    ggtitle('APV vs Survey Effort')
  )
  dev.off()

  png(paste0('../graphics/MedPV-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = MedPV, x = Distance, col = Variant)) +
    geom_line(aes(x = AvgDistance), stat = 'summary', fun = median) +
    geom_point(alpha = 0.25) +
    facet_wrap(~Scheme) +
    scale_y_log10() +
    ggtitle('Median Prediction Variance vs Distance Surveyed')
  )
  dev.off()

  png(paste0('../graphics/MedPV-effort-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = MedPV, x = Effort, col = Variant, group = Variant)) +
    geom_line(stat = 'summary', fun = median) +
    geom_point(alpha = 0.25) +
    facet_wrap(~Scheme) +
    scale_y_log10() +
    ggtitle('Median Prediction Variance vs Survey Effort')
  )
  dev.off()

  png(paste0('../graphics/MaxPV-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = MaxPV, x = Distance, col = Variant)) +
    geom_line(aes(x = AvgDistance), stat = 'summary', fun = median) +
    geom_point(alpha = 0.25) +
    facet_wrap(~Scheme) +
    scale_y_log10() +
    ggtitle('Maximum Prediction Variance vs Distance Surveyed')
  )
  dev.off()

  png(paste0('../graphics/MaxPV-effort-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = MaxPV, x = Effort, col = Variant, group = Variant)) +
    geom_line(stat = 'summary', fun = median) +
    geom_point(alpha = 0.25) +
    facet_wrap(~Scheme) +
    scale_y_log10() +
    ggtitle('Maximum Prediction Variance vs Survey Effort')
  )
  dev.off()

  png(paste0('../graphics/MSPE-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
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

  png(paste0('../graphics/MSPE-effort-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = MSPE, x = Effort, col = Variant, group = Variant)) +
    geom_line(stat = 'summary', fun = median) +
    geom_point(alpha = 0.25) +
    facet_wrap(~Scheme) +
    scale_y_log10() +
    ggtitle('MSPE vs Survey Effort')
  )
  dev.off()

  png(paste0('../graphics/MedAPE-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = MedAPE, x = Distance, col = Variant)) +
    geom_line(aes(x = AvgDistance), stat = 'summary', fun = median) +
    geom_point(alpha = 0.25) +
    facet_wrap(~Scheme) +
    scale_y_log10() +
    ggtitle('Median Absolute Prediction Error vs Distance Surveyed')
  )
  dev.off()

  png(paste0('../graphics/MedAPE-effort-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
  print(
    rect_results %>%
    filter(DataID == thisdataset) %>%
    ggplot(aes(y = MedAPE, x = Effort, col = Variant, group = Variant)) +
    geom_line(stat = 'summary', fun = median) +
    geom_point(alpha = 0.25) +
    facet_wrap(~Scheme) +
    scale_y_log10() +
    ggtitle('Median Absolute Prediction Error vs Survey Effort')
  )
  dev.off()

  png(paste0('../graphics/APV-Inhib-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
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

  png(paste0('../graphics/MSPE-Inhib-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
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

  png(paste0('../graphics/APV-Serp-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
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

  png(paste0('../graphics/MSPE-Serp-', thisdataset, '.png'), width = 9, height = 6, units = 'in', res = 600)
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
for(thisdataset in c('LGCP000004', 'Cluster000004')){
mspe_results <- rect_results %>%
  filter(DataID == thisdataset, PlanID %in% mspe_focus) %>%
  mutate(`MSPE Cluster` = ifelse(MSPE > 100, 'High', 'Low'))

pdf(paste0('../graphics/lambda-', thisdataset, '.pdf'), width = 9, height = 4)
par(mar = c(0, 0, 2, 2))
rect_datasets %>%
  filter(DataID == thisdataset) %>%
  `$`('Data') %>%
  `[[`(1) %>%
  attr('Lambda') %>%
  log %>%
  plot(main = 'Realized Log-Intensity',
       ribsep = 0.05, ribargs = list(las = 1))
rect_datasets %>%
  filter(DataID == thisdataset) %>%
  `$`('Data') %>%
  `[[`(1) %>%
  points(col = '#ffffff80', bg = '#ffffff40', pch = 21, cex = 0.5)
dev.off()

for(thisplan in mspe_focus){
  pdf(paste0('../graphics/lambda-', thisplan, '-', thisdataset, '.pdf'), width = 9, height = 4)
  par(mar = c(0, 0, 2, 2))
  thisresult <- mspe_results %>%
    filter(DataID == thisdataset, PlanID == thisplan)
  (thisresult$IntMean + inla.mesh.project(rect_R_proj, thisresult$Prediction[[1]])) %>%
    t %>%
    im(xrange = rect_R$x, yrange = rect_R$y) %>%
    plot(main = sprintf('Prediction Surface\n(MSPE = %.2f)',
                        mspe_results %>%
                        filter(DataID == thisdataset, PlanID == thisplan) %>%
                        `[`(1, 'MSPE')
                        ), ribsep = 0.05, ribargs = list(las = 1))
  plot(rect_R_mesh_tess, border = '#00000010', add = TRUE)
  allplans %>%
    filter(PlanID == thisplan) %>%
    `$`('Plan') %>%
    `[[`(1) %>%
    plot(col = 'white', add = TRUE)
  sample_ppp(
    rect_datasets %>% filter(DataID == thisdataset) %>% `$`('Data') %>% `[[`(1),
    allplans %>% filter(PlanID == thisplan) %>% `$`('Plan') %>% `[[`(1)
  ) %>%
    points(col = '#ffffff80', bg = '#ffffff40', pch = 21, cex = 0.5)
  dev.off()

  pdf(paste0('../graphics/lambdaSD-', thisplan, '-', thisdataset, '.pdf'), width = 9, height = 4)
  par(mar = c(0, 0, 2, 2))
  thisresult <- mspe_results %>%
    filter(DataID == thisdataset, PlanID == thisplan)
  plot(rect_dual_tess, border = '#80808020', do.col = TRUE,
       values = thisresult$PredictionSD[[1]], ribargs = list(las = 1),
       main = sprintf('Prediction SD of the GP\n(APV = %.2f)',
                      mspe_results %>%
                      filter(DataID == thisdataset, PlanID == thisplan) %>%
                      `[`(1, 'APV')
                      ), ribsep = 0.05)
  plot(rect_R_mesh_tess, border = '#00000010', add = TRUE)
  allplans %>%
    filter(PlanID == thisplan) %>%
    `$`('Plan') %>%
    `[[`(1) %>%
    plot(col = 'white', add = TRUE)
  sample_ppp(
    rect_datasets %>% filter(DataID == thisdataset) %>% `$`('Data') %>% `[[`(1),
    allplans %>% filter(PlanID == thisplan) %>% `$`('Plan') %>% `[[`(1)
  ) %>%
    points(col = '#ffffff80', bg = '#ffffff40', pch = 21, cex = 0.5)
  dev.off()
}
}

# Adjust the axes for the plot used in the manuscript.
thisdataset <- 'LGCP000004'
thisplan <- 'Serp000148'
mspe_results <- rect_results %>%
  filter(DataID == thisdataset, PlanID %in% mspe_focus) %>%
  mutate(`MSPE Cluster` = ifelse(MSPE > 100, 'High', 'Low'))
pdf(paste0('../graphics/lambda-', thisplan, '-', thisdataset, '.pdf'), width = 9, height = 4)
par(mar = c(0, 0, 2, 2))
thisresult <- mspe_results %>%
  filter(DataID == thisdataset, PlanID == thisplan)
(thisresult$IntMean + inla.mesh.project(rect_R_proj, thisresult$Prediction[[1]])) %>%
  t %>%
  im(xrange = rect_R$x, yrange = rect_R$y) %>%
  plot(main = sprintf('Prediction Surface\n(MSPE = %.2f)',
                      mspe_results %>%
                      filter(DataID == thisdataset, PlanID == thisplan) %>%
                      `[`(1, 'MSPE')
                      ), zlim = c(-100, 100), ribsep = 0.05, ribargs = list(las = 1))
plot(rect_R_mesh_tess, border = '#00000010', add = TRUE)
allplans %>%
  filter(PlanID == thisplan) %>%
  `$`('Plan') %>%
  `[[`(1) %>%
  plot(col = 'white', add = TRUE)
sample_ppp(
  rect_datasets %>% filter(DataID == thisdataset) %>% `$`('Data') %>% `[[`(1),
  allplans %>% filter(PlanID == thisplan) %>% `$`('Plan') %>% `[[`(1)
) %>%
  points(col = '#ffffff80', bg = '#ffffff40', pch = 21, cex = 0.5)
dev.off()


# Plot the lowest-MSPE and lowest-APV surface from each scheme.
thisdataset <- 'LGCP000004'
low_mspe <- c(
  rect_results %>%
    filter(DataID == thisdataset) %>%
    group_by(Scheme) %>%
    top_n(1, desc(MSPE)) %>%
    `$`('PlanID'),
  rect_results %>%
    filter(DataID == thisdataset) %>%
    group_by(Scheme) %>%
    top_n(1, desc(APV)) %>%
    `$`('PlanID')
) %>% unique
for(thisplan in low_mspe){
  pdf(paste0('../graphics/lambda-', thisplan, '-', thisdataset, '.pdf'), width = 9, height = 4)
  par(mar = c(0, 0, 2, 2))
  thisresult <- rect_results %>%
    filter(DataID == thisdataset, PlanID == thisplan)
  (thisresult$IntMean + inla.mesh.project(rect_R_proj, thisresult$Prediction[[1]])) %>%
    t %>%
    im(xrange = rect_R$x, yrange = rect_R$y) %>%
    plot(main = sprintf('Prediction Surface\n(MSPE = %.2f)',
                        rect_results %>%
                        filter(DataID == thisdataset, PlanID == thisplan) %>%
                        `[`(1, 'MSPE')
                        ), ribsep = 0.05, ribargs = list(las = 1))
  plot(rect_R_mesh_tess, border = '#00000010', add = TRUE)
  allplans %>%
    filter(PlanID == thisplan) %>%
    `$`('Plan') %>%
    `[[`(1) %>%
    plot(col = 'white', add = TRUE)
  sample_ppp(
    rect_datasets %>% filter(DataID == thisdataset) %>% `$`('Data') %>% `[[`(1),
    allplans %>% filter(PlanID == thisplan) %>% `$`('Plan') %>% `[[`(1)
  ) %>%
    points(col = '#ffffff80', bg = '#ffffff40', pch = 21, cex = 0.5)
  dev.off()

  pdf(paste0('../graphics/lambdaSD-', thisplan, '-', thisdataset, '.pdf'), width = 9, height = 4)
  par(mar = c(0, 0, 2, 2))
  thisresult <- rect_results %>%
    filter(DataID == thisdataset, PlanID == thisplan)
  plot(rect_dual_tess, border = '#80808020', do.col = TRUE,
       values = thisresult$PredictionSD[[1]], ribargs = list(las = 1),
       main = sprintf('Prediction SD of the GP\n(APV = %.2f)',
                      rect_results %>%
                      filter(DataID == thisdataset, PlanID == thisplan) %>%
                      `[`(1, 'APV')
                      ), ribsep = 0.05)
  plot(rect_R_mesh_tess, border = '#00000010', add = TRUE)
  allplans %>%
    filter(PlanID == thisplan) %>%
    `$`('Plan') %>%
    `[[`(1) %>%
    plot(col = 'white', add = TRUE)
  sample_ppp(
    rect_datasets %>% filter(DataID == thisdataset) %>% `$`('Data') %>% `[[`(1),
    allplans %>% filter(PlanID == thisplan) %>% `$`('Plan') %>% `[[`(1)
  ) %>%
    points(col = '#ffffff80', bg = '#ffffff40', pch = 21, cex = 0.5)
  dev.off()
}
