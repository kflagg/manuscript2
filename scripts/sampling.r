#########################
#{{{ Setup environment. #
#########################

# Load packages, create functions, create designs.
source('functions.r')


#}}}##################################
#{{{ Objects pertaining to the data. #
######################################

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

saveRDS(rect_datasets, 'rect_data.rds')


#}}}#######
#{{{ SRS. #
###########

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

# Combine all plans into one dataset.
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
saveRDS(allplans, 'rect_plans.rds')

# Create combinations of plans and data.
fit_design <- expand.grid(PlanID = allplans$PlanID, DataID = rect_datasets$DataID)

# Create a fresh cluster and prepare for fitting models.
stopCluster(cl)
cl <- makeCluster(ceiling(0.75 * detectCores()), outfile = '')
invisible(clusterEvalQ(cl, {
  library(spatstat)
  library(INLA)
  library(tibble)
  library(dplyr)
}))
clusterExport(cl, c('sample_ppp', 'model_fit', 'rect_R_spde', 'rect_R_formula', 'rect_dual_tess', 'rect_R_proj', 'rect_R_mesh', 'rect_R', 'allplans', 'fit_design', 'rect_datasets', 'rect_prior_fixed'))

# Fit models.
rect_results <- bind_rows(parLapply(cl, seq_len(nrow(fit_design)), function(r){
  thisplan <- allplans %>% filter(PlanID == fit_design$PlanID[r])
  thisdataset <- rect_datasets %>% filter(DataID == fit_design$DataID[r])
  obs_ppp <- sample_ppp(thisdataset$Data[[1]], thisplan$Plan[[1]], thisplan$xsect_radius)
  message(sprintf('Dataset %s, Plan %s, observed %d points',
                  thisdataset$DataID, thisplan$PlanID, obs_ppp$n))
  return(bind_cols(
    tibble_row(DataID = thisdataset$DataID, PlanID = thisplan$PlanID),
    Fit = model_fit(rect_R_formula, obs_ppp, rect_R_mesh, rect_dual_tess, rect_R_proj, rect_prior_fixed)
  ))}))

saveRDS(rect_results, 'rect_results.rds')

stopCluster(cl)

# vim: foldmethod=marker:
