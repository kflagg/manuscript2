# Load packages, create functions, create designs.
library(ggplot2)
theme_set(theme_classic())
theme_update(plot.title = element_text(hjust = 0.5))
source('functions.r')

stopCluster(cl)

options(scipen = 5)


# Read the data and plans.
rect_datasets <- readRDS('../data/rect_data.rds')
allplans <- readRDS('../data/rect_plans.rds')


# Parellel transects at one end.
original <- allplans %>%
  filter(PlanID == 'Serp000148') %>%
  `$`('Plan') %>%
  `[[`(1)
set.seed(64251)
aug <- serp(owin(c(1300, 1500), c(0, 700)), 1, -198, 11, XSECT_RADIUS)
aug_serp <- linnet(
  as.ppp(rbind(
    as.data.frame(vertices(original)),
    as.data.frame(vertices(aug))
  )[,c('x', 'y')], rect_R),
  as.matrix(bdiag(
    original$m,
    aug$m
  ))
)

aug_results <- bind_rows(lapply(seq_len(nrow(rect_datasets)), function(r){
  obs_ppp <- sample_ppp(rect_datasets$Data[[r]], aug_serp, XSECT_RADIUS)
  return(bind_cols(
    tibble_row(DataID = rect_datasets$DataID[r], PlanID = 'BadXsect'),
    Fit = model_fit(rect_R_formula, obs_ppp, rect_R_mesh, rect_dual_tess, rect_R_proj, rect_prior_fixed)
  ))}))


for(thisdataset in rect_datasets$DataID){
  pdf(paste0('../writeup/lambda-Aug-', thisdataset, '.pdf'), width = 9, height = 4)
  par(mar = c(0, 0, 2, 2))
  thisresult <- aug_results %>%
    filter(DataID == thisdataset)
  (thisresult$IntMean + inla.mesh.project(rect_R_proj, thisresult$Prediction[[1]])) %>%
    t %>%
    im(xrange = rect_R$x, yrange = rect_R$y) %>%
    plot(main = sprintf('Prediction Surface for %s, Augmented Design\n(MSPE = %.2f)',
                        thisdataset,
                        aug_results %>%
                        filter(DataID == thisdataset) %>%
                        `[`(1, 'MSPE')
                        ), ribsep = 0.05, ribargs = list(las = 1))
  plot(rect_R_mesh_tess, border = '#00000010', add = TRUE)
  plot(aug_serp, col = '#ffffff40', add = TRUE)
  sample_ppp(
    rect_datasets %>% filter(DataID == thisdataset) %>% `$`('Data') %>% `[[`(1),
    aug_serp
  ) %>%
    points(col = '#ffffff80', bg = '#ffffff40', pch = 21, cex = 0.5)
  dev.off()

  pdf(paste0('../writeup/lambdaSD-Aug-', thisdataset, '.pdf'), width = 9, height = 4)
  par(mar = c(0, 0, 2, 2))
  plot(rect_dual_tess, border = '#80808020', do.col = TRUE,
       values = thisresult$PredictionSD[[1]], ribargs = list(las = 1),
       main = sprintf('Prediction SD of the GP for %s, Augmented Design\n(APV = %.2f)',
                      thisdataset,
                      aug_results %>%
                      filter(DataID == thisdataset) %>%
                      `[`(1, 'APV')
                      ), ribsep = 0.05)
  plot(rect_R_mesh_tess, border = '#00000010', add = TRUE)
  plot(aug_serp, col = '#ffffff40', add = TRUE)
  sample_ppp(
    rect_datasets %>% filter(DataID == thisdataset) %>% `$`('Data') %>% `[[`(1),
    aug_serp
  ) %>%
    points(col = '#ffffff80', bg = '#ffffff40', pch = 21, cex = 0.5)
  dev.off()
}
