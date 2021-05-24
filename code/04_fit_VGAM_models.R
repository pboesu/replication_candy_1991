# fit sequential and cumulative models using VGAM
# Note that neither of the `vglm` models converged when using the cloglog link.
# This behaviour has been reported to the package maintainer as a potential bug.
# The model specifications are left in the code but commented out.

library(VGAM)
budworm_table <- readr::read_csv('data/budworm_table.csv', col_types = 'dddddddddd')

# fitting the  cumulative model with cloglog link fails - this has been reported to the package maintainer
# budworm_cumulative_cloglog_vglm <- try(vglm(cbind(stage1, stage2, stage3, stage4, stage5, stage6, stage7) ~ ddeg,
#                                         cumulative(reverse = FALSE, link = clogloglink(bvalue = 1e-8), parallel = TRUE), coefstart = c(3,6,8,10,14,17,-0.04), data = budworm_table))# fails


budworm_cumulative_logit_vglm <- vglm(cbind(stage1, stage2, stage3, stage4, stage5, stage6, stage7) ~ ddeg,
                                      cumulative(reverse = FALSE, link = "logitlink", parallel = TRUE), data = budworm_table) # reproduces parameter estimates in table 2 of candy

# save parameter estimates
candy_vglm_cm_estimates <- as.data.frame(rbind(unname(coef(budworm_cumulative_logit_vglm)), rep(NA,7)))
candy_vglm_cm_estimates$model <- "cumulative"
candy_vglm_cm_estimates$link <- c("logit","cloglog")
candy_vglm_cm_estimates$fit <- c("R \\verb+vglm+")
candy_vglm_cm_estimates$eqn <- c(NA,NA)
saveRDS(candy_vglm_cm_estimates, "outputs/candy_vglm_cm_estimates.RDS")
# create predictions to compare with candy original fit
# pred_data <- expand.grid(ddeg = 0:800, stage = 1:7, total = 100)
# pred_data$alpha_index = pred_data$stage + 1

# fit sequential model
budworm_sratio_logit_vglm <- vglm(cbind(stage1, stage2, stage3, stage4, stage5, stage6, stage7) ~ ddeg,
                                  sratio(link = "logitlink", parallel = FALSE), data = budworm_table)# this function call will succeed but produce warnings about replacing working weights to prevent numerical underflow.

# fitting the sratio model with clogloglink fails - this has been reported to the package maintainer
# budworm_sratio_cloglog_vglm <- try(vglm(cbind(stage1, stage2, stage3, stage4, stage5, stage6, stage7) ~ ddeg_cent,
#                                        sratio(link = "clogloglink", parallel = FALSE), data = budworm_table, control = vglm.control(checkwz = TRUE)))#

# save parameter estimates
candy_vglm_sm_estimates <- as.data.frame(rbind(unname(coef(budworm_sratio_logit_vglm)[1:6]),unname(coef(budworm_sratio_logit_vglm)[7:12]),rep(NA,6),rep(NA,6)))
candy_vglm_sm_estimates$model <- "sequential"
candy_vglm_sm_estimates$link <- c("logit","logit","cloglog","cloglog")
candy_vglm_sm_estimates$fit <- c("R \\verb+vglm+")
candy_vglm_sm_estimates$par <- c("$\\beta_{0j}$","$\\beta_{1j}$","$\\beta_{0j}$","$\\beta_{1j}$")
saveRDS(candy_vglm_sm_estimates, "outputs/candy_vglm_sm_estimates.RDS")
