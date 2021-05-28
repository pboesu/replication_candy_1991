# fit sequential and cumulative models using VGAM
# Note that neither of the `vglm` models converged when using the cloglog link.
# This behaviour has been reported to the package maintainer as a potential bug.
# The model specifications are left in the code but commented out.

library(VGAM)
budworm_table <- readr::read_csv('data/budworm_table.csv', col_types = 'dddddddddd')

# # fitting the  cumulative model with cloglog link results in a false convergence - this has been reported to the package maintainer
# budworm_cumulative_cloglog_vglm <- vglm(cbind(stage1, stage2, stage3, stage4, stage5, stage6, stage7) ~ ddeg,
#                                         cumulative(reverse = FALSE, link = clogloglink(), parallel = TRUE),
#                                         coefstart = c(3.32,5.85,7.71,10.02,13.46,17,-0.03),
#                                         data = budworm_table,
#                                         control = vglm.control(checkwz = TRUE, trace = TRUE, stepsize = 1, maxit = 100, epsilon = 1e-10))# fails
# budworm_cumulative_cloglog_vglm_noinit <- vglm(cbind(stage1, stage2, stage3, stage4, stage5, stage6, stage7) ~ ddeg, cumulative(reverse = FALSE, link = clogloglink(), parallel = TRUE), data = budworm_table)

budworm_cumulative_logit_vglm <- vglm(cbind(stage1, stage2, stage3, stage4, stage5, stage6, stage7) ~ ddeg,
                                      cumulative(reverse = FALSE, link = "logitlink", parallel = TRUE), data = budworm_table) # reproduces parameter estimates in table 2 of candy

# save parameter estimates
candy_vglm_cm_estimates <- as.data.frame(rbind(unname(coef(budworm_cumulative_logit_vglm)), rep(NA,7)))
candy_vglm_cm_estimates$model <- "cumulative"
candy_vglm_cm_estimates$link <- c("logit","cloglog")
candy_vglm_cm_estimates$fit <- c("R \\verb+vglm+")
candy_vglm_cm_estimates$eqn <- c(NA,NA)
saveRDS(candy_vglm_cm_estimates, "outputs/candy_vglm_cm_estimates.RDS")

# fit sequential model
budworm_sratio_logit_vglm <- vglm(cbind(stage1, stage2, stage3, stage4, stage5, stage6, stage7) ~ ddeg,
                                  sratio(link = "logitlink", parallel = FALSE), data = budworm_table)# this function call will succeed but produce warnings about replacing working weights to prevent numerical underflow.

# # fitting the sratio model with clogloglink fails - this has been reported to the package maintainer
# budworm_sratio_cloglog_vglm <- try(vglm(cbind(stage1, stage2, stage3, stage4, stage5, stage6, stage7) ~ ddeg_cent, sratio(link = clogloglink(), parallel = FALSE), coefstart = coef(budworm_sratio_logit_vglm), data = budworm_table, control = vglm.control(checkwz = TRUE)))#

# save parameter estimates
candy_vglm_sm_estimates <- as.data.frame(rbind(unname(coef(budworm_sratio_logit_vglm)[1:6]),unname(coef(budworm_sratio_logit_vglm)[7:12]),rep(NA,6),rep(NA,6)))
candy_vglm_sm_estimates$model <- "sequential"
candy_vglm_sm_estimates$link <- c("logit","logit","cloglog","cloglog")
candy_vglm_sm_estimates$fit <- c("R \\verb+vglm+")
candy_vglm_sm_estimates$par <- c("$\\beta_{0j}$","$\\beta_{1j}$","$\\beta_{0j}$","$\\beta_{1j}$")
saveRDS(candy_vglm_sm_estimates, "outputs/candy_vglm_sm_estimates.RDS")
