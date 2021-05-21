#NB: package ordinal fits this without problems.... 
library(ordinal)
budworm_individuals <- readr::read_csv('data/budworm_individuals.csv', col_types = "dddd")

budworm_cm_cloglog_ordinal <- clm(factor(stage) ~ ddeg, data = budworm_individuals, link = 'cloglog')
summary(budworm_cm_cloglog_ordinal)
budworm_cm_logit_ordinal <- clm(factor(stage) ~ ddeg, data = budworm_individuals, link = 'logit')
summary(budworm_cm_logit_ordinal)
#save parameter estimates
candy_ordinal_cm_estimates <- as.data.frame(rbind(unname(coef(budworm_cm_logit_ordinal)), unname(coef(budworm_cm_cloglog_ordinal))))
candy_ordinal_cm_estimates$model <- "cumulative"
candy_ordinal_cm_estimates$link <- c("logit","cloglog")
candy_ordinal_cm_estimates$fit <- c("R \\verb+clm+")
candy_ordinal_cm_estimates$eqn <- c(NA,NA)
saveRDS(candy_ordinal_cm_estimates, "outputs/candy_ordinal_cm_estimates.RDS")
#try to calculate the poisson deviance for this model