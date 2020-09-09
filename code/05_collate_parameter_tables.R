#create table of parameter estimates for cumulative model fits
library(dplyr)

#predictions and comparisons
candy_cm_logit_par = c(5.49,5.49+3.90,5.49+6.74,5.49+10.21,5.49+15.77,5.49+21.76,-0.0457)#transcribed from Candy 1991 Table 2 and alphas added up to match original parameterisation
candy_cm_cloglog_par = c(3.32,3.32+2.53,3.32+4.40,3.32+6.71,3.32+10.14,3.32+14.20,-0.0307)#transcribed from Candy 1991 Table 2
#candy_cm_logit_dennis_par = c(101.0,101.0+71.2,101.0+121.7,101.0+186.2,101.0+289.9,101.0+400.3,-0.8416)#transcribed from Candy 1991 Table 2

candy_estimates <- as.data.frame(rbind(candy_cm_logit_par,candy_cm_cloglog_par))
candy_estimates$model <- "cumulative"
candy_estimates$link <- c("logit","cloglog")
candy_estimates$fit <- c("Original \\citep{candy1991modeling}")
candy_estimates$eqn <- c("\\ref{eq:candy_cm_count_form}")

#load other estimates
candy_nll_estimates <- readRDS("outputs/candy_nll_estimates.RDS")
candy_ordinal_cm_estimates <- readRDS("outputs/candy_ordinal_cm_estimates.RDS")
candy_vglm_cm_estimates <- readRDS("outputs/candy_vglm_cm_estimates.RDS")

#rearrange and filter to create table for cumulative model with constant variance
cm_constant_var_table <- rbind(candy_estimates,candy_nll_estimates,candy_ordinal_cm_estimates,candy_vglm_cm_estimates) %>%
  filter(V1 < 50 | is.na(V1)) %>%
  arrange(desc(link),eqn) %>% 
  select(-eqn, -model) %>%
  knitr::kable(digits = 2, format = 'latex', align = 'c', row.names = FALSE,
               col.names = c(paste('$\\alpha_',1:6,'$', sep=''),'$\\beta$','Link','Method'),
               linesep = c('', '', '', '\\addlinespace'), escape = FALSE, booktabs = TRUE)


cat(cm_constant_var_table, file = 'outputs/cm_const_var_table.tex')
