# supplementary script to assess sensitivity of results to starting values
library(ggplot2)
library(dplyr)

# load data
budworm_counts <- readr::read_csv('data/budworm_counts.csv', col_types = "dddd")
budworm_counts$occ <- as.numeric(as.factor(budworm_counts$ddeg)) # encode sampling occasion
budworm_counts$alpha_index = budworm_counts$stage + 1 # encode an index variable to facilitate fitting the linear predictor
budworm_counts <- group_by(budworm_counts, ddeg) %>%
  arrange(ddeg,stage) %>% 
  mutate(cumulative = cumsum(count)) %>% # calculate cumulative counts across stages for cumulative likelihood
  mutate(n_star = ifelse(stage ==1, total, total - lag(cumulative)), # calculate N*_ij for GLM estimation of sequential model
         stage_factor = as.factor(stage)) 

# likelihood function without numerical thresholds against underflow 
poisson_nll_cm_candy_inf <- function(par, data, linkinv = gtools::inv.logit){
  alpha <- c(-10, par[1:6], 5000)
  beta = par[7]
  pred_n <- rep(NA, nrow(data))
  for (i in 1:nrow(data)){
    pred_n[i] <- with(data,
                      total[i] * (linkinv(alpha[alpha_index[i]] + beta*(ddeg[i])) - linkinv(alpha[alpha_index[i]-1] + beta*(ddeg[i])))
    )
  }
  nll = sum(-1*dpois(data$count, pred_n, log=T))
  return(nll)
}

# function to repeatedly fit the cumlative logit model with randomised starting values
cm_cand_sens <- function(iter, progress = TRUE, linkinv = gtools::inv.logit){
  if (progress == TRUE & (iter %% 10 == 0)) print(paste(Sys.time(), iter))
  
  # draw random initial values
  alpha_start <- cumsum(runif(6,1,20))
  beta_start <- runif(1, -1,1)
  # estimate model parameters
  cm_candy_nll <- optim(par = c(alpha_start,beta_start), poisson_nll_cm_candy, data = budworm_counts, hessian = TRUE, control = list(trace=0), method = 'BFGS', linkinv = linkinv)
  # check that converged parameters actually have a finite likelihood
  if (is.finite(poisson_nll_cm_candy_inf(par = cm_candy_nll$par, data = budworm_counts, linkinv = linkinv))){
    sens_out <- tibble(init = c(alpha_start,beta_start),
                       estimate = cm_candy_nll$par,
                       parameter = c(paste('a', 1:6, sep = ''), 'b'),
                       convergence = cm_candy_nll$convergence,
                       iter = iter,
                       value = cm_candy_nll$value) 
    } else { # else return an indicator of false convergence
    sens_out <- tibble(init = c(alpha_start,beta_start),
    convergence = -1,
    parameter = c(paste('a', 1:6, sep = ''), 'b'),
    iter = iter)
    }
 return(sens_out)
  }

set.seed(1027)
# simulate and estimate using 2x300 sets of initial values
# this can take about 15 minutes on a desktop computer
sens_list <- lapply(1:300, cm_cand_sens)
sens_list_cloglog <- lapply(1:300, cm_cand_sens, linkinv = function(x){VGAM::clogloglink(x, inverse = TRUE)})
simulation_results <- bind_rows(bind_rows(sens_list) %>% mutate(link = 'logit'),
                                bind_rows(sens_list_cloglog) %>% mutate(link = 'cloglog'))

# reference parameter estimates from Candy 1991
candy_cm_logit_pars = tibble(ref_estimate = c(5.49,5.49+3.90,5.49+6.74,5.49+10.21,5.49+15.77,5.49+21.76,-0.0457), parameter = c(paste('a', 1:6, sep = ''), 'b'))

# calculate CV and correlation between start and estimate for converged estimates
simulation_results %>% group_by(parameter, link) %>% summarize(range = paste(round(range(init)), collapse = ' - '), cv = sd(estimate, na.rm = TRUE) / mean(estimate, na.rm = TRUE), corr_init = cor(init, estimate, use = 'complete'), corr_p = cor.test(init, estimate, use = 'complete')$p.value)
simulation_results %>% filter(value < 95 & link == 'logit' | value < 100 & link == 'cloglog') %>% group_by(parameter, link) %>% summarize(range = paste(round(range(init)), collapse = ' - '), cv = sd(estimate, na.rm = TRUE) / mean(estimate, na.rm = TRUE)*100, corr_init = cor(init, estimate, use = 'complete'), corr_p = cor.test(init, estimate, use = 'complete')$p.value)


simulation_results %>%  tidyr::pivot_wider(names_from = parameter, values_from = c(estimate,init)) -> pivoted_sens 

# 
simulation_results %>% ggplot(aes(x = init, y = estimate, col = log(value))) + geom_point() + facet_wrap(link~parameter, scales = 'free',nrow = 2)
# filter estimates for 'good fits' to better show variation within the bulk of the results
pivoted_sens %>% filter(value < 95 & link == 'logit') -> pivoted_sens_good_logit
pivoted_sens %>% filter(value < 100 & link == 'cloglog') -> pivoted_sens_good_cloglog
# pairwise plots of parameter estimates against initial values, and of initial values against each other
my_cols <- c("#E7B800","#00AFBB")  
pdf('figures/fig3_initial_value_sensitivity_logit.pdf', width = 8, height = 8)
par(mfrow = c(7,7),
    mar = c(2,2,1,1),
    mgp = c(1.1,0,0),
    new = F,
    tcl = 0.2,
    cex.lab = 1.1)
for (i in 1:7){
  for (j in 1:7){
    if(i >= j){plot(pivoted_sens_good_logit[,c(11+j,4+i)], col = 'grey30', pch = 19,  cex = 0.5, xlim = range(pivoted_sens[,11+j]))} else {
      plot(pivoted_sens[pivoted_sens$link=='logit',c(11+j,11+i)], col = my_cols[pivoted_sens[pivoted_sens$link=='logit',]$convergence+2], pch = 19,  cex = 0.5 )
    }
  }
}
dev.off()
pdf('figures/fig3_initial_value_sensitivity_cloglog.pdf', width = 8, height = 8)
par(mfrow = c(7,7),
    mar = c(2,2,1,1),
    mgp = c(1.1,0,0),
    new = F,
    tcl = 0.2,
    cex.lab = 1.1)
for (i in 1:7){
  for (j in 1:7){
    if(i >= j){plot(pivoted_sens_good_cloglog[,c(11+j,4+i)], col = 'grey30', pch = 19,  cex = 0.5, xlim = range(pivoted_sens[,11+j]))} else {
      plot(pivoted_sens[pivoted_sens$link=='cloglog',c(11+j,11+i)], col = my_cols[pivoted_sens[pivoted_sens$link=='cloglog',]$convergence+2], pch = 19,  cex = 0.5 )
    }
  }
}
dev.off()
