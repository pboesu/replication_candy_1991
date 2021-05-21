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
candy_cm_pars = tibble(ref_estimate = c(c(5.49,5.49+3.90,5.49+6.74,5.49+10.21,5.49+15.77,5.49+21.76,-0.0457),
                                        c(3.32,3.32+2.53,3.32+4.40,3.32+6.71,3.32+10.14,3.32+14.20,-0.0307)),
                       parameter = rep(c(paste('a', 1:6, sep = ''), 'b'),2),
                       link = c(rep('logit', 7),rep('cloglog',7)))

# calculate CV and correlation between start and estimate for converged estimates
simulation_results %>% mutate(good = value < 95 & link == 'logit' | value < 100 & link == 'cloglog') %>% group_by(link) %>% summarize(finite_nll = sum(convergence == 0) / n(), good_nll_all = sum(good, na.rm = T)/n(), good_nll = sum(good, na.rm = T)/sum(convergence == 0))
simulation_results  %>%
  mutate(good = value < 95 & link == 'logit' | value < 100 & link == 'cloglog') %>%
  group_by(parameter, link) %>% 
  summarize(sim_range = paste(round(range(init),digits = 1), collapse = ' - '),
            conv_range = paste(round(range(init[convergence == 0], na.rm = TRUE),digits = 1), collapse = ' - ')) %>% arrange(desc(link), parameter) -> results_unfiltered

simulation_results %>% filter(value < 95 & link == 'logit' | value < 100 & link == 'cloglog') %>% group_by(parameter, link) %>% summarize(filt_range = paste(round(range(init),digits = 1), collapse = ' - '),cv = sd(estimate, na.rm = TRUE) / abs(mean(estimate, na.rm = TRUE))*100, corr_init = cor(init, estimate, use = 'complete'), corr_p = cor.test(init, estimate, use = 'complete')$p.value) %>% arrange(desc(link), parameter) -> results_filtered

# assemble table for manuscript
parameter_sens_table <- bind_cols(results_unfiltered, results_filtered[,3:6]) %>%
  knitr::kable(digits = 3, format = 'latex', align = 'c', row.names = FALSE,
               col.names = c('Par.', 'Link','Sim. range','Conv. range','Filtered range','CV (\\%)','$\\rho^{o}$','$P_\\rho$'),
               linesep = c('', '', '','', '', '', '\\addlinespace'), escape = FALSE, booktabs = TRUE)
cat(parameter_sens_table, file = 'outputs/parameter_sens_table.tex')


simulation_results %>%  tidyr::pivot_wider(names_from = parameter, values_from = c(estimate,init)) -> pivoted_sens 

# make a lookup table of parameter names for printing with plotmath and latex
par_formats = tibble(parameter_simple = c(paste('a', 1:6, sep = ''), 'b'),
                     parameter_plotmath = c(paste('alpha[',1:6,']', sep=''),'beta'),
                     parameter_latex = c(paste('$\\alpha_',1:6,'$', sep=''),'$\\beta$'),
                     init_latex = c(paste('$\\alpha^{o}_',1:6,'$', sep=''),'$\\beta^o$'),
                     est_latex = c(paste('$\\hat{\\alpha}_',1:6,'$', sep=''),'$\\hat{\\beta}$'))


simulation_results %>% left_join(par_formats, by = c('parameter' = 'parameter_simple')) %>% ggplot(aes(x = init, y = estimate, col = log(value))) + geom_point() + facet_wrap(link~parameter_plotmath, scales = 'free',nrow = 2, labeller = label_parsed) + geom_hline(aes(yintercept = ref_estimate), data = left_join(candy_cm_pars, par_formats, by = c('parameter' = 'parameter_simple')), col = "#E7B800", lty = 2, lwd = 0.7) + theme_classic() + theme(legend.position = 'bottom') + labs(color = 'log(negative log-likelihood)') + xlab('initial value')
ggsave('figures/figS1_initial_value_sensitivity_unfiltered.pdf', width = 10, height = 5.5)
# filter estimates for 'good fits' to better show variation within the bulk of the results
pivoted_sens %>% filter(value < 95 & link == 'logit') -> pivoted_sens_good_logit
pivoted_sens %>% filter(value < 100 & link == 'cloglog') -> pivoted_sens_good_cloglog
# pairwise plots of parameter estimates against initial values, and of initial values against each other
my_cols <- c("#E7B800","#00AFBB")  
pdf('figures/figS2_initial_value_sensitivity_logit.pdf', width = 8, height = 8)
par(mfrow = c(7,7),
    mar = c(2.3,2,0.5,1),
    mgp = c(1.2,0,0),
    new = F,
    tcl = 0.2,
    cex.lab = 1.1,
    cex.axis = 0.8)
for (i in 1:7){
  for (j in 1:7){
    if(i >= j){plot(pivoted_sens_good_logit[,c(11+j,4+i)], col = 'grey30',
                    pch = 19,  cex = 0.5,
                    xlim = range(pivoted_sens[,11+j]),
                    xlab = latex2exp::TeX(par_formats$init_latex[i]),
                    ylab = '')
      mtext(latex2exp::TeX(par_formats$est_latex[i]), side = 2, line = 1, las = 1, cex = 0.75)} else {
      plot(pivoted_sens[pivoted_sens$link=='logit',c(11+j,11+i)],
           col = my_cols[pivoted_sens[pivoted_sens$link=='logit',]$convergence+2],
           pch = 19,  cex = 0.5,
           xlab = latex2exp::TeX(par_formats$init_latex[i]),
           ylab = '')
        mtext(latex2exp::TeX(par_formats$init_latex[i]), side = 2, line = 1, las = 1, cex = 0.8)
    }
  }
}
dev.off()
pdf('figures/figS3_initial_value_sensitivity_cloglog.pdf', width = 8, height = 8)
par(mfrow = c(7,7),
    mar = c(2.3,2,0.5,1),
    mgp = c(1.2,0,0),
    new = F,
    tcl = 0.2,
    cex.lab = 1.1,
    cex.axis = 0.8)
for (i in 1:7){
  for (j in 1:7){
    if(i >= j){plot(pivoted_sens_good_cloglog[,c(11+j,4+i)], col = 'grey30',
                    pch = 19,  cex = 0.5,
                    xlim = range(pivoted_sens[,11+j]),
                    xlab = latex2exp::TeX(par_formats$init_latex[i]),
                    ylab = '')
      mtext(latex2exp::TeX(par_formats$est_latex[i]), side = 2, line = 1, las = 1, cex = 0.75)} else {
        plot(pivoted_sens[pivoted_sens$link=='cloglog',c(11+j,11+i)],
             col = my_cols[pivoted_sens[pivoted_sens$link=='cloglog',]$convergence+2],
             pch = 19,  cex = 0.5,
             xlab = latex2exp::TeX(par_formats$init_latex[i]),
             ylab = '')
        mtext(latex2exp::TeX(par_formats$init_latex[i]), side = 2, line = 1, las = 1, cex = 0.8)
      }
  }
}
dev.off()
