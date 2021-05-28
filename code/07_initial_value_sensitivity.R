# script to assess sensitivity of results to starting values
# this script fits two models 300 times for random sets of initial values
# computations can take about 15 minutes on a desktop computer
library(ggplot2)
library(dplyr)
# load definitions of likelihood functions and related simulation functions
source('code/00_function_definitions.R')

# load data
budworm_counts <- readr::read_csv('data/budworm_counts.csv', col_types = "dddd")
budworm_counts$occ <- as.numeric(as.factor(budworm_counts$ddeg)) # encode sampling occasion
budworm_counts$alpha_index = budworm_counts$stage + 1 # encode an index variable to facilitate fitting the linear predictor
budworm_counts <- group_by(budworm_counts, ddeg) %>%
  arrange(ddeg,stage) %>% 
  mutate(cumulative = cumsum(count)) %>% # calculate cumulative counts across stages for cumulative likelihood
  mutate(n_star = ifelse(stage ==1, total, total - lag(cumulative)), # calculate N*_ij for GLM estimation of sequential model
         stage_factor = as.factor(stage)) 

set.seed(1027)
# simulate and estimate using 2x300 sets of initial values
# this can take about 15 minutes on a desktop computer
message('running 300 simulations')
sens_list <- lapply(1:300, cm_cand_sens)
message('running 300 simulations')
sens_list_cloglog <- lapply(1:300, cm_cand_sens, linkinv = function(x){VGAM::clogloglink(x, inverse = TRUE)})
# Some fits will produce warnings about producing NaNs in dpois(). This warning occurs, because the code explicitly checks for false model convergence under infinite likelihoods.

simulation_results <- bind_rows(bind_rows(sens_list) %>% mutate(link = 'logit'),
                                bind_rows(sens_list_cloglog) %>% mutate(link = 'cloglog'))

# write out optimal starting parameters
simulation_results %>% filter(link == 'cloglog') %>% filter(iter == iter[which.min(value)]) %>% pull(init) -> init_cloglog
simulation_results %>% filter(link == 'logit') %>% filter(iter == iter[which.min(value)]) %>% pull(init) -> init_logit
saveRDS(init_cloglog, 'outputs/init_cloglog.RDS')
saveRDS(init_logit, 'outputs/init_logit.RDS')

# reference parameter estimates from Candy 1991
candy_cm_pars = tibble(ref_estimate = c(c(5.49,5.49+3.90,5.49+6.74,5.49+10.21,5.49+15.77,5.49+21.76,-0.0457),
                                        c(3.32,3.32+2.53,3.32+4.40,3.32+6.71,3.32+10.14,3.32+14.20,-0.0307)),
                       parameter = rep(c(paste('a', 1:6, sep = ''), 'b'),2),
                       link = c(rep('logit', 7),rep('cloglog',7)))

# summarise convergence success
simulation_results %>%
  mutate(good = value < 93 & link == 'logit' | value < 100 & link == 'cloglog') %>%
  group_by(link) %>%
  summarize(finite_nll = sum(convergence == 0) / n(),
            good_nll_all = sum(good, na.rm = T)/n(),
            good_nll = sum(good, na.rm = T)/sum(convergence == 0),
            .groups = 'drop')

# calculate CV and correlation between start and estimate for converged estimates
simulation_results  %>%
  mutate(good = value < 95 & link == 'logit' | value < 100 & link == 'cloglog') %>%
  group_by(parameter, link) %>% 
  summarize(sim_range = paste(round(range(init),digits = 2), collapse = ' - '),
            conv_range = paste(round(range(init[convergence == 0], na.rm = TRUE),digits = 1), collapse = ' - '),
            .groups = 'drop') %>%
  arrange(desc(link), parameter) -> results_unfiltered

simulation_results %>%
  filter(value < 93 & link == 'logit' | value < 100 & link == 'cloglog') %>%
  group_by(parameter, link) %>%
  summarize(filt_range = paste(round(range(init),digits = 1), collapse = ' - '),
            cv = sd(estimate, na.rm = TRUE) / abs(mean(estimate, na.rm = TRUE))*100,
            corr_init = cor(init, estimate, use = 'complete'),
            corr_p = pmin(1,cor.test(init, estimate, use = 'complete')$p.value*7), #bonferroni correct p-values
            .groups = 'drop') %>%
  arrange(desc(link), parameter) -> results_filtered

# make a lookup table of parameter names for printing with plotmath and latex
par_formats = tibble(parameter_simple = c(paste('a', 1:6, sep = ''), 'b'),
                     parameter_plotmath = c(paste('alpha[',1:6,']', sep=''),'beta'),
                     parameter_latex = c(paste('$\\alpha_',1:6,'$', sep=''),'$\\beta$'),
                     init_latex = c(paste('$\\alpha^{o}_',1:6,'$', sep=''),'$\\beta^o$'),
                     est_latex = c(paste('$\\hat{\\alpha}_',1:6,'$', sep=''),'$\\hat{\\beta}$'))


# assemble table for manuscript
parameter_sens_table <- bind_cols(results_unfiltered, results_filtered[,3:6]) %>%
  left_join(select(par_formats, parameter_simple, parameter_latex), by = c('parameter' = 'parameter_simple')) %>%
  relocate(parameter_latex, before = parameter) %>%
  select(-parameter) %>%
  knitr::kable(digits = 3, format = 'latex', align = 'c', row.names = FALSE,
               col.names = c('Par.', 'Link','Sim. range','Conv. range','Filtered range','CV (\\%)','$\\rho^{o}$','$P_\\rho$'),
               linesep = c('', '', '','', '', '', '\\addlinespace'), escape = FALSE, booktabs = TRUE)
cat(parameter_sens_table, file = 'outputs/parameter_sens_table.tex')

# pivot results for easier plotting
simulation_results %>%
  tidyr::pivot_wider(names_from = parameter, values_from = c(estimate,init)) -> pivoted_sens 

# make a lookup table of parameter names for printing with plotmath and latex
par_formats = tibble(parameter_simple = c(paste('a', 1:6, sep = ''), 'b'),
                     parameter_plotmath = c(paste('alpha[',1:6,']', sep=''),'beta'),
                     parameter_latex = c(paste('$\\alpha_',1:6,'$', sep=''),'$\\beta$'),
                     init_latex = c(paste('$\\alpha^{o}_',1:6,'$', sep=''),'$\\beta^o$'),
                     est_latex = c(paste('$\\hat{\\alpha}_',1:6,'$', sep=''),'$\\hat{\\beta}$'))

# plot and save Figure S1
simulation_results %>%
  left_join(par_formats, by = c('parameter' = 'parameter_simple')) %>%
  ggplot(aes(x = init, y = estimate, col = log(value))) +
  geom_point(na.rm = TRUE) +
  facet_wrap(link~parameter_plotmath, scales = 'free',nrow = 2, labeller = label_parsed) +
  geom_hline(aes(yintercept = ref_estimate),
             data = left_join(candy_cm_pars, par_formats, by = c('parameter' = 'parameter_simple')),
             col = "#E7B800",
             lty = 2,
             lwd = 0.7) +
  theme_classic() +
  theme(legend.position = 'bottom') +
  labs(color = 'log(negative log-likelihood)') +
  xlab('initial value')
ggsave('figures/figS1_initial_value_sensitivity_unfiltered.pdf', width = 10, height = 5.5)

# filter estimates for 'good fits' to better show variation within the bulk of the results
pivoted_sens %>% filter(value < 93 & link == 'logit') -> pivoted_sens_good_logit
pivoted_sens %>% filter(value < 100 & link == 'cloglog') -> pivoted_sens_good_cloglog

# pairwise plots of parameter estimates against initial values, and of initial values against each other (Figures S2 and S3)
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


# # vgam start value search
# set.seed(4567)
# vgam_sens_list <- lapply(1:10000, cm_vgam_sens)
# vgam_results <- dplyr::bind_rows(vgam_sens_list)
# hist(vgam_results$value)
# 
# bind_rows(vgam_sens_list) %>% ggplot(aes(x = init, y = estimate, col = value)) + geom_point() + facet_wrap(~parameter) + ylim(-0.5,20) + geom_hline(aes(yintercept = ref_estimate, lty = link), data = candy_cm_pars) + abline()
# vgam_results[which.max(vgam_results$value),]
# 
# # cm model likelihood surface
# set.seed(4567)
# nll_surface_list <- lapply(1:10000, cm_nll_surface)
# nll_surface <- dplyr::bind_rows(nll_surface_list)
# nll_surface %>% ggplot(aes(x = init, y = log(value))) + geom_point() + facet_wrap(~parameter, scales = 'free')
# nll_wide <- tidyr::pivot_wider(nll_surface, values_from = init, names_from = parameter, names_prefix = 'init')
# ggplot(nll_wide, aes(x = initb, y = inita1, col = log(value))) + geom_point() + scale_color_viridis_c()
# ggplot(nll_wide, aes(x = inita1, y = log(value), col =initb)) + geom_point() + scale_color_viridis_c()
# 
# set.seed(4567)#b0- -0.1
# nll_surface_list <- lapply(1:10000, cm_nll_surface, linkinv = function(x){VGAM::clogloglink(x, inverse = TRUE)})
# nll_surface <- dplyr::bind_rows(nll_surface_list)
# nll_surface %>% ggplot(aes(x = init, y = log(value))) + geom_point() + facet_wrap(~parameter, scales = 'free')
# nll_wide <- tidyr::pivot_wider(nll_surface, values_from = init, names_from = parameter, names_prefix = 'init')
# ggplot(nll_wide, aes(x = initb, y = inita1, col = log(value))) + geom_point() + scale_color_viridis_c()
# ggplot(nll_wide, aes(x = inita1, y = log(value), col =initb)) + geom_point() + scale_color_viridis_c()