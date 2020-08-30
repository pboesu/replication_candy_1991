#budworm example from Dennis et al. 1986 / Candy 1991
library(dplyr)
library(ggplot2)
library(brms)
#read raw data which has multiple columns side by side. columns are DDEG (degree days) TOT (total) LSF (life stage) NUM (count)
budworm_raw <- as.matrix(read.table("data/budworm_candy_1991_raw.txt"))
#reshape data
budworm_counts <- as.data.frame(rbind(budworm_raw[,1:4],budworm_raw[,5:8],budworm_raw[,9:12],budworm_raw[,13:16]))
names(budworm_counts) <- c('ddeg','total','stage','count')
budworm_counts %>% arrange(ddeg, stage)

budworm_table <- tidyr::pivot_wider(data = budworm_counts, names_from = stage, values_from = count, names_sort = TRUE, names_prefix = 'stage') %>% mutate(ddeg_cent = scale(ddeg, scale = FALSE))


budworm_individuals <- budworm_counts %>% tidyr::uncount(count) %>% mutate(ddeg_cent = scale(ddeg, scale = FALSE))
ggplot(budworm_counts, aes(x = ddeg, y = stage, size = count)) + geom_point()

#fit sratio model
budworm_sratio_probit <- brm(stage ~ ddeg_cent, data = budworm_individuals, family = sratio(link = 'probit'), inits = "0", chains = 2, cores = 2)
budworm_sratio_probit_cs <- brm(stage ~ cs(ddeg_cent), data = budworm_individuals, family = sratio(link = 'probit'), inits = "0", chains = 2, cores = 2)
#try to see if parameter estimates match Candy 1991 who seems to have fitted category specific effects
budworm_sratio_logit <- brm(stage ~ ddeg, data = budworm_individuals, family = sratio(link = 'logit'), inits = "0", chains = 2, cores = 2)
budworm_sratio_logit_cs <- brm(stage ~ cs(ddeg), data = budworm_individuals, family = sratio(link = 'logit'), inits = "0", chains = 2, cores = 2) # these parameter estimates broadly match the results in Candy et al. 1991
budworm_sratio_cloglog <- brm(stage ~ ddeg_cent, data = budworm_individuals, family = sratio(link = 'cloglog'), inits = "0", chains = 2, cores = 2)
budworm_sratio_cloglog_cs <- brm(stage ~ cs(ddeg), data = budworm_individuals, family = sratio(link = 'cloglog'), inits = "0", chains = 2, cores = 2)# again broadly in agreement with Candy et al. 1991

library(VGAM)
(budworm_sratio_logit_vglm <- vglm(cbind(stage1, stage2, stage3, stage4, stage5, stage6, stage7) ~ ddeg,
             sratio(link = "logitlink", parallel = FALSE), data = budworm_table))

t(coef(budworm_sratio_logit_vglm, matrix = TRUE))

(budworm_sratio_cloglog_vglm <- vglm(cbind(stage1, stage2, stage3, stage4, stage5, stage6, stage7) ~ ddeg_cent,
                                   sratio(link = "clogloglink", parallel = FALSE), data = budworm_table, trace=TRUE, control = vglm.control(checkwz = TRUE, wzepsilon = 0.5)))#sratio model with clogloglink fails because of some numerical issue - but check the expected binomial sample sizes. Candy writes "Note that the n_ij for which the binomial sample size is N^*_ij is zero must be expluded from the fit by by giving them a prior weight of zero. VGAM requires positive weights but suggests using weights "such as 1e-8" to effectively exclude observations"

#Nstar_ij = N_i - n_i1 - ... - n_i(j-1)

count_matrix <- as.matrix(budworm_table[,3:9])

Nstar <- matrix(NA, nrow = nrow(budworm_table), ncol = 7)
Nstar[,1] = rowSums(count_matrix)
Nstar[,2] = rowSums(count_matrix) - count_matrix[,1]
Nstar[,3] = rowSums(count_matrix) - rowSums(count_matrix[,1:2])
Nstar[,4] = rowSums(count_matrix) - rowSums(count_matrix[,1:3])
Nstar[,5] = rowSums(count_matrix) - rowSums(count_matrix[,1:4])
Nstar[,6] = rowSums(count_matrix) - rowSums(count_matrix[,1:5])
Nstar[,7] = rowSums(count_matrix) - rowSums(count_matrix[,1:6])

Nstar[Nstar == 0] <- 1e-8

(budworm_sratio_cloglog_vglm <- vglm(cbind(stage1, stage2, stage3, stage4, stage5, stage6, stage7) ~ ddeg,
                                     sratio(link = "clogloglink", parallel = FALSE), data = budworm_table, trace=TRUE, weights = Nstar[,2:7]))#adding weight doesn't seem to help, triggers other problems
?s
  