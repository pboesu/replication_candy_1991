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

budworm_cumulative_cloglog_vglm <- vglm(cbind(stage1, stage2, stage3, stage4, stage5, stage6, stage7) ~ ddeg,
                                        cumulative(reverse = TRUE, link = clogloglink(bvalue = 1e-7), parallel = TRUE), coefstart = c(-3,-6,-8,-10,-14,-17,0.01), data = budworm_table, trace=T)


budworm_cumulative_logit_vglm <- vglm(cbind(stage1, stage2, stage3, stage4, stage5, stage6, stage7) ~ ddeg,
                                      cumulative(reverse = TRUE, link = "logitlink", parallel = TRUE), data = budworm_table, trace=T) #reproduces parameter estimates in table 2 of candy
