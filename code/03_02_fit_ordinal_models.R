#NB: package ordinal fits this without problems.... 
library(ordinal)
pm1 <- clm(factor(stage) ~ ddeg, data = budworm_individuals, link = 'cloglog')
summary(pm1)
#try to calculate the poisson deviance for this model