#fit sratio model
budworm_sratio_probit <- brm(stage ~ ddeg_cent, data = budworm_individuals, family = sratio(link = 'probit'), inits = "0", chains = 2, cores = 2)
budworm_sratio_probit_cs <- brm(stage ~ cs(ddeg_cent), data = budworm_individuals, family = sratio(link = 'probit'), inits = "0", chains = 2, cores = 2)
#try to see if parameter estimates match Candy 1991 who seems to have fitted category specific effects
budworm_sratio_logit <- brm(stage ~ ddeg, data = budworm_individuals, family = sratio(link = 'logit'), inits = "0", chains = 2, cores = 2)
budworm_sratio_logit_cs <- brm(stage ~ cs(ddeg), data = budworm_individuals, family = sratio(link = 'logit'), inits = "0", chains = 2, cores = 2) # these parameter estimates broadly match the results in Candy et al. 1991
budworm_sratio_cloglog <- brm(stage ~ ddeg_cent, data = budworm_individuals, family = sratio(link = 'cloglog'), inits = "0", chains = 2, cores = 2)
budworm_sratio_cloglog_cs <- brm(stage ~ cs(ddeg), data = budworm_individuals, family = sratio(link = 'cloglog'), inits = "0", chains = 2, cores = 2)# again broadly in agreement with Candy et al. 1991



#NB: VGAM and brms have conflicting namespaces!
budworm_cumulative_logit <- brm(stage ~ ddeg, data = budworm_individuals, family = brms::cumulative(link = 'logit'), chains = 2, cores = 2)#in good agreement with vgam estimates


cum_clog_stancode <- stancode(budworm_cumulative_cloglog)
cum_clog_standata <- standata(budworm_cumulative_cloglog)
initfun <- function(chain_id = 1) {
  list(
    b = as.array(0.04),
    Intercept = c(3,6,8,10,14,17)-20,
    b_Intercept = c(3,6,8,10,14,17),
    disc = as.array(1)
  )
}
budworm_cumulative_cloglog <- brm(stage ~ ddeg, data = budworm_individuals, family = brms::cumulative(link = 'cloglog'), inits = initfun, chains = 1, cores = 2)#fails on inits

#library(rstan)
#rstan_cum_clog <- stan(file = 'code/cum_clog.stan', data = cum_clog_standata, chains = 1, init = initfun) #despite providing initial values  stan tries to randomly generate some  -some how the issue might be linked to the initial value for the disc parameter or the b_Intercept
