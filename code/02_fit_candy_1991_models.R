library(dplyr)
#Poisson deviance fitting approach
#
# read data ddeg day-degrees, total total individuals in stage,  
#! stage stage number 1-5 instars, 6 pupa, 7 adult
#! count total individuals in stage
budworm_counts <- readr::read_csv('data/budworm_counts.csv')
budworm_counts$occ <- as.numeric(as.factor(budworm_counts$ddeg)) #encode sampling occasion
budworm_counts$alpha_index = budworm_counts$stage + 1 #encode an index variable to facilitate fitting the linear predictor
budworm_counts <- group_by(budworm_counts, ddeg) %>%
  arrange(ddeg,stage) %>% 
  mutate(cumulative = cumsum(count)) %>% #calculate cumulative counts across stages for cumulative likelihood
  mutate(n_star = ifelse(stage ==1, total, total - lag(cumulative)), #calculate N*_ij for GLM estimation of sequential model
         stage_factor = as.factor(stage)) 

#likelihood function for cumulative model with constant variance
poisson_nll_cm_candy <- function(par, data, linkinv = gtools::inv.logit){
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

#likelihood function for cumulative model with proportional variance
poisson_nll_cm_candy_dennis <- function(par, data){
  alpha <- c(-10, par[1:6], 5000)
  beta = par[7]
  pred_n <- rep(NA, nrow(data))
  for (i in 1:nrow(data)){
    pred_n[i] <- with(data,
                             total[i]*((exp(alpha[alpha_index[i]]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i]))/(1 + exp(alpha[alpha_index[i]]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i])))) - 
                                         (exp(alpha[alpha_index[i]-1]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i]))/(1 + exp(alpha[alpha_index[i]-1]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i])))))
    )}
  nll = -1*sum(dpois(data$count, pred_n, log=T))
  return(nll)
}

#get parameter estimates for cumulative models
logit_cm_candy_nll <- optim(par = c(10,20,30,40,50,60,-0.6), poisson_nll_cm_candy, data = budworm_counts, hessian = TRUE, control = list(trace=1), method = 'BFGS')
#minimizing neg log likelihood leads to different parameter estimates than minimizing deviance
logit_cm_candy_nll$par
#residual deviance
-2*(-1*logit_cm_candy_nll$value-sum(dpois(budworm_counts$count,budworm_counts$count, log = TRUE)))

cloglog_cm_candy_nll <- optim(par = c(5,8,12,19,25,30,-0.6), poisson_nll_cm_candy, data = budworm_counts, linkinv = function(x){VGAM::clogloglink(x, inverse = TRUE)}, hessian = TRUE, control = list(trace=1), method = 'BFGS')
cloglog_cm_candy_nll$par 
#residual deviance
-2*(-1*cloglog_cm_candy_nll$value-sum(dpois(budworm_counts$count,budworm_counts$count, log = TRUE)))

logit_dennis_cm_candy_nll <- optim(par = c(101.0, 71.2+101.0, 121.7+101.0, 186.2+101.0,289.9+101.0,400.3+101.0,-0.6), poisson_nll_cm_candy_dennis, data = budworm_counts, hessian = TRUE, control = list(trace=1), method = 'BFGS')
logit_dennis_cm_candy_nll$par
#std_error_from_nll_hessian(logit_cm_candy_nll)

#save parameter estimates
candy_nll_estimates <- as.data.frame(rbind(logit_cm_candy_nll$par, cloglog_cm_candy_nll$par, logit_dennis_cm_candy_nll$par))
candy_nll_estimates$model <- "cumulative"
candy_nll_estimates$link <- c("logit","cloglog","logit")
candy_nll_estimates$fit <- c("R \\verb+optim+")
candy_nll_estimates$eqn <- c("\\ref{eq:candy_cm_count_form}","\\ref{eq:candy_cm_count_form}","\\ref{eq:candy_cm_count_form}")
saveRDS(candy_nll_estimates, "outputs/candy_nll_estimates.RDS")

#get parameter estimates for sequential models
logit_sm_candy_glm <- glm(cbind(count, n_star-count) ~ stage_factor:ddeg + stage_factor-1, family = binomial(link="logit"), data = budworm_counts, subset = stage < 7)
cloglog_sm_candy_glm <- glm(cbind(count, n_star-count) ~ stage_factor:ddeg + stage_factor-1, family = binomial(link="cloglog"), data = budworm_counts, subset = stage < 7)


#save parameter estimates
candy_sm_nll_estimates <- as.data.frame(rbind(unname(coef(logit_sm_candy_glm)[1:6]),unname(coef(logit_sm_candy_glm)[7:12]), unname(coef(cloglog_sm_candy_glm)[1:6]), unname(coef(cloglog_sm_candy_glm)[7:12])))
candy_sm_nll_estimates$model <- "sequential"
candy_sm_nll_estimates$link <- c("logit","logit","cloglog","cloglog")
candy_sm_nll_estimates$fit <- c("R \\verb+glm+")
candy_sm_nll_estimates$par <- c("$\\beta_{0j}$","$\\beta_{1j}$","$\\beta_{0j}$","$\\beta_{1j}$")
saveRDS(candy_sm_nll_estimates, "outputs/candy_sm_nll_estimates.RDS")
