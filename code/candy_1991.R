library(dplyr)
#Poisson deviance fitting approach
#
# read data ddeg day-degrees, total total individuals in stage,  
#! stage stage number 1-5 instars, 6 pupa, 7 adult
#! count total individuals in stage
budworm_counts <- readr::read_csv('data/budworm_counts.csv')
budworm_counts$occ <- as.numeric(as.factor(budworm_counts$ddeg)) #encode sampling occasion
budworm_counts$alpha_index = budworm_counts$stage + 1 #encode an index variable to facilitate fitting the linear predictor
budworm_counts <- group_by(budworm_counts, ddeg) %>% mutate(cumulative = cumsum(count)) #calculate cumulative counts across stages for cumulative likelihood


#optimize likelihood directly rather than doing the IRLS approach
poisson_deviance_cm_candy <- function(par, data, linkinv = gtools::inv.logit){
  alpha <- c(-10, par[1:6], 5000)
  beta = par[7]
  pred_n <- rep(NA, nrow(data))
  for (i in 1:nrow(data)){
    pred_n[i] <- with(data,
                      total[i] * (linkinv(alpha[alpha_index[i]] + beta*(ddeg[i])) - linkinv(alpha[alpha_index[i]-1] + beta*(ddeg[i])))
                      )
  }
  pred_n <- pmax(pred_n, 1e-10)
  pd <- ifelse(data$count != 0, 2*(data$count*log(data$count/pred_n) - (data$count - pred_n)), 0)
  return(sum(pd))
}

poisson_nll_cm_candy <- function(par, data, linkinv = gtools::inv.logit){
  alpha <- c(-10, par[1:6], 5000)
  beta = par[7]
  pred_n <- rep(NA, nrow(data))
  for (i in 1:nrow(data)){
    pred_n[i] <- with(data,
                      total[i] * (linkinv(alpha[alpha_index[i]] + beta*(ddeg[i])) - linkinv(alpha[alpha_index[i]-1] + beta*(ddeg[i])))
    )
  }
  #pred_n <- pmax(pred_n, 1e-10)
  #pd <- ifelse(data$count != 0, 2*(data$count*log(data$count/pred_n) - (data$count - pred_n)), 0)
  nll = sum(-1*dpois(data$count, pred_n, log=T))
  return(nll)
}

poisson_nll_cm_candy_dennis <- function(par, data){
  alpha <- c(-1, par[1:6], 5000)
  beta = par[7]
  pred_n <- rep(NA, nrow(data))
  for (i in 1:nrow(data)){
    pred_n[i] <- with(data,
                             total[i]*((exp(alpha[alpha_index[i]]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i]))/(1 + exp(alpha[alpha_index[i]]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i])))) - 
                                         (exp(alpha[alpha_index[i]-1]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i]))/(1 + exp(alpha[alpha_index[i]-1]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i])))))
    )}
  #pred_n <- pmax(pred_n, 1e-10)
  #pd <- ifelse(data$count != 0, 2*(data$count*log(data$count/pred_n) - (data$count - pred_n)), 0)
  nll = sum(-1*dpois(data$count, pred_n, log=T))
  return(nll)
}

logit_cm_candy <- optim(par = c(45,80,125,190,250,300,-0.6), poisson_deviance_cm_candy, data = budworm_counts, hessian = TRUE, control = list(trace=1), method = 'BFGS')
logit_cm_candy_nll <- optim(par = c(10,20,30,40,50,60,-0.6), poisson_nll_cm_candy, data = budworm_counts, hessian = TRUE, control = list(trace=1), method = 'BFGS')
#minimizing neg log likelihood leads to different parameter estimates than minimizing deviance
logit_cm_candy_nll$par
logit_cm_candy$par

cloglog_cm_candy_nll <- optim(par = c(5,8,12,19,25,30,-0.6), poisson_nll_cm_candy, data = budworm_counts, linkinv = function(x){VGAM::clogloglink(x, inverse = TRUE)}, hessian = TRUE, control = list(trace=1), method = 'BFGS')
cloglog_cm_candy_nll$par 

logit_dennis_cm_candy_nll <- optim(par = c(101.0, 71.2+101.0, 121.7+101.0, 186.2+101.0,289.9+101.0,400.3+101.0,-0.6), poisson_nll_cm_candy_dennis, data = budworm_counts, hessian = TRUE, control = list(trace=1), method = 'BFGS')
logit_dennis_cm_candy_nll$par
std_error_from_nll_hessian(logit_cm_candy_nll)

#create table of parameter estimates

#create predictions to compare with candy original fit
pred_data <- expand.grid(ddeg = 0:800, stage = 1:7, total = 100)
pred_data$alpha_index = pred_data$stage + 1

##eqn 2
alpha <- c(-10, logit_dennis_cm_candy_nll$par[1:6], 5000)
beta <- logit_dennis_cm_candy_nll$par[7]
pred_n_ij<- rep(NA, nrow(pred_data))
for (i in 1:nrow(pred_data)){
  pred_n_ij[i] <- with(pred_data,
                           total[i]*((exp(alpha[alpha_index[i]]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i]))/(1 + exp(alpha[alpha_index[i]]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i])))) - 
                                       (exp(alpha[alpha_index[i]-1]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i]))/(1 + exp(alpha[alpha_index[i]-1]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i])))))
  )}
plot(pred_data$ddeg,pred_n_ij)
plot(budworm_counts$count, expected_n_ij)
abline(0,1)



#plugging in paramater estimates from candy into his eqn 2 onnly makes sense when using alpha, not alpha* values


#eqn 3
#alpha <- c(-10, 5.49,3.9,6.74, 10.21, 15.77, 21.76, 5000)
alpha <- c(-10, 5.49,3.9+5.49,6.74+5.49 , 10.21+5.49, 15.77+5.49, 21.76+5.49, 5000)
beta = -0.0457
expected_n_ij<- rep(NA, nrow(budworm_counts))
linkinv <- gtools::inv.logit #function(x) exp(x)/(1+exp(x))#hand coded one is not numerically stable
for (i in 1:nrow(budworm_counts)){
  expected_n_ij[i] <- with(budworm_counts,
                           total[i] * (linkinv(alpha[alpha_index[i]] + beta*(ddeg[i])) - linkinv(alpha[alpha_index[i]-1] + beta*(ddeg[i]))))
}
plot(expected_n_ij, budworm_counts$count)
abline(0,1)

pred_n_ij<- rep(NA, nrow(pred_data))

for (i in 1:nrow(pred_data)){
  pred_n_ij[i] <- with(pred_data,
                       total[i] * (linkinv(alpha[alpha_index[i]] + beta*(ddeg[i])) - linkinv(alpha[alpha_index[i]-1] + beta*(ddeg[i]))))
}
plot(pred_data$ddeg, pred_n_ij/pred_data$total, col = pred_data$stage)

#eq1
expected_m_ij<- rep(NA, nrow(budworm_counts))
for (i in 1:nrow(budworm_counts)){
  expected_m_ij[i] <- budworm_counts$total[i]*linkinv(alpha[budworm_counts$alpha_index[i]] + beta*(budworm_counts$ddeg[i]))}
plot(expected_m_ij)
plot(budworm_counts$ddeg, expected_m_ij)
plot(budworm_counts$cumulative, expected_m_ij)

