library(dplyr)
#Poisson deviance fitting approach
#
# read data ddeg day-degrees, total total individuals in stage,  
#! stage stage number 1-5 instars, 6 pupa, 7 adult
#! count total individuals in stage
budworm_counts <- readr::read_csv('data/budworm_counts.csv')
budworm_counts$occ <- as.numeric(as.factor(budworm_counts$ddeg)) #encode sampling occasion
budworm_counts$alpha_index = budworm_counts$stage + 1
budworm_counts <- group_by(budworm_counts, ddeg) %>% mutate(cumulative = cumsum(count))
#attach(budworm_counts) #TODO clean up handling of variable scope later
yvar <- count
stage <- factor(stage) #This might not be required if stage is only called in the boolean expressions
stage_ofac <- factor(stage, ordered = TRUE)
M2 <- function(){
  DR=rep(1, nrow(bud_budworm_counts))
}
M3 <- function(fitted_values){
  variance_function = fitted_values
  return(variance_function)
}
M4 <- function(){
  deviance_increment=2*(yvar*log(yvar/fitted_values)-(yvar-fitted_values)) 
  return(deviance_increment)
}  
MEXT <- function(){
  #extract parameter values #not sure how to translate that meaningfully
}
M1 <- function(parameters, CUL, CLT, X){
  #$calc %L=%ne(%PL,O) $swi %L MEXT $ if length(parameters!=0) #use macro MEXT?
  LP1=parameters[7]*X + #beta*z in eq2
  parameters[1]*(1+CUL)/X +#alpha1/sqrt(t)
  parameters[2]*(stage==2)/X +#alpha2/
  parameters[3]*(stage==3)/X + #
  parameters[4]*(stage==4)/X + # 
  parameters[5]*(stage==5)/X + #
  parameters[6]*(stage==6)/X 
  LP1=LP1*(LP1<20) + 21*(LP1>20) #why truncate LP1 to 21?
  
  LP2=parameters[7]*X 
  LP2=LP2 + parameters[1]*(1+CLT)/X + parameters[2]*(stage==3)/X 
  LP2=LP2 + parameters[3]*(stage==4)/X + parameters[4]*(stage==5)/X  
  LP2=LP2 + parameters[5]*(stage==6)/X + parameters[6]*(stage==7)/X 
  LP2=LP2*(LP2<20) + 20*(LP2>20) #why truncate to 20?
  
  #apply link function?
  F1=exp(LP1)/(1+exp(LP1)) 
  F2=exp(LP2)/(1+exp(LP2)) 
  
  #link function with squared denominator(why?)
  FD1=exp(LP1)/(1+exp(LP1))**2 
  FD2=exp(LP2)/(1+exp(LP2))**2 
  
  WC2=total*(FD1*(stage==2)-FD2*(stage==3))/X#why divide by X here?
  WC3=total*(FD1*(stage==3)-FD2*(stage==4))/X
  WC4=total*(FD1*(stage==4)-FD2*(stage==5))/X 
  WC5=total*(FD1*(stage==5)-FD2*(stage==6))/X 
  WC6=total*(FD1*(stage==6)-FD2*(stage==7))/X 
  WGM=total*(FD1*(1+CUL)-FD2*(1+CLT))/X 
  WB=total*(FD1-FD2)*X 
  
  linear_predictor=parameters[1]*WGM+parameters[2]*WC2+parameters[3]*WC3 
  linear_predictor=linear_predictor+parameters[4]*WC4+parameters[5]*WC5 
  linear_predictor=linear_predictor+parameters[6]*WC6+parameters[7]*WB 
  
  fitted_values=total*(F1-F2) 
  fitted_values=fitted_values*(fitted_values>0) + 0.0001*(fitted_values<0) 

}
#$own M1 M2 M3 M4 $
#! initial estimates of cut-point parameters alpha_i
# ! and regression parameter P
parameters= c(45,80,125,190,250,300,-0.6)
X=sqrt(ddeg)
#specify abitrarily small (CLT) and large (CUL)
#starting and finishing points for X
CLT=-(stage==1) #I am not clear if CLT and CUL are supposed to be scalars or vectors 
CUL=10*(stage==7)



# 
#      $use M1 $use M4 $
#      ! deviance using initial estimates
#    $calc %S=%cu(deviance_increment) $print %S $
#                       355.3
#                    $cycle 8 1 $
#                       ! fit the model
#                     $fit -%GM+WGM+WC2+WC3+WC4+WC5+WC6+WB $disp me $
#                       deviance = 79.995 at cycle 1
#                     deviance = 38.343 at cycle 2
#                     deviance = 34.883 at cycle 3
#                     deviance = 34.841 at cycle 4
#                     deviance = 34.841 at cycle 5
#                     d.f. = 77
#                    
#                     Current model:
#                      
#                     number of units is 84 
#                    
#                     y-variate count
#                     weight *
#                       offset *
#                      
#                     probability distribution is defined via the macros 
#                     MI, M2, M3 and M4
#                     scale parameter is to be estimated by the mean deviance 
#                     
#                     terms = WGM + WC2 + WC3 + WC4 + WC5 + WC6 + WB
#                    
#                         estimate  s.e.  parameter
#                     1 101.0 4.025 WGM
#                     2 71.22 4.770 WC2
#                     3 121.6 5.285 WC3
#                     4 186.2 6.677 WC4
#                     5 289.9 9.555 WC5
#                     6 400.3 13.35 WC6
#                     7 -0.8416 0.02689 WB
#                     scale parameter taken as 0.4525
#                    
#                     ! Note that GLIM's degrees of freedom for the deviance 
#  ! is incorrect since the constraint that the fitted 
#  ! values sum to.the sample size at each sampling
#  ! occasion, produced using CLT (=alpha_0) and CUL (=alpha_a), 
#  ! has not been taken into account. The correct degrees 
#  ! of freedom here is 84-7-12=65.
#  ! Also, if a multinomial error structure is assumed 
#  ! then the standard errors in the above table should be 
#  ! divided by sqrt(0.4525)
#  $return
#  ? $stop



#Candy model based on eqn 3 from scratch
alpha_index <- as.numeric(stage) + 1
alpha <- c(-10, 101,71.2, 121, 186.0, 289.9, 400.3, 5000)
alpha <- c(-10, 5.49,3.9,6.74, 10.21, 15.77, 21.76, 5000)
beta <- -0.8416
beta <- -0.0457
expected <- rep(NA, length(total))
for (i in 1:length(total)){
  LP1 <- alpha[alpha_index[i]] + beta*sqrt(ddeg[i])
  LP1 <- ifelse(LP1 > 20, 21, LP1)
  FP1 <- logitlink(LP1, inverse = TRUE)
  LP2 <- alpha[alpha_index[i]-1] + beta*sqrt(ddeg[i])
  LP2 <- ifelse(LP2 > 20, 20, LP2 )
  FP2 <- logitlink(LP2, inverse = TRUE)
  expected[i] <- total[i]*(FP1 - FP2)
}

count
plot(count ~ expected)
cbind(ddeg,stage,count,round(expected))
sum(count == 0)
sum(expected == 0)
head(budworm_counts)


pred_data <- expand.grid(ddeg = 0:800, stage = 1:7, total = 100)
pred_data$alpha_index = pred_data$stage + 1

##eqn 2
#alpha = c(-10, 120.0,	204.7-120,	264.6-120,	341.3-120,	464.5-120,	595.7-120,5000)
alpha = c(-5000, 101.0, 71.2, 121.7, 186.2,289.9,400.3,5000)
alpha = c(-5000, 101.0, 71.2+101.0, 121.7+101.0, 186.2+101.0,289.9+101.0,400.3+101.0,5000)
#beta = -0.8008975
beta = - 0.8416
expected_m_ij<- rep(NA, nrow(budworm_counts))
for (i in 1:nrow(budworm_counts)){
  expected_m_ij[i] <- with(budworm_counts,
       total[i]*exp(alpha[alpha_index[i]]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i]))/(1 + exp(alpha[alpha_index[i]]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i]))
       ))}
plot(expected_m_ij)
plot(budworm_counts$cumulative, expected_m_ij)
abline(0,1)

pred_m_ij<- rep(NA, nrow(pred_data))
for (i in 1:nrow(pred_data)){
  pred_m_ij[i] <- with(pred_data,
                           total[i]*exp(alpha[alpha_index[i]]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i]))/(1 + exp(alpha[alpha_index[i]]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i]))
                           ))}
plot(pred_m_ij)
plot(pred_data$ddeg,pred_m_ij, col = pred_data$stage)


expected_n_ij<- rep(NA, nrow(budworm_counts))
for (i in 1:nrow(budworm_counts)){
  expected_n_ij[i] <- with(budworm_counts,
                           total[i]*((exp(alpha[alpha_index[i]]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i]))/(1 + exp(alpha[alpha_index[i]]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i])))) - 
                                       (exp(alpha[alpha_index[i]-1]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i]))/(1 + exp(alpha[alpha_index[i]-1]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i])))))
                           )}
plot(expected_n_ij)
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

poisson_deviance_cm_candy(par = c(45,80,125,190,250,300,-0.6), data = budworm_counts)
poisson_nll_cm_candy(par = c(45,80,125,190,250,300,-0.6), data = budworm_counts)
logit_cm_candy <- optim(par = c(45,80,125,190,250,300,-0.6), poisson_deviance_cm_candy, data = budworm_counts, hessian = TRUE, control = list(trace=1), method = 'BFGS')
logit_cm_candy_nll <- optim(par = c(10,20,30,40,50,60,-0.6), poisson_nll_cm_candy, data = budworm_counts, hessian = TRUE, control = list(trace=1), method = 'BFGS')
logit_cm_candy_nll$par

cloglog_cm_candy_nll <- optim(par = c(5,8,12,19,25,30,-0.6), poisson_nll_cm_candy, data = budworm_counts, linkinv = function(x){VGAM::clogloglink(x, inverse = TRUE)}, hessian = TRUE, control = list(trace=1), method = 'BFGS')
cloglog_cm_candy_nll$par 

poisson_nll_cm_candy_dennis(par = c(45,80,125,190,250,300,-0.6), data = budworm_counts)
logit_dennis_cm_candy_nll <- optim(par = c(101.0, 71.2+101.0, 121.7+101.0, 186.2+101.0,289.9+101.0,400.3+101.0,-0.6), poisson_nll_cm_candy_dennis, data = budworm_counts, hessian = TRUE, control = list(trace=1), method = 'BFGS')
logit_dennis_cm_candy_nll$par

