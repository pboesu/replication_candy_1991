P <- (stage==1)*1/(1+exp(-(A1-ddeg)/sqrt(BB*ddeg)))+
  (stage==2)*1/(1+exp(-(A2-ddeg)/sqrt(BB*ddeg)))-1/(1+exp(-(A1-ddeg)/sqrt(BB*ddeg)))+
  (stage==3)*1/(1+exp(-(A3-ddeg)/sqrt(BB*ddeg)))-1/(1+exp(-(A2-ddeg)/sqrt(BB*ddeg)))+
  (stage==4)*1/(1+exp(-(A4-ddeg)/sqrt(BB*ddeg)))-1/(1+exp(-(A3-ddeg)/sqrt(BB*ddeg)))+
  (stage==5)*1/(1+exp(-(A5-ddeg)/sqrt(BB*ddeg)))-1/(1+exp(-(A4-ddeg)/sqrt(BB*ddeg)))+
  (stage==6)*1/(1+exp(-(A6-ddeg)/sqrt(BB*ddeg)))-1/(1+exp(-(A5-ddeg)/sqrt(BB*ddeg)))+
  (stage==7)*1/(1+exp((A6-ddeg)/sqrt(BB*ddeg)))



nls(count ~ total*(stage==1)*1/(1+exp(-(A1-ddeg)/sqrt(BB*ddeg)))+
  (stage==2)*1/(1+exp(-(A2-ddeg)/sqrt(BB*ddeg)))-1/(1+exp(-(A1-ddeg)/sqrt(BB*ddeg)))+
  (stage==3)*1/(1+exp(-(A3-ddeg)/sqrt(BB*ddeg)))-1/(1+exp(-(A2-ddeg)/sqrt(BB*ddeg)))+
  (stage==4)*1/(1+exp(-(A4-ddeg)/sqrt(BB*ddeg)))-1/(1+exp(-(A3-ddeg)/sqrt(BB*ddeg)))+
  (stage==5)*1/(1+exp(-(A5-ddeg)/sqrt(BB*ddeg)))-1/(1+exp(-(A4-ddeg)/sqrt(BB*ddeg)))+
  (stage==6)*1/(1+exp(-(A6-ddeg)/sqrt(BB*ddeg)))-1/(1+exp(-(A5-ddeg)/sqrt(BB*ddeg)))+
  (stage==7)*1/(1+exp((A6-ddeg)/sqrt(BB*ddeg))), data = budworm_counts, start = list(A1 = 150, A2 = 230, A3 = 280, A4 = 330, A5 = 440, A6 = 580, BB = 2), lower = c(100,150,230,280,330,440,0.5), upper = c(600,600,600,600,600,600,3), trace = TRUE, algorithm = 'port', control = nls.control(minFactor = 1/4096, maxiter = 200)) 

#optim
# wfn<-function(par, Cfl, Tsoil, Pw){
#   K<-par[1]
#   a<-par[2]
#   b<-par[3]
#   c<-par[4]
#   z<-par[5]
#   m<-length(Pw)
#   res<-rep(NA,m)
#   for (i in 1:m){
#     tmp<-1
#     if (Pw[i] <= z) tmp<-exp(-a*Tsoil[i])
#     ff<-K*tmp*exp(-b*Pw[i])+c
#     res[i]<-Cfl[i]-ff
#   }
#   as.numeric(crossprod(res))
# }
#a1<-optim(st,wfn, control=list(trace=1), Cfl=Cfl, Tsoil=Tsoil, Pw=Pw)
weighted_residual <- function(par, count, total, stage, ddeg){
  A1 <- par[1]
  A2 <- par[2]
  A3 <- par[3]
  A4 <- par[4]
  A5 <- par[5]
  A6 <- par[6]
  BB <- par[7]
  pred <- weight <- rep(NA, length(total))
  for (i in seq_along(total)){
    pred[i] <- (stage[i]==1)*1/(1+exp(-(A1-ddeg[i])/sqrt(BB*ddeg[i])))+
    (stage[i]==2)*(1/(1+exp(-(A2-ddeg[i])/sqrt(BB*ddeg[i])))-1/(1+exp(-(A1-ddeg[i])/sqrt(BB*ddeg[i]))))+
    (stage[i]==3)*(1/(1+exp(-(A3-ddeg[i])/sqrt(BB*ddeg[i])))-1/(1+exp(-(A2-ddeg[i])/sqrt(BB*ddeg[i]))))+
    (stage[i]==4)*(1/(1+exp(-(A4-ddeg[i])/sqrt(BB*ddeg[i])))-1/(1+exp(-(A3-ddeg[i])/sqrt(BB*ddeg[i]))))+
    (stage[i]==5)*(1/(1+exp(-(A5-ddeg[i])/sqrt(BB*ddeg[i])))-1/(1+exp(-(A4-ddeg[i])/sqrt(BB*ddeg[i]))))+
    (stage[i]==6)*(1/(1+exp(-(A6-ddeg[i])/sqrt(BB*ddeg[i])))-1/(1+exp(-(A5-ddeg[i])/sqrt(BB*ddeg[i]))))+
    (stage[i]==7)*(1/(1+exp((A6-ddeg[i])/sqrt(BB*ddeg[i]))))
  }
  pred <- pmax(pred, 1e-7)
  res <- (count - total*pred)^2
  weight <- 1/(total*pred)
  return(sum(weight*res))
}
weighted_residual(par = c(A1 = 150, A2 = 230, A3 = 280, A4 = 330, A5 = 440, A6 = 580, BB = 2), count = budworm_counts$count, total = budworm_counts$total, stage = budworm_counts$stage, ddeg = budworm_counts$ddeg)
weighted_residual(par = c(A1 = 120, A2 = 204, A3 = 264, A4 = 342, A5 = 465, A6 = 599, BB = 1.5), count = budworm_counts$count, total = budworm_counts$total, stage = budworm_counts$stage, ddeg = budworm_counts$ddeg)

a1 <- optim(par = c(A1 = 150, A2 = 230, A3 = 280, A4 = 330, A5 = 440, A6 = 580, BB = 1), weighted_residual, count = budworm_counts$count, total = budworm_counts$total, stage = budworm_counts$stage, ddeg = budworm_counts$ddeg, hessian = TRUE, control = list(trace=1), method = 'BFGS')


#optimize likelihood directly rather than doing the IRLS approach
neg_log_lik <- function(par, count, total, stage, ddeg){
  A1 <- par[1]
  A2 <- par[2]
  A3 <- par[3]
  A4 <- par[4]
  A5 <- par[5]
  A6 <- par[6]
  BB <- par[7]
  pred <- weight <- rep(NA, length(total))
  for (i in seq_along(total)){
    pred[i] <- (stage[i]==1)*1/(1+exp(-(A1-ddeg[i])/sqrt(BB*ddeg[i])))+
      (stage[i]==2)*(1/(1+exp(-(A2-ddeg[i])/sqrt(BB*ddeg[i])))-1/(1+exp(-(A1-ddeg[i])/sqrt(BB*ddeg[i]))))+
      (stage[i]==3)*(1/(1+exp(-(A3-ddeg[i])/sqrt(BB*ddeg[i])))-1/(1+exp(-(A2-ddeg[i])/sqrt(BB*ddeg[i]))))+
      (stage[i]==4)*(1/(1+exp(-(A4-ddeg[i])/sqrt(BB*ddeg[i])))-1/(1+exp(-(A3-ddeg[i])/sqrt(BB*ddeg[i]))))+
      (stage[i]==5)*(1/(1+exp(-(A5-ddeg[i])/sqrt(BB*ddeg[i])))-1/(1+exp(-(A4-ddeg[i])/sqrt(BB*ddeg[i]))))+
      (stage[i]==6)*(1/(1+exp(-(A6-ddeg[i])/sqrt(BB*ddeg[i])))-1/(1+exp(-(A5-ddeg[i])/sqrt(BB*ddeg[i]))))+
      (stage[i]==7)*(1/(1+exp((A6-ddeg[i])/sqrt(BB*ddeg[i]))))
  }
  pred <- pmax(pred, 1e-8)
  nll <- -1*sum(count*log(pred))
  return(nll)
}
a2 <- optim(par = c(A1 = 150, A2 = 230, A3 = 280, A4 = 330, A5 = 440, A6 = 580, BB = 3), neg_log_lik, count = budworm_counts$count, total = budworm_counts$total, stage = budworm_counts$stage, ddeg = budworm_counts$ddeg, hessian = TRUE, control = list(trace=1), method = 'BFGS')


predicted_proportion <- function(par, stage, ddeg){
  A1 <- par[1]
  A2 <- par[2]
  A3 <- par[3]
  A4 <- par[4]
  A5 <- par[5]
  A6 <- par[6]
  BB <- par[7]
  pred <- weight <- rep(NA, length(stage))
  for (i in seq_along(stage)){
    pred[i] <- (stage[i]==1)*1/(1+exp(-(A1-ddeg[i])/sqrt(BB*ddeg[i])))+
      (stage[i]==2)*(1/(1+exp(-(A2-ddeg[i])/sqrt(BB*ddeg[i])))-1/(1+exp(-(A1-ddeg[i])/sqrt(BB*ddeg[i]))))+
      (stage[i]==3)*(1/(1+exp(-(A3-ddeg[i])/sqrt(BB*ddeg[i])))-1/(1+exp(-(A2-ddeg[i])/sqrt(BB*ddeg[i]))))+
      (stage[i]==4)*(1/(1+exp(-(A4-ddeg[i])/sqrt(BB*ddeg[i])))-1/(1+exp(-(A3-ddeg[i])/sqrt(BB*ddeg[i]))))+
      (stage[i]==5)*(1/(1+exp(-(A5-ddeg[i])/sqrt(BB*ddeg[i])))-1/(1+exp(-(A4-ddeg[i])/sqrt(BB*ddeg[i]))))+
      (stage[i]==6)*(1/(1+exp(-(A6-ddeg[i])/sqrt(BB*ddeg[i])))-1/(1+exp(-(A5-ddeg[i])/sqrt(BB*ddeg[i]))))+
      (stage[i]==7)*(1/(1+exp((A6-ddeg[i])/sqrt(BB*ddeg[i]))))
  }
  return(pred)}

#predictions
newdata <- expand.grid(stage = 1:7, ddeg = 50:700)
newdata$model_predictions <- predicted_proportion(a2$par, newdata$stage, newdata$ddeg)

ggplot(newdata, aes(x = ddeg, y = model_predictions, col = factor(stage), group = stage)) + geom_line() + geom_point(data = budworm_counts, aes(y = count/total)) + theme_classic()

ggplot(newdata, aes(x = ddeg, y = model_predictions, group = stage)) + geom_line() + geom_point(data = budworm_counts, aes(y = count/total)) + facet_wrap(~stage) + theme_classic()
