budworm_counts <- readr::read_csv('data/budworm_counts.csv')

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


#predictions and comparisons
dennis_par = c(A1 = 121.08, A2 = 204.36, A3 = 264.41, A4 = 342.473, A5 = 465.620, A6 = 599.57, BB = 1.559)
sas_par = c(120.0, 204.7, 264.6, 341.3, 464.5, 595.7, 1.4119)

#figure 2
logistic_pdf <- function(s, t, bb){
  exp((s-t)/(sqrt(bb*t)))/(sqrt(bb*t)*(1+exp((s-t)/(sqrt(bb*t))))^2)
}
s <- 0:800
#redraw fig 2 using paramter estimates from Kemp
pdf('figures/dennis_fig2.pdf', width = 10, height = 7)
plot(s, logistic_pdf(s, 50, 1.559), type='l', ylab = "Probaility density", ylim = c(0,0.031), xaxs = 'i', yaxs='i')
for (t in c(150,250,350,450,550,650)){
  lines(s, logistic_pdf(s, t, 1.559))
}
abline(v = dennis_par[1:6], lty=2)
#overplot replication results
lines(s, logistic_pdf(s, 50, a2$par[7]), type='l', col = '#fc8d59')
for (t in c(150,250,350,450,550,650)){
  lines(s, logistic_pdf(s, t, a2$par[7]), col = "#fc8d59")
}
abline(v = a2$par[1:6], lty=2, col = "#fc8d59")
legend('topright',lty = 1, col = c('black','#fc8d59'), legend = c('Kemp et al. 1986', 'Replication'), bty = 'n')
text(dennis_par[1:6]+10, 0.030, latex2exp::TeX(c('a_1','a_2','a_3','a_4','a_5','a_6')))
text(c(50,150,250,350,450,550,650),c(0.0305,0.018,0.014,0.012,0.011,0.01,0.009), paste('t =',c(50,150,250,350,450,550,650)))
dev.off()
#redraw fig 3
pdf('figures/dennis_fig3.pdf', width = 10, height = 7)
plot(s,predicted_proportion(dennis_par,rep(1, length(s)),s), type = 'l', ylab = 'Proportion in stage', xlab = 'Time (degree-days)')
for (a in 2:7){
  lines(s,predicted_proportion(dennis_par,rep(a, length(s)),s))
}
lines(s,predicted_proportion(a2$par,rep(1, length(s)),s), col = "#fc8d59")
for (a in 2:7){
  lines(s,predicted_proportion(a2$par,rep(a, length(s)),s), col = "#fc8d59")
}
legend('top',lty = 1, col = c('black','#fc8d59'), legend = c('Kemp et al. 1986', 'Replication'), bty = 'n')
text(c(90,160,233,300,400,530,690), c(0.98,0.91,0.7,0.75,0.89,0.86,0.98), 1:7)
dev.off()