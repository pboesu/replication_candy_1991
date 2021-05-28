# function definitions

# likelihood function for cumulative model with constant variance
poisson_nll_cm_candy <- function(par, data, linkinv = gtools::inv.logit){
  alpha <- c(-10, par[1:6], 5000)
  beta = par[7]
  pred_n <- rep(NA, nrow(data))
  for (i in 1:nrow(data)){
    pred_n[i] <- with(data,
                      total[i] * (linkinv(alpha[alpha_index[i]] + beta*(ddeg[i])) - linkinv(alpha[alpha_index[i]-1] + beta*(ddeg[i])))
    )
  }
  # avoid numerical underflow during fitting
  pred_n <- pmax(pred_n, .Machine$double.xmin)
  nll = sum(-1*dpois(data$count, pred_n, log=T))
  return(nll)
}

# likelihood function for cumulative model with proportional variance
poisson_nll_cm_candy_dennis <- function(par, data){
  alpha <- c(-10, par[1:6], 5000)
  beta = par[7]
  pred_n <- rep(NA, nrow(data))
  for (i in 1:nrow(data)){
    pred_n[i] <- with(data,
                      total[i]*((exp(alpha[alpha_index[i]]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i]))/(1 + exp(alpha[alpha_index[i]]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i])))) - 
                                  (exp(alpha[alpha_index[i]-1]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i]))/(1 + exp(alpha[alpha_index[i]-1]/sqrt(ddeg[i]) + beta*sqrt(ddeg[i])))))
    )}
  # avoid numerical underflow during fitting
  pred_n <- pmax(pred_n, .Machine$double.xmin)
  nll = -1*sum(dpois(data$count, pred_n, log=T))
  return(nll)
}

# likelihood function for the cumulative model with constant variance without numerical thresholds against underflow 
poisson_nll_cm_candy_inf <- function(par, data, linkinv = gtools::inv.logit){
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

# function to repeatedly fit the cumulative logit model with randomised starting values
cm_cand_sens <- function(iter, progress = TRUE, linkinv = gtools::inv.logit){
  if (progress == TRUE & (iter %% 10 == 0)) print(paste(Sys.time(), '- iteration', iter))
  
  # draw random initial values
  alpha_start <- cumsum(runif(6,1,20))
  beta_start <- runif(1, -1,1)
  # estimate model parameters
  cm_candy_nll <- optim(par = c(alpha_start,beta_start), poisson_nll_cm_candy, data = budworm_counts, hessian = TRUE, control = list(trace=0, parscale = c(rep(1,6), 0.01)), method = 'BFGS', linkinv = linkinv)
  # check that converged parameters actually have a finite likelihood
  if (is.finite(poisson_nll_cm_candy_inf(par = cm_candy_nll$par, data = budworm_counts, linkinv = linkinv))){
    sens_out <- tibble::tibble(init = c(alpha_start,beta_start),
                       estimate = cm_candy_nll$par,
                       parameter = c(paste('a', 1:6, sep = ''), 'b'),
                       convergence = cm_candy_nll$convergence,
                       iter = iter,
                       value = cm_candy_nll$value) 
  } else { # else return an indicator of false convergence
    sens_out <- tibble::tibble(init = c(alpha_start,beta_start),
                       convergence = -1,
                       parameter = c(paste('a', 1:6, sep = ''), 'b'),
                       iter = iter)
  }
  return(sens_out)
}

# function to draw the likelihood surface for the cumulative model with constant variance 
cm_nll_surface <- function(iter, progress = TRUE, linkinv = gtools::inv.logit){
  if (progress == TRUE & (iter %% 1000 == 0)) print(paste(Sys.time(), '- iteration', iter))
  # draw random initial values
  alpha_start <- cumsum(c(runif(1,0,9),runif(5,0.1,7)))
  beta_start <- runif(1, -0.1,0)
  # estimate model parameters
  nll = poisson_nll_cm_candy(par = c(alpha_start,beta_start), data = budworm_counts, linkinv = linkinv)
  # check that converged parameters actually have a finite likelihood
  sens_out <- tibble::tibble(init = c(alpha_start,beta_start),
                             parameter = c(paste('a', 1:6, sep = ''), 'b'),
                             iter = iter,
                             value = nll)
  return(sens_out)
}

# function to repeatedly fit the VGAM cumulative logit model with randomised starting values
cm_vgam_sens <- function(iter, progress = TRUE, link = clogloglink()){
  if (progress == TRUE & (iter %% 200 == 0)) print(paste(Sys.time(), '- iteration', iter))
  
  # draw random initial values
  alpha_start <- cumsum(runif(6,0.1,5))
  beta_start <- runif(1, -0.1,-0.01)
  # estimate model parameters
  cm_vgam_nll <- try(vglm(cbind(stage1, stage2, stage3, stage4, stage5, stage6, stage7) ~ ddeg_cent, cumulative(reverse = FALSE, link = clogloglink(), parallel = TRUE), coefstart = c(alpha_start, beta_start), data = budworm_table, control = vglm.control(checkwz = TRUE)))#
  # check that model didn't error out
  if (!inherits(cm_vgam_nll, 'try-error')){
    sens_out <- tibble::tibble(init = c(alpha_start,beta_start),
                       estimate = coef(cm_vgam_nll),
                       parameter = c(paste('a', 1:6, sep = ''), 'b'),
                       convergence = 0,
                       iter = iter,
                       value = logLik(cm_vgam_nll)) 
  } else { # else return an indicator of false convergence
    sens_out <- tibble::tibble(init = c(alpha_start,beta_start),
                       convergence = -1,
                       parameter = c(paste('a', 1:6, sep = ''), 'b'),
                       iter = iter)
  }
  return(sens_out)
}

# likelihood definition for cumulative model with proportional variance. This is a quite literal translation of the original SAS code.
nll_cm_dennis <- function(par, count, total, stage, ddeg){
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

# alternative implementation of the cumulative model with proportional variance that uses the R function stats::plogis 
nll_cm_dennis_plogis <- function(par, count, total, stage, ddeg){
  A1 <- par[1]
  A2 <- par[2]
  A3 <- par[3]
  A4 <- par[4]
  A5 <- par[5]
  A6 <- par[6]
  BB <- par[7]
  pred <- weight <- rep(NA, length(total))
  for (i in seq_along(total)){
    pred[i] <- (stage[i]==1)*stats::plogis(A1, location = ddeg[i], scale = sqrt(BB*ddeg[i]))+
      (stage[i]==2)*(stats::plogis(A2, location = ddeg[i], scale = sqrt(BB*ddeg[i]))-stats::plogis(A1, location = ddeg[i], scale = sqrt(BB*ddeg[i])))+
      (stage[i]==3)*(stats::plogis(A3, location = ddeg[i], scale = sqrt(BB*ddeg[i]))-stats::plogis(A2, location = ddeg[i], scale = sqrt(BB*ddeg[i])))+
      (stage[i]==4)*(stats::plogis(A4, location = ddeg[i], scale = sqrt(BB*ddeg[i]))-stats::plogis(A3, location = ddeg[i], scale = sqrt(BB*ddeg[i])))+
      (stage[i]==5)*(stats::plogis(A5, location = ddeg[i], scale = sqrt(BB*ddeg[i]))-stats::plogis(A4, location = ddeg[i], scale = sqrt(BB*ddeg[i])))+
      (stage[i]==6)*(stats::plogis(A6, location = ddeg[i], scale = sqrt(BB*ddeg[i]))-stats::plogis(A5, location = ddeg[i], scale = sqrt(BB*ddeg[i])))+
      (stage[i]==7)*(stats::plogis(A6, location = ddeg[i], scale = sqrt(BB*ddeg[i]), lower.tail = FALSE))
  }
  pred <- pmax(pred, 1e-8)
  nll <- -1*sum(count*log(pred))
  return(nll)
}

# predictions for the cumulative model with proportional variance
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

# definition of logistic PDF
logistic_pdf <- function(s, t, bb){
  exp((s-t)/(sqrt(bb*t)))/(sqrt(bb*t)*(1+exp((s-t)/(sqrt(bb*t))))^2)
}
