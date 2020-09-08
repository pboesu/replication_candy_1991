#helper functions
add_se_vcov_nll_hessian <- function(fit){
  fit$vcov <- solve(fit$hessian) # var-cov matrix
  fit$se <- sqrt(diag(fit$vcov))     # standard errors
  return(fit)
}
