#
# Log-likelihood for trade model
#
loglike_trade <- function(pars,data,het=FALSE){
  Y1 <- data$Y1  # dep var for outcome equation
  X  <- data$X   # regressors for outcome equation
  Z0 <- data$Z0  # Z for selection equation (no trade)
  beta_m <- matrix(0,nrow=21,ncol=1)
  beta_m[-11,] <- pars[1:20]
  #
  beta_z <- matrix(pars[21:41])
  #
  delta <- pars[42]
  sigma_m <- exp(pars[43])
  rho <- (exp(pars[44])-1)/(1+exp(pars[44]))
  #
  if(het){
    varpi <- pars[45:46]
    sigma_z <- exp(varpi[1]*X[,2] + varpi[2]*X[,2]^2)
    sigma_z0 <- exp(varpi[1]*Z0[,2] + varpi[2]*Z0[,2]^2)
  }else{
    sigma_z <- sigma_z0 <- 1
  }
  sigma_v <- sqrt(sigma_m^2+2*rho*delta*sigma_m*sigma_z+(delta*sigma_z)^2)
  r_Sigma <- (rho*sigma_m*sigma_z+(sigma_z^2)*delta)/(sigma_v*sigma_z)
  #
  if(any(is.na(r_Sigma))==TRUE || any(is.nan(r_Sigma))==TRUE){
    return(1E10)
  }
  if(any(r_Sigma < -1)==TRUE || any(r_Sigma > 1)==TRUE){
    return(1E10)
  }
  #
  term_1 <- .Internal(pnorm(Z0%*%beta_z / sigma_z0,0,1,FALSE,TRUE))
  #
  term_2_1 <- -log(sigma_v)
  resid <- Y1-X%*%(beta_m+delta*beta_z)
  #
  term_2_2 <- -0.5*(resid/sigma_v)^2- log(sqrt(2*pi))
  #
  term_2_3 <- .Internal(pnorm( (X%*%beta_z + (sigma_z / sigma_v)*r_Sigma*resid) / (sqrt(1-r_Sigma^2) * sigma_z) ,0,1,TRUE,TRUE))
  #
  term_2 <- term_2_1 + term_2_2 + term_2_3
  #
  ret <- sum(term_1) + sum(term_2)
  #
  # return negative log likelihood
  if(is.na(ret)==TRUE||is.nan(ret)==TRUE||is.infinite(ret)==TRUE){
    return(1E10)
  }else{
    return(-ret)
  }
}
#
# ML estimates and conventional standard errors
#
mle <- function(initial_vals,data,het){
  #
  n_pars <- length(initial_vals)
  #
  opt <- optim(initial_vals,loglike_trade,data=data,het=het,method="Nelder-Mead",control=list(maxit=10000),hessian=FALSE)
  opt <- optim(opt$par,loglike_trade,data=data,het=het,method="BFGS",control=list(maxit=10000),hessian=FALSE)
  opt <- optim(opt$par,loglike_trade,data=data,het=het,method="SANN",control=list(maxit=10000),hessian=FALSE)
  opt <- optim(opt$par,loglike_trade,data=data,het=het,method="Nelder-Mead",control=list(maxit=10000),hessian=FALSE)
  opt <- optim(opt$par,loglike_trade,data=data,het=het,method="BFGS",control=list(maxit=10000),hessian=FALSE)
  opt <- optim(opt$par,loglike_trade,data=data,het=het,method="SANN",control=list(maxit=10000),hessian=FALSE)
  opt <- optim(opt$par,loglike_trade,data=data,het=het,method="BFGS",control=list(maxit=10000),hessian=TRUE)
  #
  mle <- opt$par
  cov <- solve(opt$hessian)
  ses <- sqrt(diag(cov))
  #
  return(list(mle=mle,cov=cov,ses=ses))
}
#
# ML estimates without standard errors
#
mle_short <- function(initial_vals,data,het){
  #
  n_pars <- length(initial_vals)
  #
  opt <- optim(initial_vals,loglike_trade,data=data,het=het,method="Nelder-Mead",control=list(maxit=10000),hessian=FALSE)
  opt <- optim(initial_vals,loglike_trade,data=data,het=het,method="BFGS",control=list(maxit=10000),hessian=FALSE)
  #
  mle <- opt$par
  #
  return(mle)
}
#
# profile MLE
#
profile_mle <- function(fixed_pars,fixed_ix,initial_vals,data,het){
  #
  initial_vals_profile <- initial_vals[-fixed_ix]
  #
  opt <- optim(initial_vals_profile,loglike_trade_profile,fixed_pars=fixed_pars,fixed_ix=fixed_ix,data=data,het=het,method="Nelder-Mead",control=list(maxit=10000),hessian=FALSE)
  opt <- optim(opt$par,loglike_trade_profile,fixed_pars=fixed_pars,fixed_ix=fixed_ix,data=data,het=het,method="BFGS",control=list(maxit=10000),hessian=FALSE)
  mle <- opt$par
  #
  return(mle)
}
#
# profile MLE objective function
#
loglike_trade_profile <- function(free_pars,fixed_pars,fixed_ix,data,het){
  #
  pars <- numeric(length(fixed_pars)+length(free_pars))
  pars[-fixed_ix] <- free_pars
  pars[fixed_ix] <- fixed_pars
  #
  ret <- loglike_trade(pars,data,het)
  #
  return(ret)
}
#
# profile LL
# 
profile_LL <- function(fixed_par,fixed_ix,initial_vals,data,het){
  #
  pars_profile_mle <- profile_mle(fixed_par,fixed_ix,initial_vals,data,het)
  #
  pars <- numeric(length(initial_vals))
  pars[-fixed_ix] <- pars_profile_mle
  pars[fixed_ix] <- fixed_par
  #
  ret <- -loglike_trade(pars,data,het)
  #
  return(ret)
}
#
# END