#
# MLE for Example 2
#
mle <- function(initial_vals,mdata_short,lower_bounds=NULL,upper_bounds=NULL){
    #
    initial_vals <- par_trans(initial_vals,lower_bounds,upper_bounds)
    #
    opt <- nlm(loglike_short_true,p=initial_vals,mdata_short=mdata_short,lower_bounds=lower_bounds,upper_bounds=upper_bounds)
    opt <- nlm(loglike_short_true,p=opt$estimate,mdata_short=mdata_short,lower_bounds=lower_bounds,upper_bounds=upper_bounds)
    #
    MLE <- par_inv_trans(opt$estimate,lower_bounds,upper_bounds)
    #
    return(MLE)
}
#
# Short log likelihood for example 2
#
loglike_short <- function(pars,mdata_short,transformed=TRUE,lower_bounds=NULL,upper_bounds=NULL){
  #
  if(transformed==TRUE){
    pars <- (lower_bounds + upper_bounds*exp(pars))/(1+exp(pars))
    if(any(is.nan(pars))){pars[is.nan(pars)]<-upper_bounds[is.nan(pars)]}
  }
  #
  pi_vec <- par_to_prob(pars)
  #
  ret <- sum(mdata_short[1:4]*log(pi_vec))
  if(is.nan(ret)){
    ret <- -1000000
  }
  # return MINUS the log-likelihood
  return(-ret)
}
#
# Short log likelihood for example 2
#
loglike_short_true <- function(pars,mdata_short,lower_bounds=NULL,upper_bounds=NULL){
  #
  pars <- (lower_bounds + upper_bounds*exp(pars))/(1+exp(pars))
  if(any(is.nan(pars))){pars[is.nan(pars)]<-upper_bounds[is.nan(pars)]}
  if(any(is.infinite(pars))){pars[is.infinite(pars)]<-upper_bounds[is.infinite(pars)]}
  #
  pi_vec <- par_to_prob(pars)
  #
  ret <- sum(mdata_short[1:4]*log(pi_vec))
  if(is.nan(ret)){
    ret <- -1000000
  }
  # return MINUS the log-likelihood
  return(-ret)
}
#
# MLE with Delta_1 fixed
#
mle_Delta_1_fixed <- function(initial_vals,Delta_1,mdata_short,lower_bounds=NULL,upper_bounds=NULL){
  #
  Delta_1 <- par_trans(Delta_1,lower_bounds[1],upper_bounds[1])
  initial_vals <- par_trans(initial_vals,lower_bounds[-(1)],upper_bounds[-(1)])
  #
  opt <- nlm(loglike_Delta_1_fixed,p=initial_vals,Delta_1=Delta_1,mdata_short=mdata_short,lower_bounds=lower_bounds,upper_bounds=upper_bounds)
  #
  MLE <- par_inv_trans(opt$estimate,lower_bounds[-(1)],upper_bounds[-(1)])
  #
  return(MLE)
}
#
# LL with Delta_1 fixed
#
loglike_Delta_1_fixed <- function(pars,Delta_1,mdata_short,lower_bounds=NULL,upper_bounds=NULL){
  #
  pars <- c(Delta_1,pars)
  ret <- loglike_short_true(pars,mdata_short,lower_bounds,upper_bounds)
  #
  return(ret)
}
#
# Profile LL for Delta_1
# 
profile_LL_Delta_1 <- function(Delta_1,pars_mle,mdata_short,lower_bounds=NULL,upper_bounds=NULL){
  #
  pars_hat_mu <- mle_Delta_1_fixed(pars_mle[-1],Delta_1,mdata_short,lower_bounds,upper_bounds)
  #
  ret <- -loglike_short(c(Delta_1,pars_hat_mu),mdata_short,FALSE)
  #
  return(ret)
}
#
# MLE with beta_1 fixed
#
mle_beta_1_fixed <- function(initial_vals,beta_1,mdata_short,lower_bounds=NULL,upper_bounds=NULL){
  #
  beta_1 <- par_trans(beta_1,lower_bounds[3],upper_bounds[3])
  initial_vals <- par_trans(initial_vals,lower_bounds[-(3)],upper_bounds[-(3)])
  #
  opt <- nlm(loglike_beta_1_fixed,p=initial_vals,beta_1=beta_1,mdata_short=mdata_short,lower_bounds=lower_bounds,upper_bounds=upper_bounds)
  #
  MLE <- par_inv_trans(opt$estimate,lower_bounds[-(3)],upper_bounds[-(3)])
  
  #
  return(MLE)
}
#
# LL with beta_1 fixed
#
loglike_beta_1_fixed <- function(pars,beta_1,mdata_short,lower_bounds=NULL,upper_bounds=NULL){
  #
  pars_full <- numeric(6)
  pars_full[3] <- beta_1
  pars_full[-3] <- pars
  ret <- loglike_short_true(pars_full,mdata_short,lower_bounds,upper_bounds)
  #
  return(ret)
}
#
# Profile LL for beta_1
# 
profile_LL_beta_1 <- function(beta_1,pars_mle,mdata_short,lower_bounds=NULL,upper_bounds=NULL){
  #
  pars_hat_mu <- mle_beta_1_fixed(pars_mle[-3],beta_1,mdata_short,lower_bounds,upper_bounds)
  #
  pars_full <- numeric(6)
  pars_full[3] <- beta_1
  pars_full[-3] <- pars_hat_mu
  ret <- -loglike_short(pars_full,mdata_short,FALSE)
  #
  return(ret)
}
#
#END