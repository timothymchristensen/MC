#
# MLE for empirical game
#
mle <- function(initial_vals,mdata_short,lower_bounds=NULL,upper_bounds=NULL,ssame){
    #
    initial_vals <- par_trans(initial_vals,lower_bounds,upper_bounds)
    #
    opt <- nlm(loglike_short_true,p=initial_vals,mdata_short=mdata_short,lower_bounds=lower_bounds,upper_bounds=upper_bounds,ssame=ssame)
    opt <- nlm(loglike_short_true,p=opt$estimate,mdata_short=mdata_short,lower_bounds=lower_bounds,upper_bounds=upper_bounds,ssame=ssame)
    #
    MLE <- par_inv_trans(opt$estimate,lower_bounds,upper_bounds)
    #
    return(MLE)
}
#
# Short log likelihood for empirical game
#
loglike_short <- function(pars,mdata_short,transformed=TRUE,lower_bounds=NULL,upper_bounds=NULL,ssame){
  #
  if(transformed==TRUE){
    pars <- (lower_bounds + upper_bounds*exp(pars))/(1+exp(pars))
    if(any(is.nan(pars))){pars[is.nan(pars)]<-upper_bounds[is.nan(pars)]}
    if(any(is.infinite(pars))){pars[is.infinite(pars)]<-upper_bounds[is.infinite(pars)]}
  }
  #
  pi_vec <- par_to_prob(pars,ssame)
  #
  ret <- sum(mdata_short[mdata_short>0]*log(pi_vec[mdata_short>0]))
  if(is.nan(ret)){
    ret <- -10000000
  }
  # return MINUS the log-likelihood
  return(-ret)
}
#
# Short log likelihood for empirical game
#
loglike_short_true <- function(pars,mdata_short,lower_bounds=NULL,upper_bounds=NULL,ssame){
  #
  pars <- (lower_bounds + upper_bounds*exp(pars))/(1+exp(pars))
  if(any(is.nan(pars))){pars[is.nan(pars)]<-upper_bounds[is.nan(pars)]}
  if(any(is.infinite(pars))){pars[is.infinite(pars)]<-upper_bounds[is.infinite(pars)]}
  #
  pi_vec <- par_to_prob(pars,ssame)
  #
  ret <- sum(mdata_short[mdata_short>0]*log(pi_vec[mdata_short>0]))
  if(is.nan(ret)|is.infinite(ret)){
    ret <- -10000000
  }
  # return MINUS the log-likelihood
  return(-ret)
}
#
# MLE with subvector fixed
#
mle_par_fixed <- function(initial_vals,par_fixed,fixed_ix,mdata_short,lower_bounds=NULL,upper_bounds=NULL,ssame){
  #
  par_fixed <- par_trans(par_fixed,lower_bounds[fixed_ix],upper_bounds[fixed_ix])
  initial_vals <- par_trans(initial_vals,lower_bounds[-fixed_ix],upper_bounds[-fixed_ix])
  #
  opt <- nlm(loglike_par_fixed,p=initial_vals,par_fixed=par_fixed,fixed_ix=fixed_ix,mdata_short=mdata_short,lower_bounds=lower_bounds,upper_bounds=upper_bounds,ssame=ssame)
  opt <- nlm(loglike_par_fixed,p=opt$estimate,par_fixed=par_fixed,fixed_ix=fixed_ix,mdata_short=mdata_short,lower_bounds=lower_bounds,upper_bounds=upper_bounds,ssame=ssame)
  #
  MLE <- par_inv_trans(opt$estimate,lower_bounds[-fixed_ix],upper_bounds[-fixed_ix])
  #
  return(MLE)
}
#
# LL with subvector fixed
#
loglike_par_fixed <- function(pars,par_fixed,fixed_ix,mdata_short,lower_bounds=NULL,upper_bounds=NULL,ssame){
  #
  pars_full <- numeric(length(par_fixed)+length(pars))
  pars_full[fixed_ix] <- par_fixed
  pars_full[-fixed_ix] <- pars
  ret <- loglike_short_true(pars_full,mdata_short,lower_bounds,upper_bounds,ssame)
  #
  return(ret)
}
#
# Profile LL for subvector
# 
profile_LL <- function(par_fixed,fixed_ix,pars_mle,mdata_short,lower_bounds=NULL,upper_bounds=NULL,ssame){
  #
  pars_hat_mu <- mle_par_fixed(pars_mle[-fixed_ix],par_fixed,fixed_ix,mdata_short,lower_bounds,upper_bounds,ssame)
  #
  pars_full <- numeric(length(pars_mle))
  pars_full[fixed_ix] <- par_fixed
  pars_full[-fixed_ix] <- pars_hat_mu
  #
  ret <- -loglike_short(pars_full,mdata_short,FALSE,NULL,NULL,ssame)
  #
  return(ret)
}
#
#END