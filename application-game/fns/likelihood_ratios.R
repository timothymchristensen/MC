#
# Calculate critical values for procedure 1
#
p1_crit_weighted <- function(pars_mle,mcmc_draws,mdata_short,ssame,weights,alpha_vec,cl){
  #
  n_alpha <- length(alpha_vec)
  B <- nrow(mcmc_draws)
  mcmc_LL <- numeric(B)
  #
  mle_LL <- -loglike_short(pars_mle,mdata_short,FALSE,NULL,NULL,ssame)
  #
  mcmc_LL <- -parRapply(cl,mcmc_draws,loglike_short,mdata_short=mdata_short,transformed=FALSE,lower_bounds=NULL,upper_bounds=NULL,ssame=ssame)
  #
  LL_diff <- abs(mle_LL - mcmc_LL)
  #
  ret <- wquantile(LL_diff,weights,alpha_vec)
  #
  return(ret)
}
#
# Weighted quantile
#
wquantile <- function(x,weights,p){
  #
  ix <- order(x)
  x_sorted <- x[ix]
  weights_sorted <- weights[ix]
  cum_weights <- (cumsum(weights/sum(weights)))
  #
  ret <- numeric(length(p))
  for(i in 1:length(p)){
    ret[i] <- x_sorted[sum(p[i]>cum_weights)+1]
  }
  #
  return(ret)
}
#
# Percentile CS
#
percentile_CS_weighted <- function(chain,weights,alpha_vec){
  #
  lower <- wquantile(chain,weights,(1-alpha_vec)/2)
  upper <- wquantile(chain,weights,1-(1-alpha_vec)/2)
  #
  return(list(lower=lower,upper=upper,length=upper-lower))
}
#
# Check length of CS for subvectors by profiling
#
subvec_length <- function(fixed_ix,pars_mle,mdata_short,lower_bounds,upper_bounds,ssame,cl,xi){
  #
  mle_LL <- -loglike_short(pars_mle,mdata_short,FALSE,NULL,NULL,ssame)
  #
  # first check on discrete grid
  mgrid <- seq(lower_bounds[fixed_ix],upper_bounds[fixed_ix],length.out=500)
  LL_mgrid <- parRapply(cl,matrix(mgrid),profile_LL,fixed_ix,pars_mle=pars_mle,mdata_short=mdata_short,lower_bounds=lower_bounds,upper_bounds=upper_bounds,ssame=ssame)
  #
  LL_diff <- (abs(mle_LL - LL_mgrid))
  #
  CS_length <- matrix(NA,3,length(xi))
  #
  for(j in 1:length(xi)){
    ix <- (LL_diff <= xi[j])
    upper_lim <- max(mgrid[ix])
    lower_lim <- min(mgrid[ix])
    #
    #
    CS_length[1,j] <- lower_lim
    CS_length[2,j] <- upper_lim
    CS_length[3,j] <- upper_lim - lower_lim
  }
  #
  return(CS_length)
}
#
# Calculate critical values for procedure 2
#
p2_crit_weighted <- function(fixed_ix,pars_mle,mcmc_draws,mdata_short,lower_bounds,upper_bounds,ssame,weights,alpha_vec,cl){
  #
  n_alpha <- length(alpha_vec)
  B <- nrow(mcmc_draws)
  profile_qlr <- numeric(B)
  #
  mle_LL <- -loglike_short(pars_mle,mdata_short,FALSE,NULL,NULL,ssame)
  #
  profile_qlr <- mle_LL-parRapply(cl,mcmc_draws,p2_crit_fn,fixed_ix=fixed_ix,pars_mle=pars_mle,mdata_short=mdata_short,lower_bounds=lower_bounds,upper_bounds=upper_bounds,ssame)
  if(any(is.nan(profile_qlr))){profile_qlr[is.nan(profile_qlr)]<--1000000}
  ret <- wquantile(profile_qlr,weights,alpha_vec)
  #
  return(ret)
}
#
# Inner part for profile QLR 
#
p2_crit_fn <- function(pars_temp,fixed_ix,pars_mle,mdata_short,lower_bounds,upper_bounds,ssame){
  #
  ptemp <- par_to_prob(pars_temp,ssame)
  initial_vals <- par_trans(pars_temp,lower_bounds,upper_bounds)[-fixed_ix]
  #
  lower_bd  <- find_int_lower(p0=ptemp,fixed_par=pars_temp[fixed_ix],fixed_ix,initial_vals,lower_bounds,upper_bounds,ssame,tol=1e-5)
  upper_bd  <- find_int_upper(p0=ptemp,fixed_par=pars_temp[fixed_ix],fixed_ix,initial_vals,lower_bounds,upper_bounds,ssame,tol=1e-5)
  #
  LL_lower <- profile_LL(lower_bd,fixed_ix,pars_mle,mdata_short,lower_bounds,upper_bounds,ssame)
  LL_upper <- profile_LL(upper_bd,fixed_ix,pars_mle,mdata_short,lower_bounds,upper_bounds,ssame)
  #
  ret <- min(LL_lower,LL_upper)
  #
  return(ret)
}
#
# END









