#
# Calculate critical values for procedure 1
#
p1_crit_weighted <- function(pars_mle,mcmc_draws,mdata_short,weights,alpha_vec,cl){
  #
  n_alpha <- length(alpha_vec)
  B <- nrow(mcmc_draws)
  mcmc_LL <- numeric(B)
  #
  mle_LL <- -loglike_short(pars_mle,mdata_short,FALSE)
  #
  mcmc_LL <- -parRapply(cl,mcmc_draws,loglike_short,mdata_short=mdata_short,transformed=FALSE)
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
# Check coverage for the full vector identified set
#
full_check <- function(pars_mle,pars_0,mdata_short,xi){
  #
  mle_LL <- -loglike_short(pars_mle,mdata_short,FALSE)
  mle_0  <- -loglike_short(pars_0,mdata_short,FALSE)
  #
  LL_diff <- abs(mle_LL - mle_0)
  #
  ret <-  (LL_diff <= xi)*1
  #
  return(ret)
}
#
# Check coverage for the subvector identified set for Delta_1
#
subvec_Delta_1_cvg <- function(pars_mle,M_I,mdata_short,lower_bounds,upper_bounds,cl,xi){
  #
  mle_LL <- -loglike_short(pars_mle,mdata_short,FALSE,lower_bounds,upper_bounds)
  #
  mle_M_I <- parRapply(cl,matrix(M_I),profile_LL_Delta_1,pars_mle=pars_mle,mdata_short=mdata_short,lower_bounds=lower_bounds,upper_bounds=upper_bounds)
  #
  LL_diff <- max(abs(mle_LL - mle_M_I))
  #
  ret <-  (LL_diff <= xi)*1
  #
  return(ret)
}
#
# Check coverage for the subvector identified set for beta_1
#
subvec_beta_1_cvg <- function(pars_mle,M_I,mdata_short,lower_bounds,upper_bounds,cl,xi){
  #
  mle_LL <- -loglike_short(pars_mle,mdata_short,FALSE,lower_bounds,upper_bounds)
  #
  mle_M_I <- parRapply(cl,matrix(M_I),profile_LL_beta_1,pars_mle=pars_mle,mdata_short=mdata_short,lower_bounds=lower_bounds,upper_bounds=upper_bounds)
  #
  LL_diff <- max(abs(mle_LL - mle_M_I))
  #
  ret <-  (LL_diff <= xi)*1
  #
  return(ret)
}
#
# Percentile CS
#
percentile_CS_weighted <- function(chain,weights,alpha_vec,M_I=NULL){
  #
  lower <- wquantile(chain,weights,(1-alpha_vec)/2)
  upper <- wquantile(chain,weights,1-(1-alpha_vec)/2)
  #
  if(is.null(M_I)==TRUE){
    return(list(lower=lower,upper=upper,length=upper-lower))
  }else{
    cvg_check <- numeric(length(alpha_vec))
    for(i in 1:length(alpha_vec)){
      cvg_check[i] <- (lower[i]<=min(M_I)&&upper[i]>=max(M_I))*1
    }
    #
    return(list(lower=lower,upper=upper,length=upper-lower,cvg_check=cvg_check))
  }
}
#
# Check length of CS for subvector Delta_1 by profiling
#
subvec_Delta_1_length <- function(pars_mle,mdata_short,lower_bounds,upper_bounds,cl,xi){
  #
  mle_LL <- -loglike_short(pars_mle,mdata_short,FALSE,lower_bounds,upper_bounds)
  #
  # first check on discrete grid
  mgrid <- seq(-1.6,-0.01,length.out=400)
  LL_mgrid <- parRapply(cl,matrix(mgrid),profile_LL_Delta_1,pars_mle=pars_mle,mdata_short=mdata_short,lower_bounds=lower_bounds,upper_bounds=upper_bounds)
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
    # check that not artificially short
    if(upper_lim == max(mgrid)){
      #
      diff <- abs(mle_LL-profile_LL_Delta_1(upper_lim,pars_mle,mdata_short,lower_bounds,upper_bounds))
      while(upper_lim < upper_bounds[1] && diff < xi[j]){
        diff <- abs(mle_LL-profile_LL_Delta_1(upper_lim,pars_mle,mdata_short,lower_bounds,upper_bounds))
        upper_lim <- upper_lim + 0.001
      }
    }
    #
    if(lower_lim == min(mgrid)){
      #
      diff <- abs(mle_LL-profile_LL_Delta_1(lower_lim,pars_mle,mdata_short,lower_bounds,upper_bounds))
      while(lower_lim > lower_bounds[1] && diff < xi[j]){
        diff <- abs(mle_LL-profile_LL_Delta_1(lower_lim,pars_mle,mdata_short,lower_bounds,upper_bounds))
        lower_lim <- lower_lim - 0.001
      }
    }
    #
    CS_length[1,j] <- lower_lim
    CS_length[2,j] <- upper_lim
    CS_length[3,j] <- upper_lim - lower_lim
  }
  #
  return(CS_length)
}
#
# Check length of CS for subvector beta_1 by profiling
#
subvec_beta_1_length <- function(pars_mle,mdata_short,lower_bounds,upper_bounds,cl,xi){
  #
  mle_LL <- -loglike_short(pars_mle,mdata_short,FALSE,lower_bounds,upper_bounds)
  #
  # first check on discrete grid
  mgrid <- seq(-0.2,1,length.out=400)
  LL_mgrid <- parRapply(cl,matrix(mgrid),profile_LL_beta_1,pars_mle=pars_mle,mdata_short=mdata_short,lower_bounds=lower_bounds,upper_bounds=upper_bounds)
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
    # check that not artificially short
    if(upper_lim == max(mgrid)){
      #
      diff <- abs(mle_LL-profile_LL_beta_1(upper_lim,pars_mle,mdata_short,lower_bounds,upper_bounds))
      while(upper_lim < upper_bounds[3] && diff < xi[j]){
        diff <- abs(mle_LL-profile_LL_beta_1(upper_lim,pars_mle,mdata_short,lower_bounds,upper_bounds))
        upper_lim <- upper_lim + 0.001
      }
    }
    #
    if(lower_lim == min(mgrid)){
      #
      diff <- abs(mle_LL-profile_LL_beta_1(lower_lim,pars_mle,mdata_short,lower_bounds,upper_bounds))
      while(lower_lim > lower_bounds[3] && diff < xi[j]){
        diff <- abs(mle_LL-profile_LL_beta_1(lower_lim,pars_mle,mdata_short,lower_bounds,upper_bounds))
        lower_lim <- lower_lim - 0.001
      }
    }
    #
    CS_length[1,j] <- lower_lim
    CS_length[2,j] <- upper_lim
    CS_length[3,j] <- upper_lim - lower_lim
  }
  #
  return(CS_length)
}
#
# Calculate critical values for procedure 2 for Delta_1
#
p2_crit_Delta_1_weighted <- function(pars_mle,mcmc_draws,mdata_short,weights,alpha_vec,cl){
  #
  n_alpha <- length(alpha_vec)
  B <- nrow(mcmc_draws)
  profile_qlr <- numeric(B)
  #
  mle_LL <- -loglike_short(pars_mle,mdata_short,FALSE)
  initial_vals <- par_trans(pars_mle[-1],lower_bounds[-1],upper_bounds[-1])
  #
  profile_LL <- parRapply(cl,mcmc_draws,p2_crit_Delta_1_fn,initial_vals=initial_vals,pars_mle=pars_mle,mdata_short=mdata_short,lower_bounds=lower_bounds,upper_bounds=upper_bounds)
  if(any(is.nan(profile_LL))){profile_LL[is.nan(profile_LL)]<--1000000}
  ret <- wquantile(mle_LL-profile_LL,weights,alpha_vec)
  #
  return(ret)
}
#
# Inner part for profile QLR 
#
p2_crit_Delta_1_fn <- function(pars_temp,initial_vals,pars_mle,mdata_short,lower_bounds,upper_bounds){
  #
  ptemp <- par_to_prob(pars_temp)
  #
  lower_bd  <- find_Delta_1_lower(p0=ptemp,fixed_ix=1,initial_vals,lower_bounds,upper_bounds,tol=1e-7)
  upper_bd  <- find_Delta_1_upper(p0=ptemp,initial_vals,lower_bounds,upper_bounds,tol=1e-7)
  #
  LL_lower <- profile_LL_Delta_1(lower_bd,pars_mle,mdata_short,lower_bounds,upper_bounds)
  LL_upper <- profile_LL_Delta_1(upper_bd,pars_mle,mdata_short,lower_bounds,upper_bounds)
  #
  ret <- min(LL_lower,LL_upper)
  #
  return(ret)
}
#
# Calculate critical values for procedure 2 for beta_1
#
p2_crit_beta_1_weighted <- function(pars_mle,mcmc_draws,mdata_short,weights,alpha_vec,cl){
  #
  n_alpha <- length(alpha_vec)
  B <- nrow(mcmc_draws)
  profile_qlr <- numeric(B)
  #
  mle_LL <- -loglike_short(pars_mle,mdata_short,FALSE)
  initial_vals <- par_trans(pars_mle[-3],lower_bounds[-3],upper_bounds[-3])
  #
  profile_LL <- parRapply(cl,mcmc_draws,p2_crit_beta_1_fn,initial_vals=initial_vals,pars_mle=pars_mle,mdata_short=mdata_short,lower_bounds=lower_bounds,upper_bounds=upper_bounds)
  if(any(is.nan(profile_LL))){profile_LL[is.nan(profile_LL)]<--1000000}
  ret <- wquantile(mle_LL-profile_LL,weights,alpha_vec)
  #
  return(ret)
}
#
# Inner part for profile QLR 
#
p2_crit_beta_1_fn <- function(pars_temp,initial_vals,pars_mle,mdata_short,lower_bounds,upper_bounds){
  #
  ptemp <- par_to_prob(pars_temp)
  #
  lower_bd  <- find_beta_1_lower(p0=ptemp,fixed_ix=3,initial_vals,lower_bounds,upper_bounds,tol=1e-7)
  upper_bd  <- find_beta_1_upper(p0=ptemp,fixed_ix=3,initial_vals,lower_bounds,upper_bounds,tol=1e-7)
  #
  LL_lower <- profile_LL_beta_1(lower_bd,pars_mle,mdata_short,lower_bounds,upper_bounds)
  LL_upper <- profile_LL_beta_1(upper_bd,pars_mle,mdata_short,lower_bounds,upper_bounds)
  #
  ret <- min(LL_lower,LL_upper)
  #
  return(ret)
}
#
# END









