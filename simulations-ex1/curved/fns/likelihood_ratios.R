#
# Calculate critical values for procedure 1
#
p1_crit <- function(pars_mle,mcmc_draws,mdata,alpha_vec){
  #
  n_alpha <- length(alpha_vec)
  B <- nrow(mcmc_draws)
  mcmc_LL <- numeric(B)
  #
  mle_LL <- -loglike(pars_mle,mdata,FALSE)
  #
  for(j in 1:B){
    mcmc_LL[j] <- -loglike(mcmc_draws[j,],mdata,FALSE)
  }
  #
  LL_diff <- abs(mle_LL - mcmc_LL)
  #
  ret <- as.numeric(quantile(LL_diff,alpha_vec))
  #
  return(ret)
}
#
# Calculate critical values for procedure 2
#
p2_crit <- function(pars_mle,mcmc_draws,mdata,alpha_vec){
  #
  n_alpha <- length(alpha_vec)
  B <- nrow(mcmc_draws)
  profile_qlr <- numeric(B)
  #
  mle_LL <- -loglike(pars_mle,mdata,FALSE)
  #
  for(j in 1:B){
    #
    pars_temp <- mcmc_draws[j,]
    #
    lower_bd  <- pars_temp[1]-pars_temp[2]*(1-pars_temp[3])
    upper_bd  <- pars_temp[1]+(1-pars_temp[2])*(1-pars_temp[3])
    #
    mle_lower <- mle_mu_fixed(lower_bd,mdata)
    mle_upper <- mle_mu_fixed(upper_bd,mdata)
    #
    qlr_lower <- mle_LL+loglike(c(lower_bd,mle_lower),mdata,FALSE)
    qlr_upper <- mle_LL+loglike(c(upper_bd,mle_upper),mdata,FALSE)
    profile_qlr[j] <- max(abs(c(qlr_lower,qlr_upper)))
  }
  #
  ret <- as.numeric(quantile(profile_qlr,alpha_vec))
  #
  return(ret)
}
#
# Check coverage of for the full vector identified set
#
full_check <- function(pars_mle,pars_0,mdata,xi){
  #
  mle_LL <- -loglike(pars_mle,mdata,FALSE)
  mle_0  <- -loglike(pars_0,mdata,FALSE)
  #
  LL_diff <- abs(mle_LL - mle_0)
  #
  ret <-  (LL_diff <= xi)*1
  #
  return(ret)
}
#
# Check coverage of for subvector by profiling
#
subvec_mu_cvg <- function(pars_mle,M_I,mdata,xi){
  #
  mle_LL <- -loglike(pars_mle,mdata,FALSE)
  #
  mle_M_I <- numeric(length(M_I))
  for(i in 1:length(M_I)){
    pars_hat_mu <- mle_mu_fixed(M_I[i],mdata)
    mle_M_I[i] <- -loglike(c(M_I[i],pars_hat_mu),mdata,FALSE)
  }
  #
  LL_diff <- max(abs(mle_LL - mle_M_I))
  #
  ret <-  (LL_diff <= xi)*1
  #
  return(ret)
}
#
# Check length of CS for subvector mu by profiling
#
subvec_mu_length <- function(pars_mle,mdata,xi){
  #
  mle_LL <- -loglike(pars_mle,mdata,FALSE)
  #
  # first check on discrete grid
  mgrid <- seq(0.2,0.8,length.out=200)
  LL_mgrid <- numeric(length(mgrid))
  for(i in 1:length(mgrid)){
    #
    pars_hat_mu <- mle_mu_fixed(mgrid[i],mdata)
    LL_mgrid[i] <- -loglike(c(mgrid[i],pars_hat_mu),mdata,FALSE)
  }
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
      pars_hat_mu <- mle_mu_fixed(max(mgrid),mdata)
      diff <- abs(mle_LL+loglike(c(max(mgrid),pars_hat_mu),mdata,FALSE))
      while(upper_lim < 1 && diff < xi[j]){
        pars_hat_mu <- mle_mu_fixed(upper_lim,mdata)
        diff <- abs(mle_LL+loglike(c(upper_lim,pars_hat_mu),mdata,FALSE))
        upper_lim <- upper_lim + 0.002
      }
    }
    #
    if(lower_lim == min(mgrid)){
      #
      pars_hat_mu <- mle_mu_fixed(min(mgrid),mdata)
      diff <- abs(mle_LL+loglike(c(min(mgrid),pars_hat_mu),mdata,FALSE))
      while(lower_lim > 0 && diff < xi[j]){
        pars_hat_mu <- mle_mu_fixed(lower_lim,mdata)
        diff <- abs(mle_LL+loglike(c(lower_lim,pars_hat_mu),mdata,FALSE))
        lower_lim <- lower_lim - 0.002
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
# Percentile CS
#
percentile_CS <- function(chain,alpha_vec,M_I=NULL){
  #
  lower <- as.numeric(quantile(chain,(1-alpha_vec)/2))
  upper <- as.numeric(quantile(chain,1-(1-alpha_vec)/2))
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
#END