#
# Calculate critical values for procedure 1
#
p1_crit_weighted <- function(pars_hat,mcmc_draws,mdata_short,W_mat,weights,alpha_vec){
  #
  n_alpha <- length(alpha_vec)
  B <- nrow(mcmc_draws)
  n <- mdata_short[3]
  mcmc_crit <- numeric(B)
  #
  gmm_crit <- gmm_objfn(pars_hat,mdata_short,W_mat,FALSE)
  #
  for(j in 1:B){
    mcmc_crit[j] <- gmm_objfn(mcmc_draws[j,],mdata_short,W_mat,FALSE)
  }
  #
  LL_diff <- abs(2*n*(mcmc_crit - gmm_crit))
  #
  ret <- wquantile(LL_diff,weights,alpha_vec)
  #
  return(ret)
}
#
# Check coverage of for the full vector identified set
#
full_check <- function(pars_hat,pars_0,mdata_short,W_mat,xi){
  #
  gmm_hat <-  gmm_objfn(pars_hat,mdata_short,W_mat,FALSE)
  gmm_0  <-  gmm_objfn(pars_0,mdata_short,W_mat,FALSE)
  n <- mdata_short[3]
  #
  LL_diff <- as.numeric(2*n*abs(gmm_0 - gmm_hat))
  #
  ret <-  (LL_diff <= xi)*1
  #
  return(ret)
}
#
# Check coverage of for subvector by profiling
#
subvec_mu_cvg <- function(pars_hat,M_I,mdata_short,W_mat,xi){
  #
  gmm_hat <-  gmm_objfn(pars_hat,mdata_short,W_mat,FALSE)
  #
  n <- mdata_short[3]
  #
  gmm_M_I <- numeric(length(M_I))
  if(mdata_short[2]<1){
    for(i in 1:length(M_I)){
      pars_hat_mu <- gmm_mu_subvec(M_I[i],mdata_short,W_mat)
      gmm_M_I[i] <- gmm_objfn(c(M_I[i],pars_hat_mu),mdata_short,W_mat,FALSE)
    }
  }else{
    for(i in 1:length(M_I)){
      pars_hat_mu <- c(0.5,1)
      gmm_M_I[i] <- gmm_objfn(c(M_I[i],pars_hat_mu),mdata_short,W_mat,FALSE)
    }
  }
  #
  LL_diff <- max(2*n*(gmm_M_I - gmm_hat))
  #
  ret <-  (LL_diff <= xi)*1
  #
  return(ret)
}
#
# Check length of CS for subvector mu by profiling
#
subvec_mu_length <- function(pars_hat,mdata_short,W_mat,xi){
  #
  gmm_hat <- gmm_objfn(pars_hat,mdata_short,W_mat,FALSE)
  n <- mdata_short[3]
  pars_hat_mu <- c(0.5,1)
  #
  # first check on discrete grid
  mgrid <- seq(0.2,0.8,length.out=200)
  gmm_mgrid <- numeric(length(mgrid))
  for(i in 1:length(mgrid)){
    #
    if(mdata_short[2]<1){pars_hat_mu <- gmm_mu_subvec(mgrid[i],mdata_short,W_mat)}
    gmm_mgrid[i] <- gmm_objfn(c(mgrid[i],pars_hat_mu),mdata_short,W_mat,FALSE)
  }
  gmm_diff <- 2*n*(abs(gmm_mgrid - gmm_hat))
  #
  CS_length <- matrix(NA,3,length(xi))
  #
  for(j in 1:length(xi)){
    ix <- (gmm_diff <= xi[j])
    upper_lim <- max(mgrid[ix])
    lower_lim <- min(mgrid[ix])
    #
    # check that not artificially short
    if(upper_lim == max(mgrid)){
      #
      if(mdata_short[2]<1){pars_hat_mu <- gmm_mu_subvec(max(mgrid),mdata_short,W_mat)}
      diff <- 2*n*(gmm_objfn(c(max(mgrid),pars_hat_mu),mdata_short,W_mat,FALSE)-gmm_hat)
      while(upper_lim < 1 && diff < xi[j]){
        if(mdata_short[2]<1){pars_hat_mu <- gmm_mu_subvec(upper_lim,mdata_short,W_mat)}
        diff <- 2*n*(gmm_objfn(c(upper_lim,pars_hat_mu),mdata_short,W_mat,FALSE)-gmm_hat)
        upper_lim <- upper_lim + 0.002
      }
    }
    #
    if(lower_lim == min(mgrid)){
      #
      if(mdata_short[2]<1){pars_hat_mu <- gmm_mu_subvec(min(mgrid),mdata_short,W_mat)}
      diff <- 2*n*(gmm_objfn(c(max(mgrid),pars_hat_mu),mdata_short,W_mat,FALSE)-gmm_hat)
      while(lower_lim > 0 && diff < xi[j]){
        if(mdata_short[2]<1){pars_hat_mu <- gmm_mu_subvec(lower_lim,mdata_short,W_mat)}
        diff <- 2*n*(gmm_objfn(c(lower_lim,pars_hat_mu),mdata_short,W_mat,FALSE)-gmm_hat)
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
# weighted quantile
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
# Calculate critical values for procedure 2 via parrelization
#
p2_crit_par_weighted <- function(pars_hat,mcmc_draws,mdata_short,W_mat,weights,alpha_vec,cl){
  #
  n_alpha <- length(alpha_vec)
  B <- nrow(mcmc_draws)
  profile_qlr <- numeric(B)
  #
  gmm_hat <- gmm_objfn(pars_hat,mdata_short,W_mat,FALSE)
  #
  if(mdata_short[2]<1){
    profile_qlr <- parRapply(cl,mcmc_draws,p2_crit_par_fn,mdata_short=mdata_short,W_mat=W_mat,gmm_hat=gmm_hat)
  }else{
    profile_crit <- parRapply(cl,mcmc_draws,gmm_objfn,mdata_short=mdata_short,W_mat=W_mat,transformed=FALSE)
    profile_qlr <- 2*mdata_short[3]*(profile_crit-gmm_hat)
  }
  #
  ret <- wquantile(profile_qlr,weights,alpha_vec)
  # 
  return(ret)
}
#
# Inner part for profile QLR 
#
p2_crit_par_fn <- function(pars_temp,mdata_short,W_mat,gmm_hat){
  #
  lower_bd  <- pars_temp[1]-pars_temp[2]*(1-pars_temp[3])
  upper_bd  <- pars_temp[1]+(1-pars_temp[2])*(1-pars_temp[3])
  n <- mdata_short[3]
  #
  par_hat_lower <- gmm_mu_subvec(lower_bd,mdata_short,W_mat)
  par_hat_upper <- gmm_mu_subvec(upper_bd,mdata_short,W_mat)
  #
  qlr_lower <- 2*n*(gmm_objfn(c(lower_bd,par_hat_lower),mdata_short,W_mat,FALSE)-gmm_hat)
  qlr_upper <- 2*n*(gmm_objfn(c(upper_bd,par_hat_upper),mdata_short,W_mat,FALSE)-gmm_hat)
  #
  ret <- max(abs(c(qlr_lower,qlr_upper)))
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
#END