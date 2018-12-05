#
# Calculate critical values for procedure 1
#
p1_crit_weighted <- function(pars_mle,mcmc_draws,mdata,weights,alpha_vec){
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
  ret <- wquantile(LL_diff,weights,alpha_vec)
  #
  return(ret)
}
#
# Calculate critical values for procedure 2
#
p2_crit_weighted <- function(pars_mle,mcmc_draws,mdata,weights,alpha_vec){
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
  ret <- wquantile(profile_qlr,weights,alpha_vec)
  #
  return(ret)
}
#
# Calculate critical values for procedure 2 via parrelization
#
p2_crit_par_weighted <- function(pars_mle,mcmc_draws,mdata,weights,alpha_vec,cl){
  #
  n_alpha <- length(alpha_vec)
  B <- nrow(mcmc_draws)
  profile_qlr <- numeric(B)
  #
  mle_LL <- -loglike(pars_mle,mdata,FALSE)
  #
  profile_qlr <- parRapply(cl,mcmc_draws,p2_crit_par_fn,mdata=mdata,mle_LL=mle_LL,mle_mu_fixed_par=mle_mu_fixed_par,loglike=loglike,loglike_constr=loglike_constr)
  #
  ret <- wquantile(profile_qlr,weights,alpha_vec)
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
#END