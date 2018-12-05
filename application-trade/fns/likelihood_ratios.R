#
# Calculate critical values for procedure 1
#
p1_crit_weighted <- function(pars_mle,mcmc_draws,data,het,weights,alpha_vec,cl){
  #
  n_alpha <- length(alpha_vec)
  B <- nrow(mcmc_draws)
  mcmc_LL <- numeric(B)
  #
  mle_LL <- -loglike_trade(pars_mle,data,het)
  #
  mcmc_LL <- -apply(mcmc_draws,1,loglike_trade,data=data,het=het)
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
p2_crit_weighted <- function(fixed_ix,mle_out,mcmc_draws,data,het,weights,alpha_vec,cl){
  #
  n_alpha <- length(alpha_vec)
  B <- nrow(mcmc_draws)
  profile_qlr <- numeric(B)
  #
  mle_LL <- -loglike_trade(mle_out$mle,data,het)
  #
  profile_qlr <- mle_LL-parRapply(cl,mcmc_draws,p2_crit_fn,fixed_ix=fixed_ix,mle_out=mle_out,data=data,het=het)
  if(any(is.nan(profile_qlr))){profile_qlr[is.nan(profile_qlr)]<--1000000}
  ret <- wquantile(profile_qlr,weights,alpha_vec)
  #
  return(ret)
}
#
# Inner part for profile QLR 
#
p2_crit_fn <- function(pars_temp,fixed_ix,mle_out,data,het){
  #
  initial_vals <- mle_out$mle
  step <- max(0.1*mle_out$ses[fixed_ix],0.01)
  #
  lower_bd  <- find_int_lower(fixed_ix,pars_temp,step,data,het)
  upper_bd  <- find_int_upper(fixed_ix,pars_temp,step,data,het)
  #
  if(upper_bd != lower_bd){
    LL_lower <- profile_LL(lower_bd,fixed_ix,pars_mle,data,het)
    LL_upper <- profile_LL(upper_bd,fixed_ix,pars_mle,data,het)
    ret <- min(LL_lower,LL_upper)
  }else{
    ret <- profile_LL(lower_bd,fixed_ix,initial_vals,data,het)
  }
  #
  return(ret)
}
#
# Check length of CS for subvectors by profiling
#
subvec_length <- function(fixed_ix,mle_out,data,het,cl,xi){
  #
  mle_LL <- -loglike_trade(mle_out$mle,data,het)
  #
  # check on discrete grid
  scale <- 8
  mgrid <- seq(mle_out$mle[fixed_ix]-scale*mle_out$ses[fixed_ix],mle_out$mle[fixed_ix]+scale*mle_out$ses[fixed_ix],length.out=400)
  #
  LL_mgrid <- numeric(length(mgrid))
  for(i in 1:length(mgrid)){
    LL_mgrid[i] <- profile_LL(mgrid[i],fixed_ix,initial_vals=mle_out$mle,data=data,het=het)
    cat("iteration: ",i," index: ",fixed_ix,"\n")
  }
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
  return(list(CS_length=CS_length,mgrid=mgrid))
}
#
# END









