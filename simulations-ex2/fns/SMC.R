#
# SMC functions
#
SMC <- function(mdata_short,lower_bounds,upper_bounds,keep,scalepar,cl){
  #
  # Draw from flat prior
  Delta1 <- runif(keep)*(upper_bounds[1]-lower_bounds[1])+lower_bounds[1]
  Delta2 <- runif(keep)*(upper_bounds[2]-lower_bounds[2])+lower_bounds[2]
  beta1  <- runif(keep)*(upper_bounds[3]-lower_bounds[3])+lower_bounds[3]
  beta2  <- runif(keep)*(upper_bounds[4]-lower_bounds[4])+lower_bounds[4]
  rho    <- runif(keep)*(upper_bounds[5]-lower_bounds[5])+lower_bounds[5]
  s      <- runif(keep)*(upper_bounds[6]-lower_bounds[6])+lower_bounds[6]
  #
  # Initialize
  draws_old <- cbind(Delta1,Delta2,beta1,beta2,rho,s)
  draws_old <- t(apply(draws_old,1,function(x) par_trans(x,lower_bounds,upper_bounds)))
  w_old <- matrix(1,keep,1)
  phi_old <- 0
  #
  lambda <- scalepar$lambda
  Nphi <- scalepar$Nphi
  Nmh <- scalepar$Nmh
  Nblocks <- scalepar$Nblocks
  scale_cov <- scalepar$scale_cov
  scale_cov_mat <- diag(6)
  #
  AR <- ESS <- numeric(Nphi)
  n <- 2
  #
  # Recursion
  while(n <= Nphi){
    #
    # Correction step
    phi_new <- ((n-1)/(Nphi-1))^lambda
    phi_diff <- phi_new-phi_old
    #
    v <- exp(-phi_diff*parRapply(cl,draws_old,loglike_short_true,mdata_short=mdata_short,lower_bounds=lower_bounds,upper_bounds=upper_bounds))
    #
    w_new <- v*w_old
    w_new <- keep*w_new/sum(w_new)
    #
    # Selection step
    ess <- keep/(sum(w_new^2)/keep)
    #
    if(ess < keep/2){
      cw_new <- cumsum(w_new/sum(w_new))
      w_new <- matrix(1,keep,1)
      ix_new <- parRapply(cl,w_new,multinomial_resample,cum_weights=cw_new)
      draws_new <- draws_old[ix_new,]
    }else{
      draws_new <- draws_old
      w_new <- w_new
    }
    #
    # Mutation step: propagate via RWMH
    if(n > 2){
      scale_cov <- scale_rule(scale_cov,accept_rate)
    }
    mh_out <- parRapply(cl,draws_new,MH_resample,mdata_short=mdata_short,phi=phi_new,scale_cov=scale_cov,Nmh=Nmh,lower_bounds=lower_bounds,upper_bounds=upper_bounds)
    mh_out <- t(matrix(mh_out,ncol=keep))
    #
    draws_new <- mh_out[,(1:6)]
    accept_rate <- sum(mh_out[,7]/keep)
    #
    # Update
    draws_old <- draws_new
    phi_old <- phi_new
    w_old <- w_new
    #
    AR[n] <- accept_rate
    ESS[n] <- ess
    #
    n <- n+1
  }
  #
  draws <- t(apply(draws_old,1,function(x) par_inv_trans(x,lower_bounds,upper_bounds)))
  #
  return(list(Draws=draws,weights=w_old,accept_ratio=accept_rate,AR=AR,ESS=ESS))
}
#
# Multinomial resampler for selection step
#
multinomial_resample <- function(x,cum_weights){
  #
  u <- runif(1)
  j <- sum(u>cum_weights)+1
  return(j)
}
#
# Propose/accept/reject step for MH step
#
MH_resample <- function(pars,mdata_short,phi,scale_cov,Nmh,lower_bounds=NULL,upper_bounds=NULL){
  #
  pars <- matrix(pars)
  lp_pars <- PKernel_temp(c(pars),mdata_short,phi,lower_bounds,upper_bounds)
  #
  j <- 1
  acceptance <- 0
  #
  while(j <= Nmh){
    proposal <- pars + scale_cov*matrix(rnorm(6))
    #
    lp_prop <- PKernel_temp(c(proposal),mdata_short,phi,lower_bounds,upper_bounds)
    #
    if(is.nan(exp(lp_prop-lp_pars))==FALSE){
      if(runif(1) < exp(lp_prop-lp_pars)){
        pars <- proposal
        acceptance <- acceptance + 1
        lp_pars <- lp_prop
      }
    }else{
      pars <- pars
      acceptance <- acceptance + 0
      lp_pars <- lp_pars
    }
    #
    j <- j+1
  }
  acceptance <- acceptance/Nmh
  #
  return(c(pars,acceptance))
}
#
# Tempered kernel for MH step
#
PKernel_temp <- function(pars,mdata_short,phi,lower_bounds=NULL,upper_bounds=NULL){
  #
  jacobian <- sum(log(upper_bounds - lower_bounds) + pars - 2*log(1+exp(pars)))
  ret <- -phi*loglike_short_true(pars,mdata_short,lower_bounds,upper_bounds) + jacobian
  # returns the posterior evaluated at (pars)
  return(ret)
}
#
# Adaptive scaling rule for MH step
# 
scale_rule <- function(c,rate){
  #
  x  <- rate
  fx <- 0.95 + 0.10*(exp(16*(x-0.35)))/(1+exp(16*(x-0.35)))
  #
  return(c*fx)
}
#
# END

