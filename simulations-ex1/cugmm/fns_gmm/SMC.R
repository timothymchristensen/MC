#
# SMC functions
#
SMC <- function(mdata_short,W_mat,trans_pars,keep,scalepar,cl){
  #
  # Draw from flat prior
  beta <- runif(keep)
  if(mdata_short[2]<1){
    gamma <- runif(keep)
  }else{
    gamma <- numeric(keep)+1
  }
  mu <- runif(keep)*gamma + beta*(1-gamma)
  #
  # Initialize
  draws_old <- cbind(mu,beta,gamma)
  draws_old <- log(draws_old/(1-draws_old))
  w_old <- matrix(1,keep,1)
  phi_old <- 0
  #
  lambda <- scalepar$lambda
  Nphi <- scalepar$Nphi
  Nmh <- scalepar$Nmh
  scale_cov <- scalepar$scale_cov
  #
  AR <- ESS <- numeric(Nphi)
  #
  # Recursion
  for(n in 2:Nphi){
    #
    # Correction step
    phi_new <- ((n-1)/(Nphi-1))^lambda
    phi_diff <- phi_new-phi_old
    #
    v <- exp(-phi_diff*mdata_short[3]*parRapply(cl,draws_old,gmm_objfn,mdata_short=mdata_short,W_mat=W_mat,transformed=TRUE))
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
    mh_out <- parRapply(cl,draws_new,MH_resample,mdata_short=mdata_short,W_mat=W_mat,phi=phi_new,scale_cov=scale_cov,Nmh=Nmh)
    mh_out <- t(matrix(mh_out,ncol=keep))
    #
    draws_new <- mh_out[,(1:3)]
    accept_rate <- sum(mh_out[,4]/keep)
    #
    # Update
    draws_old <- draws_new
    phi_old <- phi_new
    w_old <- w_new
    #
    AR[n] <- accept_rate
    ESS[n] <- ess
  }
  #
  draws <- exp(draws_old)/(1+exp(draws_old))
  draws[is.nan(draws)] <- 1
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
MH_resample <- function(pars,mdata_short,W_mat,phi,scale_cov,Nmh){
  #
  pars <- matrix(pars)
  lp_pars <- PKernel_temp(c(pars),mdata_short,W_mat,phi)
  #
  j <- 1
  acceptance <- 0
  #
  while(j <= Nmh){
    #
    if(mdata_short[2]<1){
      proposal <- pars + scale_cov*matrix(rnorm(3))
    }else{
      proposal <-  pars + rbind(scale_cov*matrix(rnorm(2)),0)
    }
    #
    lp_prop <- PKernel_temp(c(proposal),mdata_short,W_mat,phi)
    #
    if(is.nan(lp_prop) || is.na(lp_prop)){
      lp_prop <- -1000000
    }
    #
    if(runif(1) < exp(lp_prop-lp_pars)){
      pars <- proposal
      acceptance <- acceptance + 1
      lp_pars <- lp_prop
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
PKernel_temp <- function(pars,mdata_short,W_mat,phi){
  #
  # flat prior?
  prior <- 0
  # should also add a Jacobian term if we're using 'transformed' parameters
  if(mdata_short[2]<1){
    jacobian <- sum(pars - 2*log(1+exp(pars)))
  }else{
    jacobian <- sum(pars[1:2] - 2*log(1+exp(pars[1:2])))
  }
  #
  ret <- -phi*mdata_short[3]*gmm_objfn(pars,mdata_short,W_mat,TRUE) + prior + jacobian
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
