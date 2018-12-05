#
# SMC functions
#
SMC_curved <- function(mdata_short,trans_pars,keep,scalepar,loglike_short,PKernel_temp,cl){
  #
  # Draw from curved prior
  beta <- rbeta(keep,3,8)
  gamma <- rbeta(keep,8,1)
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
    v <- exp(-phi_diff*parRapply(cl,draws_old,loglike_short,mdata_short=mdata_short))
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
    mh_out <- parRapply(cl,draws_new,MH_resample,mdata_short=mdata_short,phi=phi_new,scale_cov=scale_cov,Nmh=Nmh,TRUE,PKernel_temp=PKernel_temp,loglike_short=loglike_short)
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
# Multinomial resampler for selection step
#
multinomial_resample_2 <- function(x,cum_weights){
  #
  u <- runif(1)
  j <- 1
  while(cum_weights[j]<=u){
    j <- j+1
  }
  return(j)
}
#
# Propose/accept/reject step for MH step
#
MH_resample <- function(pars,mdata_short,phi,scale_cov,Nmh,transformed=TRUE,PKernel_temp=PKernel_temp,loglike_short=loglike_short){
  #
  pars <- matrix(pars)
  lp_pars <- PKernel_temp(c(pars),mdata_short,phi,transformed,loglike_short)
  #
  j <- 1
  acceptance <- 0
  #
  while(j <= Nmh){
    proposal <- pars + scale_cov*matrix(rnorm(3))
    #
    lp_prop <- PKernel_temp(c(proposal),mdata_short,phi,transformed,loglike_short)
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
PKernel_temp <- function(pars,mdata_short,phi,transformed=TRUE,loglike_short=loglike_short){
  #
  mu <- beta <- gamma <- 0
  #
  if(transformed==TRUE){
    mu <- exp(pars[1])/(1+exp(pars[1]))    # Pr(Y_i = 1)
    beta  <- exp(pars[2])/(1+exp(pars[2])) # Pr(Y_i = 1 | D_i = 0)
    gamma <- exp(pars[3])/(1+exp(pars[3])) # Pr(D_i = 1)
  }else{
    mu <- pars[1]     # Pr(Y_i = 1)
    beta  <- pars[2]  # Pr(Y_i = 1 | D_i = 0)
    gamma <- pars[3]  # Pr(D_i = 1)
  }
  # flat prior?
  prior <- log(dbeta(beta,3,8)) + log(dbeta(gamma,8,1))
  # should also add a Jacobian term if we're using 'transformed' parameters
  jacobian <- 0
  if(transformed==TRUE){
    jacobian <- sum(pars - 2*log(1+exp(pars)))
  }
  #
  ret <- -phi*loglike_short(pars,mdata_short,transformed) + prior + jacobian
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
# short log-likelihood for Example 1
#
loglike_short <- function(pars,mdata_short,transformed=TRUE){
  #
  mu <- beta <- gamma <- 0
  #
  if(transformed==TRUE){
    mu    <- exp(pars[1])/(1+exp(pars[1]))    # Pr(Y_i = 1)
    beta  <- exp(pars[2])/(1+exp(pars[2]))    # Pr(Y_i = 1 | D_i = 0)
    gamma <- exp(pars[3])/(1+exp(pars[3]))    # Pr(D_i = 1)
  }else{
    mu    <- pars[1]    # Pr(Y_i = 1)
    beta  <- pars[2]    # Pr(Y_i = 1 | D_i = 0)
    gamma <- pars[3]    # Pr(D_i = 1)
  }
  # for convenience, define
  rho <- mu - beta*(1-gamma)
  #
  D  <- mdata_short[1]
  YD <- mdata_short[2]
  n  <- mdata_short[3]
  #
  Term1 <- YD*log(rho)
  #
  Term2 <- 0
  if(sum(D-YD)!=0){
    Term2 <- (D - YD)*log(gamma - rho)
  }
  #
  Term3 <- 0
  if(sum(D)!=n){
    Term3 <- (n - D)*log(1 - gamma)
  }
  # sum the loglikelihood terms
  ret <- sum(Term1 + Term2 + Term3)
  #
  # return MINUS the log-likelihood (as we're running a minimization routine)
  return(-ret)
}

  
