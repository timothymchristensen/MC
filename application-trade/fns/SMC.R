#
# SMC functions
#
SMC <- function(data,het,prior_par,keep,scalepar,mle_out,cl){
  #
  # Initial draw
  pars_mle <- mle_out$mle
  s_mle <- chol(mle_out$cov)
  i_mle <- solve(s_mle)
  draws_old <- replicate(44+2*het,rnorm(keep))%*%s_mle+t(replicate(keep,pars_mle))
  #
  # Initialize
  w_old <- matrix(1,keep,1)
  phi_old <- 0
  accept_rate <- 0
  #
  lambda <- scalepar$lambda
  Nphi <- scalepar$Nphi
  Nmh <- scalepar$Nmh
  Nblocks <- scalepar$Nblocks
  scale_cov <- scalepar$scale_cov
  scale_cov_mat <- diag(ncol(draws_old))
  #
  AR <- ESS <- numeric(Nphi)
  #
  # Recursion
  n <- 2
  while(n<=Nphi){
    #
    # Correction step
    phi_new <- ((n-1)/(Nphi-1))^lambda
    phi_diff <- phi_new-phi_old
    #
    ll <- parRapply(cl,draws_old,loglike_trade,data,het)
    lp <- apply(draws_old,1,log_prior,prior_par)
    li <- log_initial(draws_old,pars_mle,i_mle)
    v <- exp(phi_diff*(lp-ll-li))
    w_new <- v*w_old
    w_new <- keep*w_new/sum(w_new)
    w_new[is.na(w_new)] <- 0
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
    scale_cov_mat <- cov.wt(x=draws_new,wt=as.numeric(w_new))$cov
    if(n > 2){
      scale_cov <- scale_rule(scale_cov,accept_rate)
    }
    mh_out <- parRapply(cl,draws_new,MH_resample_blocked,data=data,het=het,prior_par=prior_par,phi=phi_new,pars_mle=pars_mle,i_mle=i_mle,scale_cov=scale_cov,scale_cov_mat=scale_cov_mat,Nmh=Nmh,Nblocks=Nblocks)
    mh_out <- t(matrix(mh_out,ncol=keep))
    #
    draws_new <- mh_out[,(1:(ncol(mh_out)-1))]
    accept_rate <- sum(mh_out[,ncol(mh_out)]/keep)
    #
    # Update
    draws_old <- draws_new
    phi_old <- phi_new
    w_old <- w_new
    #
    AR[n] <- accept_rate
    ESS[n] <- ess
    #
    cat("Iteration:",n,"ESS",ess,"Avg log like",mean(ll),"AR",accept_rate,"\n")
    #
    n <- n+1
  }
  #
  draws <- draws_new
  #
  return(list(Draws=draws,weights=w_old,accept_ratio=accept_rate,AR=AR,ESS=ESS))
}
#
# Multinomial resampler for selection step
#
multinomial_resample <- function(x,cum_weights){
  #
  u <- runif(1)
  j <- 1
  while(cum_weights[j]<=u){
    j <- j+1
  }
  return(j)
}
#
# Tempered kernel for MH step
#
PKernel_trade <- function(pars,data,het,prior_par,phi,pars_mle,i_mle){
  #
  # prior contribution
  prior <- log_prior(pars,prior_par)
  #
  # initializing contribution
  initial <- log_initial_pk(pars,pars_mle,i_mle)
  #
  ret <- phi*(prior - loglike_trade(pars,data,het)) + (1-phi)*initial
  # returns the posterior evaluated at (pars)
  return(ret)
}
#
# Adaptive scaling rule for MH step
# 
scale_rule <- function(c,rate){
  #
  x  <- rate
  fx <- 0.95 + 0.1*(exp(16*(x-0.35)))/(1+exp(16*(x-0.35)))
  #
  return(c*fx)
}
#
# Propose/accept/reject step for block MH step
#
MH_resample_blocked <- function(pars,data,het,prior_par,phi,pars_mle,i_mle,scale_cov,scale_cov_mat,Nmh,Nblocks){
  #
  pars <- matrix(pars)
  npars <- nrow(pars)
  lp_pars <- PKernel_trade(c(pars),data,het,prior_par,phi,pars_mle,i_mle)
  #
  j <- 1
  acceptance <- 0
  #
  # block indices
  block_ix <- floor(runif(length(pars))*Nblocks)+1
  #
  while(j <= Nmh){
    #
    for(i in 1:Nblocks){
      #
      if(sum(block_ix==i)>0){
        #
        ix <- (block_ix==i)
        #
        proposal <- pars
        proposal[ix] <- pars[ix] + scale_cov*t(chol(scale_cov_mat[ix,ix]))%*%matrix(rnorm(sum(ix)))
        #
        lp_prop <- PKernel_trade(c(proposal),data,het,prior_par,phi,pars_mle,i_mle)
        #
        if(is.nan(lp_prop) || is.na(lp_prop) || is.infinite(lp_prop)){
          lp_prop <- -1000000
        }
        #
        if(runif(1) < exp(lp_prop-lp_pars)){
          pars[ix] <- proposal[ix]
          acceptance <- acceptance + 1
          lp_pars <- lp_prop
        }else{
          pars <- pars
          acceptance <- acceptance + 0
          lp_pars <- lp_pars
        }
      }
    }
    #
    j <- j+1
  }
  acceptance <- acceptance/(Nmh*Nblocks)
  #
  return(c(pars,acceptance))
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
# log prior
#
log_prior <- function(pars,prior_par){
  #
  prior <- sum(log(pnorm(pars/prior_par)/prior_par))
  #
  return(prior)
}
#
# log of initializing distribution 
#
log_initial <- function(pars_mat,pars_mle,i_mle){
  #
  prior <- rowSums(log(pnorm((pars_mat-matrix(1,nrow(pars_mat),1)%*%t(pars_mle))%*%i_mle)))+log(det(i_mle))
  #
  return(prior)
}
#
# log of initializing distribution for posterior kernel
#
log_initial_pk <- function(pars,pars_mle,i_mle){
  #
  prior <- sum(log(pnorm(t(pars-pars_mle)%*%i_mle)))+log(det(i_mle))
  #
  return(prior)
}
#
# END
