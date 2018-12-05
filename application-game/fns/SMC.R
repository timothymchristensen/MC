#
# SMC functions
#
SMC <- function(mdata_short,lower_bounds,upper_bounds,ssame,keep,scalepar,cl){
  #
  # Draw from flat prior
  Delta1 <- runif(keep)*(upper_bounds[1]-lower_bounds[1])+lower_bounds[1]
  Delta2 <- runif(keep)*(upper_bounds[2]-lower_bounds[2])+lower_bounds[2]
  beta11 <- runif(keep)*(upper_bounds[3]-lower_bounds[3])+lower_bounds[3]
  beta12 <- runif(keep)*(upper_bounds[4]-lower_bounds[4])+lower_bounds[4]
  beta13 <- runif(keep)*(upper_bounds[5]-lower_bounds[5])+lower_bounds[5]
  beta1  <- cbind(beta11,beta12,beta13)
  beta21 <- runif(keep)*(upper_bounds[6]-lower_bounds[6])+lower_bounds[6]
  beta22 <- runif(keep)*(upper_bounds[7]-lower_bounds[7])+lower_bounds[7]
  beta23 <- runif(keep)*(upper_bounds[8]-lower_bounds[8])+lower_bounds[8]
  beta2  <- cbind(beta21,beta22,beta23)
  rho    <- runif(keep)*(upper_bounds[9]-lower_bounds[9])+lower_bounds[9]
  if(ssame){
    s <- runif(keep)*(upper_bounds[10]-lower_bounds[10])+lower_bounds[10]
    draws_old <- cbind(Delta1,Delta2,beta1,beta2,rho,s)
  }else{
    s000 <- runif(keep)*(upper_bounds[10]-lower_bounds[10])+lower_bounds[10]
    s001 <- runif(keep)*(upper_bounds[11]-lower_bounds[11])+lower_bounds[11]
    s010 <- runif(keep)*(upper_bounds[12]-lower_bounds[12])+lower_bounds[12]
    s100 <- runif(keep)*(upper_bounds[13]-lower_bounds[13])+lower_bounds[13]
    s011 <- runif(keep)*(upper_bounds[14]-lower_bounds[14])+lower_bounds[14]
    s101 <- runif(keep)*(upper_bounds[15]-lower_bounds[15])+lower_bounds[15]
    s110 <- runif(keep)*(upper_bounds[16]-lower_bounds[16])+lower_bounds[16]
    s111 <- runif(keep)*(upper_bounds[17]-lower_bounds[17])+lower_bounds[17]
    s    <- cbind(s000,s001,s010,s100,s011,s101,s110,s111)
    draws_old <- cbind(Delta1,Delta2,beta1,beta2,rho,s)
  }
  #
  # Initialize
  draws_old <- t(apply(draws_old,1,function(x) par_trans(x,lower_bounds,upper_bounds)))
  w_old <- matrix(1,keep,1)
  phi_old <- 0
  #
  lambda <- scalepar$lambda
  Nphi <- scalepar$Nphi
  Nmh <- scalepar$Nmh
  Nblocks <- scalepar$Nblocks
  scale_cov <- scalepar$scale_cov
  scale_cov_mat <- diag(ncol(draws_old))
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
    v <- exp(-phi_diff*parRapply(cl,draws_old,loglike_short_true,mdata_short=mdata_short,lower_bounds=lower_bounds,upper_bounds=upper_bounds,ssame=ssame))
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
      scale_cov_mat <- cov.wt(x=draws_new,wt=as.numeric(w_new))$cov
    }
    mh_out <- parRapply(cl,draws_new,MH_resample_blocked,mdata_short=mdata_short,phi=phi_new,scale_cov_mat=scale_cov_mat,scale_cov=scale_cov,Nmh=Nmh,Nblocks=Nblocks,lower_bounds=lower_bounds,upper_bounds=upper_bounds,ssame)
    mh_out <- t(matrix(mh_out,ncol=keep))
    #
    draws_new <- mh_out[,1:ncol(draws_old)]
    accept_rate <- sum(mh_out[,ncol(draws_old)+1]/keep)
    #
    # Update
    draws_old <- draws_new
    phi_old <- phi_new
    w_old <- w_new
    #
    cat(n,ess,accept_rate,"\n")
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
MH_resample <- function(pars,mdata_short,phi,scale_cov,Nmh,lower_bounds=NULL,upper_bounds=NULL,ssame){
  #
  pars <- matrix(pars)
  lp_pars <- PKernel_temp(c(pars),mdata_short,phi,lower_bounds,upper_bounds,ssame)
  #
  j <- 1
  acceptance <- 0
  #
  while(j <= Nmh){
    proposal <- pars + scale_cov*matrix(rnorm(pars))
    #
    lp_prop <- PKernel_temp(c(proposal),mdata_short,phi,lower_bounds,upper_bounds,ssame)
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
# Propose/accept/reject step for block MH step
#
MH_resample_blocked <- function(pars,mdata_short,phi,scale_cov_mat,scale_cov,Nmh,Nblocks,lower_bounds=NULL,upper_bounds=NULL,ssame){
  #
  pars <- matrix(pars)
  lp_pars <- PKernel_temp(c(pars),mdata_short,phi,lower_bounds,upper_bounds,ssame)
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
      if(any(block_ix==i)){
        #
        ix <- (block_ix==i)
        #
        proposal <- pars
        proposal[ix] <- pars[ix] + scale_cov*t(chol(scale_cov_mat[ix,ix]))%*%matrix(rnorm(sum(ix)))
        #
        lp_prop <- PKernel_temp(c(proposal),mdata_short,phi,lower_bounds,upper_bounds,ssame)
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
# Tempered kernel for MH step
#
PKernel_temp <- function(pars,mdata_short,phi,lower_bounds=NULL,upper_bounds=NULL,ssame){
  #
  jacobian <- sum(log(upper_bounds - lower_bounds) + pars - 2*log(1+exp(pars)))
  ret <- -phi*loglike_short_true(pars,mdata_short,lower_bounds,upper_bounds,ssame) + jacobian
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

