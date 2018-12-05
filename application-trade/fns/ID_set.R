#
# KL distance objective function to compute identified set
#
klobj <- function(pars0,pars1,Xmat,het,u,wu){
  #
  beta_m0 <- beta_m1 <- matrix(0,nrow=21,ncol=1)
  beta_m0[-11,] <- pars0[1:20]
  beta_m1[-11,] <- pars1[1:20]
  #
  beta_z0 <- matrix(pars0[21:41])
  beta_z1 <- matrix(pars1[21:41])
  #
  delta0 <- pars0[42]
  delta1 <- pars1[42]
  #
  sigma_m0 <- exp(pars0[43])
  sigma_m1 <- exp(pars1[43])
  #
  rho0 <- (exp(pars0[44])-1)/(1+exp(pars0[44]))
  rho1 <- (exp(pars1[44])-1)/(1+exp(pars1[44]))
  #
  if(het){
    varpi0 <- pars0[45:46]
    varpi1 <- pars1[45:46]
    #
    sigma_z0 <- exp(varpi0[1]*Xmat[,2] + varpi0[2]*Xmat[,2]^2)
    sigma_z1 <- exp(varpi1[1]*Xmat[,2] + varpi1[2]*Xmat[,2]^2)
  }else{
    sigma_z0 <- sigma_z1 <- sigma_z00 <- sigma_z01 <- 1
  }
  sigma_v0 <- sqrt(sigma_m0^2+2*rho0*delta0*sigma_m0*sigma_z0+(delta0*sigma_z0)^2)
  sigma_v1 <- sqrt(sigma_m1^2+2*rho1*delta1*sigma_m1*sigma_z1+(delta1*sigma_z1)^2)
  #
  r_Sigma0 <- (rho0*sigma_m0*sigma_z0+(sigma_z0^2)*delta0)/(sigma_v0*sigma_z0)
  r_Sigma1 <- (rho1*sigma_m1*sigma_z1+(sigma_z1^2)*delta1)/(sigma_v1*sigma_z1)
  #
  if(any(is.na(r_Sigma0))==TRUE || any(is.nan(r_Sigma0))==TRUE || any(is.na(r_Sigma1))==TRUE || any(is.nan(r_Sigma1))==TRUE){
    return(1E10)
  }
  if(any(r_Sigma0 <= -1)==TRUE || any(r_Sigma0 >= 1)==TRUE || any(r_Sigma1 <= -1)==TRUE || any(r_Sigma1 >= 1)==TRUE){
    return(1E10)
  }
  #
  p00 <- .Internal(pnorm(Xmat%*%beta_z0 / sigma_z0,0,1,FALSE,FALSE))
  p01 <- .Internal(pnorm(Xmat%*%beta_z1 / sigma_z1,0,1,FALSE,FALSE))
  term_1 <- p00*log(p00/p01)
  #
  p10 <- .Internal(pnorm(as.numeric(Xmat%*%beta_z0) / sigma_z0,0,1,TRUE,FALSE))
  term_2 <- p10*log(sigma_v1/sigma_v0)
  #
  term_3_0 <- (Xmat%*%beta_z0)/sigma_z0
  term_3_1 <- (Xmat%*%beta_z1)/sigma_z1
  #
  # use change of variables: u = (m_{ij} - X_{ij}(\beta_m + \delta \beta_z))/\sigma_v(X_{ij}) under true parameters
  # transformed variable under alternate parameters
  u1 <- apply(u,2,function(x) (x*sigma_v0 + Xmat%*%(beta_m0+delta0*beta_z0-beta_m1-delta1*beta_z1))/sigma_v1)
  #
  # cdf terms under change of variables
  term_3_cdf0 <- apply(u ,2,function(x) .Internal(pnorm((term_3_0 + r_Sigma0*x)/sqrt(1-r_Sigma0^2),0,1,TRUE,FALSE)))
  term_3_cdf1 <- apply(u1,2,function(x) .Internal(pnorm((term_3_1 + r_Sigma1*x)/sqrt(1-r_Sigma1^2),0,1,TRUE,FALSE)))
  #
  # contribution from cdf terms
  term_3_a <- (term_3_cdf0*log(term_3_cdf0/term_3_cdf1))%*%wu
  #
  # pdf terms under change of variables
  term_3_pdf0 <- matrix(1,nrow(Xmat),1)%*%.Internal(dnorm(u,0,1,FALSE))
  term_3_pdf1 <- .Internal(dnorm(u1,0,1,FALSE))
  #
  # contribution from pdf terms
  term_3_b <- (term_3_cdf0*(log(term_3_pdf0/term_3_pdf1)))%*%wu
  #
  term_3 <- term_3_a + term_3_b
  #
  ret <- sum(term_1) + sum(term_2) + sum(term_3)
  #
  if(is.na(ret)==TRUE||is.nan(ret)==TRUE||is.infinite(ret)==TRUE){
    return(1E10)
  }else{
    return(ret)
  }
}
#
# KL distance objective function with some fixed parameters
#
klobj_fixed <- function(free_pars,fixed_pars,fixed_ix,pars0,Xmat,het,u,wu){
  #
  pars <- numeric(length(pars0))
  #
  pars[-fixed_ix] <- free_pars
  pars[fixed_ix] <- fixed_pars
  #
  ret <- klobj(pars0,pars,Xmat,het,u,wu)
  #
  return(ret)
}
#
# Minimized KL distance objective function with a subset of parameters fixed 
#
kldist_fixed <- function(fixed_pars,fixed_ix,pars0,initial_vals,Xmat,het){
  #
  # quadrature nodes and weights
  u <- t(c(-4.85946282833231215015516494660,
           -3.58182348355192692277623675546,
           -2.48432584163895458087625118368,
           -1.46598909439115818325066466416,
           -0.484935707515497653046233483105,
           +0.484935707515497653046233483105,
           +1.46598909439115818325066466416,
           +2.48432584163895458087625118368,
           +3.58182348355192692277623675546,
           +4.85946282833231215015516494660))
  wu <- matrix(c(0.0000043106526307182867322209547262,
                 0.000758070934312217670069636036508,
                 0.0191115805007702856047383687629,
                 0.135483702980267735563431657727,
                 0.344642334932019042875028116518,
                 0.344642334932019042875028116518,
                 0.135483702980267735563431657727,
                 0.0191115805007702856047383687629,
                 0.000758070934312217670069636036508,
                 0.0000043106526307182867322209547262))
  #
  opt <- nlm(klobj_fixed,p=initial_vals,fixed_pars=fixed_pars,fixed_ix=fixed_ix,pars0=pars0,Xmat=Xmat,het=het,u=u,wu=wu)
  #
  return(opt$minimum)
}
#
# Find lower bound of interval M(\theta)
#
find_int_lower <- function(fixed_ix,pars0,step,data,het){
  #
  tol <- 1e-7
  #
  X  <- data$X   # regressors for outcome equation
  Z0 <- data$Z0  # Z for selection equation (no trade)
  Xmat <- rbind(X,Z0)
  #
  # initialize
  lower_temp <- pars0[fixed_ix]
  initial_vals <- pars0[-fixed_ix]
  eps <- 0
  #
  while(step>0.001){
    #
    eps <- kldist_fixed(lower_temp-step,fixed_ix,pars0,initial_vals,Xmat,het)
    if(eps > tol){
      step <- step/2
    }else{
      lower_temp <- lower_temp - step
    }
  }
  return(lower_temp)
}
#
# Find upper bound of interval M(\theta)
#
find_int_upper <- function(fixed_ix,pars0,step,data,het){
  #
  tol <- 1e-7
  #
  X  <- data$X   # regressors for outcome equation
  Z0 <- data$Z0  # Z for selection equation (no trade)
  Xmat <- rbind(X,Z0)
  #
  # initialize
  upper_temp <- pars0[fixed_ix]
  initial_vals <- pars0[-fixed_ix]
  eps <- 0
  #
  while(step>0.001){
    #
    eps <- kldist_fixed(upper_temp+step,fixed_ix,pars0,initial_vals,Xmat,het)
    if(eps > tol){
      step <- step/2
    }else{
      upper_temp <- upper_temp + step
    }
  }
  #
  return(upper_temp)
}
#
# END



