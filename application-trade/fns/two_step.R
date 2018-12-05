#
# Fit via two-step procedure to initialize
#
two_step <- function(data,het=FALSE){
  #
  Y1 <- data$Y1  # dep var for outcome equation
  X  <- data$X   # regressors for outcome equation
  Z0 <- data$Z0  # Z for selection equation (no trade)
  #
  # Step 1: heretoskedastic probit
  if(het){
    initial_vals_probit <- numeric(23)
  }else{
    initial_vals_probit <- numeric(21)
  }
  Y_probit <- rbind(matrix(1,nrow(data$X),1),matrix(0,nrow(data$Z0),1))
  X_probit <- rbind(X,Z0)
  opt_het_probit <- optim(initial_vals_probit,loglike_het_probit,X=X_probit,Y=Y_probit,het=het,method="Nelder-Mead",control=list(maxit=10000),hessian=FALSE)
  opt_het_probit <- optim(opt_het_probit$par,loglike_het_probit,X=X_probit,Y=Y_probit,het=het,method="BFGS",control=list(maxit=10000),hessian=FALSE)
  opt_het_probit <- optim(opt_het_probit$par,loglike_het_probit,X=X_probit,Y=Y_probit,het=het,method="Nelder-Mead",control=list(maxit=10000),hessian=FALSE)
  probit_pars <- matrix(as.numeric(opt_het_probit$par))
  beta_z <- probit_pars[1:21]
  #
  # Step 2: outcome equation
  if(het){
    varpi <- probit_pars[22:23]
    sigma_z <- exp(varpi[1]*X[,2] + varpi[2]*X[,2]^2)
    inv_mills <- dnorm((X%*%probit_pars[1:21])/sigma_z) / pnorm((X%*%probit_pars[1:21])/sigma_z)
  }else{
    sigma_z <- 1
    inv_mills <- dnorm(X%*%probit_pars[1:21]) / pnorm(X%*%probit_pars[1:21])
  }
  #
  X_reg <- cbind(X,inv_mills)
  second_step_pars <- solve(t(X_reg)%*%X_reg)%*%t(X_reg)%*%Y1
  delta <- second_step_pars[11]/beta_z[11]
  beta_m <- second_step_pars[-22] - delta*beta_z
  beta_m <- beta_m[-11]
  #
  rho_sigma_m <- second_step_pars[22]
  #
  e_Y1 <- Y1 - X_reg%*%second_step_pars
  #
  term_1 <- var(e_Y1)
  term_2 <- rho_sigma_m^2*( mean(inv_mills^2 + inv_mills*(X%*%probit_pars[1:21])/sigma_z))
  sigma_m <- sqrt(term_1 + term_2)
  #
  rho <- rho_sigma_m/sigma_m
  #
  if(het){
    ret <- c(beta_m,beta_z,delta,log(sigma_m),log((1+rho)/(1-rho)),varpi)
  }else{
    ret <- c(beta_m,beta_z,delta,log(sigma_m),log((1+rho)/(1-rho)))
  }
  #
  return(ret)
}
#
# Heteroskedastic probit likelihood
#
loglike_het_probit <- function(het_probit_pars,X,Y,het=FALSE){
  #
  beta_z <- matrix(het_probit_pars[1:21])
  #
  if(het){
    varpi <- het_probit_pars[22:23]
    sigma_z0 <- exp(varpi[1]*X[,2] + varpi[2]*X[,2]^2)
  }else{
    sigma_z0 <- 1
  }
  #
  p <- .Internal(pnorm(X%*%beta_z / sigma_z0,0,1,TRUE,FALSE))
  #
  ret <- sum( Y*log(p) + (1-Y)*log(1-p))
  #
  # return negative log likelihood
  if(is.na(ret)==TRUE||is.nan(ret)==TRUE||is.infinite(ret)==TRUE){
    return(1E10)
  }else{
    return(-ret)
  }
}
#
# END
