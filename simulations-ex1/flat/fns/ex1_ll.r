#
# Log-likelihood functions for Example 1
#
loglike <- function(pars,mdata,transformed=TRUE){
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
  n <- nrow(mdata)
  #
  D  <- mdata[,1]
  YD <- mdata[,2]
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
    Term3 <- (1 - D)*log(1 - gamma)
  }
  # sum the loglikelihood terms
  ret <- sum(Term1 + Term2 + Term3)
  #
  # return MINUS the log-likelihood (as we're running a minimization routine)
  return(-ret)
}
#
# As above, but using known solution for MLE for gamma 
#
loglike_gamma_analytical <- function(pars,mdata,transformed=TRUE){
    #
    n <- nrow(mdata)
    #
    mu <- beta <- gamma <- 0
    #
    D  <- mdata[,1]
    YD <- mdata[,2] # Y*D
    #
    pi_00 <- sum(1-D)/n
    gamma <- 1 - pi_00
    #
    if(transformed==TRUE){
        mu <- exp(pars[1])/(1+exp(pars[1]))     # Pr(Y_i = 1)
        beta  <- exp(pars[2])/(1+exp(pars[2]))  # Pr(Y_i = 1 | D_i = 0)
    }else{
        mu <- pars[1]     # Pr(Y_i = 1)
        beta  <- pars[2]  # Pr(Y_i = 1 | D_i = 0)
    }
    # for convenience, define
    rho <- mu - beta*(1-gamma)
    #
    Term1 <- YD*log(rho)
    Term2 <- (D - YD)*log(gamma - rho)
    #
    Term3 <- 0
    if(sum(D)!=n){
        Term3 <- (1 - D)*log(1 - gamma)
    }
    # sum the loglikelihood terms
    ret <- sum(Term1 + Term2 + Term3)
    # return MINUS the log-likelihood (as we're running a minimization routine)
    return(-ret)
}
#
# LL constrained problem when we breach identified set constraints
#
loglike_constr <- function(pars,mu,mdata,transformed=TRUE){
    #
    beta <- gamma <- 0
    #
    if(transformed==TRUE){
        beta  <- exp(pars[1])/(1+exp(pars[1]))
        gamma <- exp(pars[2])/(1+exp(pars[2]))
    }else{
        beta  <- pars[1]
        gamma <- pars[2]
    }
    # for convenience, define
    rho <- mu - beta*(1-gamma)
    #
    ret <- -10000000
    #
    if(beta <= 1 && beta >= 0 && rho<=gamma && rho >= 0){
        #
        n <- nrow(mdata)
        #
        D  <- mdata[,1]
        YD <- mdata[,2] # Y*D
        #
        Term1 <- YD*log(rho)
        Term2 <- (D - YD)*log(gamma - rho)
        #
        Term3 <- 0  #tim_normalization
        if(sum(D)!=n){
            Term3 <- (1 - D)*log(1 - gamma)
        }
        # sum the loglikelihood terms
        ret_temp <- sum(Term1 + Term2 + Term3)
        #
        # check
        if(is.nan(ret_temp)==TRUE || is.infinite(abs(ret_temp))==TRUE){
          ret <- ret
        }else{
          ret <- ret_temp
        }
    }
    # return MINUS the log-likelihood (as we're running a minimization routine)
    return(-ret)
}
#
#
loglike_constr_2 <- function(pars,mu,mdata){
    #
    beta <- gamma <- 0
    #
    beta  <- pars[1]
    gamma <- pars[2]
    # for convenience, define
    rho <- mu - beta*(1-gamma)
    #
    ret <- -10000000
    #
    if(beta <= 1 && beta >= 0 && rho<=gamma && rho >= 0){
        #
        n <- nrow(mdata)
        #
        D  <- mdata[,1]
        YD <- mdata[,2] # Y*D
        #
        Term1 <- YD*log(rho)
        Term2 <- (D - YD)*log(gamma - rho)
        #
        Term3 <- 0
        if(sum(D)!=n){
            Term3 <- (1 - D)*log(1 - gamma)
        }
        # sum the loglikelihood terms
        ret <- sum(Term1 + Term2 + Term3)
    }
    # return MINUS the log-likelihood (as we're running a minimization routine)
    return(-ret)
}
#
# Log-likelihood at unconstrained MLE
#
loglike_p_hat <- function(mdata){
    #
    n <- nrow(mdata)
    #
    D  <- mdata[,1]
    YD <- mdata[,2] # Y*D
    #
    pi_11 <- sum(YD)/n
    pi_10 <- sum(D-YD)/n
    pi_00 <- sum(1-D)/n
    #
    Term1 <- YD*log(pi_11)
    Term2 <- (D - YD)*log(pi_10)
    #
    Term3 <- 0
    if(pi_00!=0){
        Term3 <- (1 - D)*log(pi_00)
    }
    #
    return(sum(Term1 + Term2 + Term3))
}
# 
#END