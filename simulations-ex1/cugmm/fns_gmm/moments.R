#
# Moment conditions for Example 1
#
moments_raw <- function(pars,mdata_short){
    #
    mu    <- pars[1]
    beta  <- pars[2]
    gamma <- pars[3]
    #
    check_par_space <- enforce_restrictions(pars,FALSE)
    #
    mmnts <- 0
    #
    if(check_par_space==1){# restrictions on the parameter space are respected
        YD <- mdata_short[1]
        D  <- mdata_short[2]
        #
        mmnts <- cbind(YD - mu + beta*(1 - gamma),D - gamma)
    }else{
        mmnts <- cbind(1e8,1e8)
    }
    #
    return(mmnts)
}
#
enforce_restrictions <- function(pars,transformed){
  #
  if(transformed==TRUE){
    pars <- exp(pars)/(1+exp(pars))
  }
  #
  mu    <- pars[1]    # Pr(Y_i = 1)
  beta  <- pars[2]    # Pr(Y_i = 1 | D_i = 0)
  gamma <- pars[3]    # Pr(D_i = 1)
  #
  ret <- 1
  #
  check_term <- mu - beta * (1 - gamma)
  #
  if(check_term < 0){
    ret <- 0
  }
  #
  if(check_term > gamma){
    ret <- 0
  }
  #
  return(ret) # = 1 if theta lies in the parameter space contraints, 0 otherwise
}
#
#END