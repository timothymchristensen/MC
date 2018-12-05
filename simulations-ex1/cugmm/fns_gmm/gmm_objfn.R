#
# GMM objective function
#
gmm_objfn <- function(pars,mdata_short,W_mat,transformed=TRUE){
  #
  if(transformed==TRUE){
    pars <- exp(pars)/(1+exp(pars))
  }
  pars[is.nan(pars)] <- 1
  #
  check_par_space <- enforce_restrictions(pars,FALSE)
  #
  ret <- 1e8
  #
  if(check_par_space==1){ # restrictions on the parameter space are respected
    #
    mu    <- pars[1]
    beta  <- pars[2]
    gamma <- pars[3]
    #
    if(mdata_short[2]==1&gamma==1){
      #
      YD <- mdata_short[1]
      g_vec <- YD - mu + beta*(1 - gamma)
      ret <- 0.5*g_vec^2*W_mat[1,1]
    }
    #
    if(mdata_short[2]<1){
      #
      YD <- mdata_short[1]
      D  <- mdata_short[2]
      #
      g_vec <- rbind(YD - mu + beta*(1 - gamma),D - gamma)
      #
      ret <- 0.5*t(g_vec)%*%W_mat%*%g_vec
    }
  }
  #
  return(ret)
}
#
# Profile GMM objective function for fixed mu
#
gmm_objfn_profile <- function(pars,mu,mdata_short,W_mat,transformed=TRUE){
  #
  if(transformed==TRUE){
    pars <- exp(pars)/(1+exp(pars))
  }
  if(is.nan(pars[2])==TRUE){pars[2] <- 1}
  #
  pars <- c(mu,pars)
  #
  check_par_space <- enforce_restrictions(pars,FALSE)
  #
  ret <- 1e8
  #
  if(check_par_space==1){ # restrictions on the parameter space are respected
    ret <- as.numeric(gmm_objfn(pars,mdata_short,W_mat,FALSE))
  }
  #
  return(ret)
}
#
#END