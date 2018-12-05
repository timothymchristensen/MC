#
# GMM estimator
#
gmm_argmin <- function(initial_vals,mdata_short,W_mat){
  #
  # analytic solution under point identification
  if(mdata_short[2]==1){
    theta_hat <- c(mdata_short[1],0.5,mdata_short[2])
  }else{
    initial_vals <- log(initial_vals/(1-initial_vals))
    #
    opt1 <- optim(initial_vals,gmm_objfn,mdata_short=mdata_short,W_mat=W_mat,transformed=TRUE,method="Nelder-Mead",control=list(maxit=10000,reltol=1e-12),hessian=FALSE)
    opt <- optim(opt1$par,gmm_objfn,mdata_short=mdata_short,W_mat=W_mat,transformed=TRUE,method="Nelder-Mead",control=list(maxit=10000,reltol=1e-12),hessian=FALSE)
    #
    theta_hat <- exp(opt$par)/(1+exp(opt$par))
    theta_hat[is.nan(theta_hat)] <- 1
  }
  #
  return(theta_hat)
}
#
# Profile GMM estimator for fixed mu 
#
gmm_mu_subvec <- function(mu,mdata_short,W_mat){
    #
    # initial values
    gamma <- mdata_short[2]
    #
    beta <- 0.5
    if(gamma!=1){
        beta <- (mu - mdata_short[1])/(1-gamma)
    }else{
      gamma <- 0.999
    }
    #
    initial_vals <- c(beta,gamma)
    #
    if(beta >= 1 || beta <= 0 || mu - beta*(1-gamma) <= 0 || mu - beta*(1-gamma) >= gamma){
      initial_vals[1] <- mu
    }
    #
    initial_vals <- log(initial_vals/(1-initial_vals))
    #
    opt1 <- optim(initial_vals,gmm_objfn_profile,mu=mu,mdata_short=mdata_short,W_mat=W_mat,transformed=TRUE,method="Nelder-Mead",control=list(maxit=10000),hessian=FALSE)
    opt <- optim(opt1$par,gmm_objfn_profile,mu=mu,mdata_short=mdata_short,W_mat=W_mat,transformed=TRUE,method="Nelder-Mead",control=list(maxit=10000),hessian=FALSE)
    #
    theta_hat_cons <- exp(opt$par)/(1+exp(opt$par))
    #
    return(theta_hat_cons)
}
#
#END