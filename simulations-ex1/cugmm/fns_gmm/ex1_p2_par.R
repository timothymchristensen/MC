#
# Calculate critical values for procedure 2
#
p2_crit_par <- function(pars_mle,mcmc_draws,mdata,alpha_vec,cl){
  #
  n_alpha <- length(alpha_vec)
  B <- nrow(mcmc_draws)
  profile_qlr <- numeric(B)
  #
  mle_LL <- -loglike(pars_mle,mdata,FALSE)
  #
  profile_qlr <- parRapply(cl,mcmc_draws,p2_crit_par_fn,mdata=mdata,mle_LL=mle_LL,mle_mu_fixed_par=mle_mu_fixed_par,loglike=loglike,loglike_constr=loglike_constr)
  #
  ret <- as.numeric(quantile(profile_qlr,alpha_vec))
  #
  return(ret)
}
#
# Inner part for profile QLR 
#
p2_crit_par_fn <- function(pars_temp,mdata,mle_LL,mle_mu_fixed_par,loglike,loglike_constr){
  #
  lower_bd  <- pars_temp[1]-pars_temp[2]*(1-pars_temp[3])
  upper_bd  <- pars_temp[1]+(1-pars_temp[2])*(1-pars_temp[3])
  #
  mle_lower <- mle_mu_fixed_par(lower_bd,mdata,loglike_constr)
  mle_upper <- mle_mu_fixed_par(upper_bd,mdata,loglike_constr)
  #
  qlr_lower <- mle_LL+loglike(c(lower_bd,mle_lower),mdata,FALSE)
  qlr_upper <- mle_LL+loglike(c(upper_bd,mle_upper),mdata,FALSE)
  #
  ret <- max(abs(c(qlr_lower,qlr_upper)))
  #
  return(ret)
}
#
mle_mu_fixed_par <- function(mu,mdata,loglike_constr,opt_method=1){
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
  gamma <- 1 - pi_00
  #
  beta <- 0.5
  if(gamma!=1){
    beta <- (mu - pi_11)/(1-gamma)
  }
  #
  opt <- 0; MLE <- c(beta,gamma)
  if(beta > 1 || beta < 0 || mu - beta*(1-gamma) < 0 || mu - beta*(1-gamma) > gamma){
    #
    initial_vals <- c(mu,gamma)
    #
    if(opt_method==1){
      initial_vals <- log(initial_vals/(1-initial_vals))
      opt <- optim(initial_vals,loglike_constr,mu=mu,mdata=mdata,transformed=TRUE,method="Nelder-Mead",control=list(maxit=10000),hessian=FALSE)
      MLE <- exp(opt$par)/(1+exp(opt$par))
    }
    #
    else{
      eval_g0 <- function(x,mu,mdata){
        beta <- x[1]
        gamma <- x[2]
        ret <- c(-mu + beta*(1-gamma),mu - beta*(1-gamma) - gamma)
        return(ret)
      }
      #
      opt <- nloptr(x0=initial_vals,eval_f=loglike_constr_2,
                    lb=c(0,0),ub=c(1,1),eval_g_ineq=eval_g0,
                    opts=list("algorithm"="NLOPT_LN_COBYLA","xtol_rel"=1.0e-8,"maxeval"=1000,
                              "tol_constraints_ineq"=c(0,0),"print_level"=2),
                    mu=mu,mdata=mdata)
      #
      MLE <- exp(opt$par)/(1+exp(opt$par))
    }
  }
  #
  return(MLE)
}
