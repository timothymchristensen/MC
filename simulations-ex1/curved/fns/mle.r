#
# Unconstrained MLE
#
mle <- function(initial_vals,mdata,loglike=NULL,transform=TRUE){
    #
    if(transform==TRUE){
        initial_vals <- log(initial_vals/(1-initial_vals))
    }
    #
    if(is.null(loglike)==TRUE){
        stop("please provide a log-likelihood function!")
    }
    #
    opt <- optim(initial_vals,loglike,mdata=mdata,transformed=transform,method="Nelder-Mead",control=list(maxit=10000),hessian=FALSE)
    #
    MLE <- 0
    if(transform==TRUE){
        MLE <- exp(opt$par)/(1+exp(opt$par))
    }else{
        MLE <- opt$par
    }
    #
    return(MLE)
}
#
# MLE for theta using analytical solution for gamma MLE
#
mle_gamma_analytical <- function(initial_vals,mdata,loglike=NULL,transform=TRUE){
    #
    if(transform==TRUE){
        initial_vals <- log(initial_vals/(1-initial_vals))
    }
    #
    if(is.null(loglike)==TRUE){
        stop("please provide a log-likelihood function!")
    }
    #
    opt <- optim(initial_vals,loglike,mdata=mdata,transformed=transform,method="Nelder-Mead",control=list(maxit=10000),hessian=FALSE)
    #
    MLE <- 0
    if(transform==TRUE){
        MLE <- exp(opt$par)/(1+exp(opt$par))
    }else{
        MLE <- opt$par
    }
    #
    D  <- mdata[,1]
    #
    pi_00 <- sum(1-D)/n
    gamma <- 1 - pi_00
    #
    MLE <- c(MLE,gamma)
    #
    return(MLE)
}
#
# Profile MLE for (beta,gamma) when mu is fixed
#
mle_mu_fixed <- function(mu,mdata,opt_method=1){
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
#
# Identical to 'mle_mu_fixed', but returns an indicator of whether constrained optimization was necessary
#
mle_mu_fixed2 <- function(mu,mdata,opt_method=1){
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
    use_constrain <- 0
    #
    opt <- 0; MLE <- c(beta,gamma)
    if(beta > 1 || beta < 0 || mu - beta*(1-gamma) < 0 || mu - beta*(1-gamma) > gamma){
        #
        use_constrain <- 1
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
    MLE <- list(MLE=MLE,use_constrain=use_constrain)
    #
    return(MLE)
}
#
#END