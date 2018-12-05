#
# Generate data for Example 1
#
gen_data <- function(pars,n=10000){
    #
    mu <- pars[1]      # Pr(Y_i = 1)
    beta  <- pars[2]   # Pr(Y_i = 1 | D_i = 0)
    gamma <- pars[3]   # Pr(D_i = 1)
    #
    delta <- (mu - beta*(1-gamma))/gamma # Pr(Y_i = 1 | D_i = 1)
    #
    D <- rbinom(n,1,gamma)
    #
    Y <- numeric(n)
    #
    for(i in 1:n){
        if(D[i]==1){
            Y[i] <- rbinom(1,1,delta)
        }else if(D[i]==0){
            Y[i] <- rbinom(1,1,beta)
        }
    }
    #
    mdata <- cbind(Y*D,D)
    #
    return(mdata)
}
#
# Transform data for Example 1
#
data_trans <- function(mdata){
  #
  mvec <- colMeans(mdata)
  mdata_short <- as.numeric(c(mvec,nrow(mdata)))
  cvmm <- cov(mdata)
  #
  W <- matrix(NA,2,2)
  W[1,1] <- mvec[2]/(mvec[1]*(mvec[2]-mvec[1]))
  W[1,2] <- W[2,1] <- -1/(mvec[2]-mvec[1])
  if(mvec[2]==1){
    W[2,2] <- Inf
  }else{
    W[2,2] <- (1-mvec[1])/((1-mvec[2])*(mvec[2]-mvec[1]))
  }
  #
  return(list(mdata_short=mdata_short,W_mat=W))
}
#
#END