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
    mdata <- cbind(D,Y*D)
    #
    return(mdata)
}
#
#END