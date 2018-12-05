#
# Generate data for Example 2
#
gen_data <- function(pars,n=10000){
    #
    Delta1 <- pars[1]
    Delta2 <- pars[2]
    beta1  <- pars[3]
    beta2  <- pars[4]
    rho    <- pars[5]
    s      <- pars[6]
    #
    pi_vec <- par_to_prob(c(Delta1,Delta2,beta1,beta2,rho,s))
    #
    pi_check <- c(0,cumsum(pi_vec))
    D_unif <- runif(n)
    #
    D_mat <- t(apply(matrix(D_unif),1,function(x) (x <= pi_check[2:5])*(x > pi_check[1:4])))
    #
    colnames(D_mat) <- c("D00", "D01", "D10", "D11")
    #
    return(D_mat)
}
#
# END