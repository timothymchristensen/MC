#
# Convert parameters to reduced-form probabilities for Example 2
#
par_to_prob <- function(pars,ssame){
  #
  Delta1 <- pars[1]
  Delta2 <- pars[2]
  beta1  <- pars[3:5]
  beta2  <- pars[6:8]
  rho    <- as.double(pars[9])
  if(ssame){
    s <- rep(pars[10],8)
  }else{
    s <- pars[10:17]
  }
  #
  lt <- as.integer(1)
  prob <- double(lt)
  #
  b1x <- beta1[1]+c(0,0,beta1[3],beta1[2],beta1[3],beta1[2],beta1[2]+beta1[3],beta1[2]+beta1[3])
  b2x <- beta2[1]+c(0,beta2[3],0,beta2[2],beta2[3],beta2[2]+beta2[3],beta2[2],beta2[2]+beta2[3])
  #
  pi_00 <- pi_01_Term1_a <- pi_01_Term1_b <- pi_01_Term3_a <- numeric(8)
  #
  for(i in 1:8){
    pi_00[i] <- pbivnormu(-b1x[i],-b2x[i],rho,lt,prob)
    pi_01_Term1_a[i] <- pbivnormu(-b1x[i]-Delta1,-b2x[i]-Delta2,rho,lt,prob)
    pi_01_Term1_b[i] <- pbivnormu(-b1x[i],-b2x[i]-Delta2,rho,lt,prob)
    pi_01_Term3_a[i] <- pbivnorml(-b1x[i],-b2x[i],-b1x[i]-Delta1,-b2x[i]-Delta2,rho,lt,prob)
  }
  #
  pi_01_Term1 <- pnorm(-b1x-Delta1) - pnorm(-b1x) - pi_01_Term1_a + pi_01_Term1_b
  pi_01_Term2 <- pnorm(-b1x) - pi_00   
  pi_01_Term3 <- (1-s)*pi_01_Term3_a
  pi_01 <- pi_01_Term1 + pi_01_Term2 + pi_01_Term3
  #
  pi_11 <- 1 - pnorm(-b1x-Delta1) - pnorm(-b2x-Delta2) + pi_01_Term1_a
  #
  pi_10 <- 1 - pi_00 - pi_01 - pi_11
  #
  ret <- c(pi_00,pi_01,pi_10,pi_11)
  #
  return(ret)
}
#
# Slightly edited pbivnorm function to reduce computations
#
pbivnormu <- function(x, y, rho,lt,prob) {
  #
  ## coerce arguments to proper types to be passed to fortran
  lower <- as.double(c(0, 0))
  infin <- as.integer(c(0, 0))
  uppera <- as.double(x)
  upperb <- as.double(y)
  
  ans <- .Fortran("PBIVNORM", prob, lower, uppera, upperb, infin, rho, lt,
                  PACKAGE="pbivnorm")[[1]]
  return(ans)
}
#
# Slightly edited pbivnorm function to reduce computations, with lower limit
#
pbivnorml <- function(xl, yl, xu, yu, rho,lt,prob) {
  #
  ## coerce arguments to proper types to be passed to fortran
  lower <- as.double(c(xl, yl))
  infin <- as.integer(c(2, 2))
  uppera <- as.double(xu)
  upperb <- as.double(yu)
  
  ans <- .Fortran("PBIVNORM", prob, lower, uppera, upperb, infin, rho, lt,
                  PACKAGE="pbivnorm")[[1]]
  return(ans)
}
#
# END