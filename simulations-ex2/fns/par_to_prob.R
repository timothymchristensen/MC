#
# Convert parameters to reduced-form probabilities for Example 2
#
par_to_prob <- function(pars){
  #
  Delta1 <- pars[1]
  Delta2 <- pars[2]
  beta1  <- pars[3]
  beta2  <- pars[4]
  rho    <- as.double(pars[5])
  s      <- pars[6]
  #
  lt <- as.integer(1)
  prob <- double(lt)
  #
  pi_00 <- pbivnormu(-beta1,-beta2,rho,lt,prob)
  #
  pi_01_Term1_a <- pbivnormu(-beta1 - Delta1,-beta2-Delta2,rho,lt,prob)
  #
  pi_01_Term1 <- pnorm(-beta1-Delta1) - pnorm(-beta1) - pi_01_Term1_a + pbivnormu(-beta1,-beta2-Delta2,rho,lt,prob)
  pi_01_Term2 <- pnorm(-beta1) - pi_00   
  pi_01_Term3 <- s*pbivnorml(-beta1,-beta2,-beta1-Delta1,-beta2-Delta2,rho,lt,prob)
  pi_01 <- pi_01_Term1 + pi_01_Term2 + pi_01_Term3
  #
  pi_11 <- 1 - pnorm(-beta1-Delta1) - pnorm(-beta2-Delta2) + pi_01_Term1_a
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