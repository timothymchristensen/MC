# 
# Construct CS for mu based on GMS test
#
GMS_mu <- function(mdata,alpha_vec){
  #
  n <- nrow(mdata)
  kappa <- sqrt(log(n))
  # 
  # Do testing at each point on a grid
  mgrid <- seq(0,1,length.out=1001)
  #
  nboot <- 1000
  MR_boot <- matrix(NA,nrow=nboot,ncol=length(mgrid))
  #
  # pre-calculate moment functions
  mn1 <- mdata[,2]
  mn2 <- mdata[,2]+(1-mdata[,1])
  m1 <- mean(mn1)
  m2 <- mean(mn2)
  s1 <- sd(mn1) 
  s2 <- sd(mn2)
  #
  # pre-calculate thresholds
  l1 <- sqrt(n)*(mgrid-m1)/(s1*kappa)
  l2 <- sqrt(n)*(m2-mgrid)/(s2*kappa)
  #
  # criterion function
  Tn <- numeric(length(mgrid))
  for(i in 1:length(mgrid)){
    Tn[i] <- n*MI_criterion(mgrid[i],m1,m2,s1,s2)
  }
  #
  T_boot <- matrix(NA,nrow=nboot,ncol=length(mgrid))
  #
  for(i in 1:nboot){
    #
    zeta <- rnorm(n)
    nu1 <- sum(zeta*(m1-mn1))/(sqrt(n)*s1)
    nu2 <- sum(zeta*(mn2-m2))/(sqrt(n)*s2)
    #
    T_boot[i,] <-  (min(nu1,0)^2)*(l1 <= 1) + (min(nu2,0)^2)*(l2 <= 1)
  }
  #
  # Determine whether or not in CS at level alpha
  CI_mat <- matrix(NA,nrow=length(alpha_vec),ncol=length(mgrid))
  lower <- upper <- numeric(length(mgrid))
  #
  for(i in 1:length(mgrid)){
    xi_crit <- quantile(T_boot[,i],alpha_vec)
    CI_mat[,i] <- matrix(Tn[i] <= xi_crit)*1
  }
  #
  # return CS
  lower <- upper <- numeric(length(alpha_vec))
  for(j in 1:length(alpha_vec)){
    lower[j] <- min(mgrid[as.logical(CI_mat[j,])])
    upper[j] <- max(mgrid[as.logical(CI_mat[j,])])
  }
  #
  return(list(lower=lower,upper=upper,length=(upper-lower)))
  #
}
# 
# Construct CS for mu based on moment inequalities using Bugni, Canay, Shi minimum resampling test
#
MR_mu <- function(mdata,alpha_vec){
  #
  n <- nrow(mdata)
  kappa <- sqrt(log(n))
  #
  # Do testing at each point on a grid
  mgrid <- seq(0,1,length.out=1001)
  #
  nboot <- 1000
  MR_boot <- matrix(NA,nrow=nboot,ncol=length(mgrid))
  #
  # pre-calculate moment functions
  mn1 <- mdata[,2]
  mn2 <- mdata[,2]+(1-mdata[,1])
  m1 <- mean(mn1)
  m2 <- mean(mn2)
  s1 <- sd(mn1) 
  s2 <- sd(mn2)
  #
  # pre-calculate thresholds
  l1 <- sqrt(n)*(mgrid-m1)/(s1*kappa)
  l2 <- sqrt(n)*(m2-mgrid)/(s2*kappa)
  #
  # criterion function
  Tn <- numeric(length(mgrid))
  for(i in 1:length(mgrid)){
    Tn[i] <- n*MI_criterion(mgrid[i],m1,m2,s1,s2)
  }
  #
  T_boot <- matrix(NA,nrow=nboot,ncol=length(mgrid))
  #
  for(i in 1:nboot){
    #
    zeta <- rnorm(n)
    nu1 <- sum(zeta*(m1-mn1))/(sqrt(n)*s1)
    nu2 <- sum(zeta*(mn2-m2))/(sqrt(n)*s2)
    #
    T_DR <- (min(nu1,0)^2)*(l1 <= 1) + (min(nu2,0)^2)*(l2 <= 1)
    T_PR <- (pmin(nu1+l1,0))^2 + (pmin(nu2+l2,0))^2
    T_boot[i,] <- pmin(T_DR,T_PR)
  }
  #
  # Determine whether or not in CS at level alpha
  CI_mat <- matrix(NA,nrow=length(alpha_vec),ncol=length(mgrid))
  lower <- upper <- numeric(length(mgrid))
  #
  for(i in 1:length(mgrid)){
    xi_crit <- quantile(T_boot[,i],alpha_vec)
    CI_mat[,i] <- matrix(Tn[i] <= xi_crit)*1
  }
  #
  # return CS
  lower <- upper <- numeric(length(alpha_vec))
  for(j in 1:length(alpha_vec)){
    lower[j] <- min(mgrid[as.logical(CI_mat[j,])])
    upper[j] <- max(mgrid[as.logical(CI_mat[j,])])
  }
  #
  return(list(lower=lower,upper=upper,length=(upper-lower)))
  #
}
#
# Moment inequality criterion
#
MI_criterion <- function(mu,m1,m2,s1,s2){
  #
  t1 <- (mu-m1)/s1
  t2 <- (m2-mu)/s2
  #
  return(min(t1,0)^2+min(t2,0)^2)
  #
}
#
# Moment inequality converage check for full set
# 
mi_cvg_check <- function(M_I,lower,upper){
  #
  lower_lim <- min(M_I)
  upper_lim <- max(M_I)
  #
  return((lower <= lower_lim)*(upper >= upper_lim))
}
#
# END

