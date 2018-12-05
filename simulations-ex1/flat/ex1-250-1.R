#
# Example 1: Missing Data via SMC algorithm and likelihood criterion function
#
rm(list=ls())
#
n <- 250  # sample size
N <- 5000  # number of replications
#
source("fns/ex1_ll.r")
source("fns/gen_data.R")
source("fns/mle.r")
source("fns/likelihood_ratios.R")
source("fns/likelihood_ratios_weighted.R")
source("fns/ex1_p2_par.R")
source("fns/SMC.R")
source("fns/ex1_gms.R")
require("parallel")
#
set.seed(1988)
#
trans_pars <- TRUE    # whether to transform the parameters to have unbounded support
#
mu    <- 0.50    # Pr(Y_i = 1)
beta  <- 0.50    # Pr(Y_i = 1 | D_i = 0)
gamma <- 1-1/sqrt(n)    # Pr(D_i = 1); rho in the paper
#
# grid over M_I to check coverage
if(gamma < 1){
  mu_seq_n <- 200  
  M_I <- seq(mu - beta*(1-gamma), mu +(1-beta)*(1-gamma),length.out=mu_seq_n)
}else{
  mu_seq_n <- 1
  M_I <- mu
}
#
# true parameters
pars_0 <- c(mu,beta,gamma)
npar <- length(pars_0)
#
initial_vals <- c(0.45,0.45,0.2) # initial values for MLE, if needed
#
scalepar <- list(Nphi=200,Nmh=1,lambda=2,scale_cov=1)
keep <- 10000   # number of draws
#
alpha_vec <- c(0.9,0.95,0.99)  # significance levels
#
# pre-allocate for acceptance ratio
accept <- numeric(N)
#
# pre-allocate for coverage checking and CS length
p1_full_check <- p1_proj_check <- p2_check  <- p3_check <- pc_check <- mi_check <- matrix(0,nrow=N,ncol=length(alpha_vec))
p1_lower <- p2_lower <- p3_lower <- pc_lower <- mi_lower <- matrix(0,nrow=N,ncol=length(alpha_vec))
p1_upper <- p2_upper <- p3_upper <- pc_upper <- mi_upper <- matrix(0,nrow=N,ncol=length(alpha_vec))
p1_length <- p2_length <- p3_length <- pc_length <- mi_length <- matrix(0,nrow=N,ncol=length(alpha_vec))
Xi_1 <- Xi_2 <- Xi_3 <- matrix(0,nrow=N,ncol=length(alpha_vec))
#
# set up parallelization 
cl <- makeCluster(getOption("cl.cores", 10),type="PSOCK")
#
for(boot_iter in 1:N){
  #
  mdata <- gen_data(pars_0,n)
  mdata_short <- c(colSums(mdata),n)
  #
  pars_mle <- mle_gamma_analytical(initial_vals[-3],mdata,loglike_gamma_analytical,trans_pars)
  #
  smc_run <- SMC(mdata_short,trans_pars,keep,scalepar,loglike_short,PKernel_temp,cl)
  #
  accept[boot_iter] <- smc_run$accept_ratio
  #
  # critical values
  # procedure 1: full set
  xi_1 <- p1_crit_weighted(pars_mle,mcmc_draws=smc_run$Draws,mdata,weights=smc_run$weights,alpha_vec)
  p1_full_check[boot_iter,] <- full_check(pars_mle,pars_0,mdata,xi=xi_1)
  #
  # procedure 1: projection for mu
  p1_proj_check[boot_iter,] <- subvec_mu_cvg(pars_mle,M_I,mdata,xi=xi_1)
  subvec_out <- subvec_mu_length(pars_mle,mdata,xi=xi_1)
  p1_lower[boot_iter,] <- subvec_out[1,]
  p1_upper[boot_iter,] <- subvec_out[2,]
  p1_length[boot_iter,] <- subvec_out[3,]
  #
  # store critical value
  Xi_1[boot_iter,] <- xi_1
  #
  # procedure 2: profiling for mu
  xi_2 <- p2_crit_par_weighted(pars_mle,smc_run$Draws,mdata,weights=smc_run$weights,alpha_vec,cl)
  #
  p2_check[boot_iter,] <- subvec_mu_cvg(pars_mle,M_I,mdata,xi=xi_2)
  subvec_out <- subvec_mu_length(pars_mle,mdata,xi=xi_2)
  p2_lower[boot_iter,] <- subvec_out[1,]
  p2_upper[boot_iter,] <- subvec_out[2,]
  p2_length[boot_iter,] <- subvec_out[3,]
  #
  # store critical value
  Xi_2[boot_iter,] <- xi_2
  #
  # procedure 3: chi-square 1 upper bound
  xi_3 <- qchisq(alpha_vec,df=1)/2
  p3_check[boot_iter,] <- subvec_mu_cvg(pars_mle,M_I,mdata,xi=xi_3)
  subvec_out <- subvec_mu_length(pars_mle,mdata,xi=xi_3)
  p3_lower[boot_iter,] <- subvec_out[1,]
  p3_upper[boot_iter,] <- subvec_out[2,]
  p3_length[boot_iter,] <- subvec_out[3,]
  #
  # store critical value
  Xi_3[boot_iter,] <- xi_3
  #
  # percentile
  pc_out <- percentile_CS_weighted(smc_run$Draws[,1],weights=smc_run$weights,alpha_vec,M_I)
  pc_check[boot_iter,] <- pc_out$cvg_check
  pc_lower[boot_iter,] <- pc_out$lower
  pc_upper[boot_iter,] <- pc_out$upper
  pc_length[boot_iter,] <- pc_out$length
  #
  # comparison with moment inequalities
  GMS_out <- GMS_mu(mdata,alpha_vec)
  mi_check[boot_iter,] <- mi_cvg_check(M_I,GMS_out$lower,GMS_out$upper)
  mi_lower[boot_iter,] <- GMS_out$lower
  mi_upper[boot_iter,] <- GMS_out$upper
  mi_length[boot_iter,] <- GMS_out$length
  #
  # print outputs
  if(boot_iter%%10==0){
    cat("Iteration: ",boot_iter,". Time: ",as.character(Sys.time()),". Chain acceptance: ",sum(accept)/boot_iter,".\n",sep="")
  }
  if(boot_iter%%100==0){
    print_p1_full <- colSums(p1_full_check[1:boot_iter,])/boot_iter
    cat("Coverage p1 full: (",print_p1_full[1],",",print_p1_full[2],",",print_p1_full[3],").\n",sep="")
    print_p1_proj <- colSums(p1_proj_check[1:boot_iter,])/boot_iter
    cat("Coverage p1 proj: (",print_p1_proj[1],",",print_p1_proj[2],",",print_p1_proj[3],").\n",sep="")
    print_p2 <- colSums(p2_check[1:boot_iter,])/boot_iter
    cat("Coverage p2:      (",print_p2[1],",",print_p2[2],",",print_p2[3],").\n",sep="")
    print_p3 <- colSums(p3_check[1:boot_iter,])/boot_iter
    cat("Coverage p3:      (",print_p3[1],",",print_p3[2],",",print_p3[3],").\n",sep="")
    print_mi <- colSums(mi_check[1:boot_iter,])/boot_iter
    cat("Coverage GM:      (",print_mi[1],",",print_mi[2],",",print_mi[3],").\n",sep="")
    
  }
}
boot_results <- list(p1_full_check=p1_full_check,
                     p1_proj_check=p1_proj_check,
                     p2_check=p2_check,
                     p3_check=p3_check,
                     pc_check=pc_check,
                     mi_check=mi_check,
                     p1_lower=p1_lower,
                     p2_lower=p2_lower,
                     p3_lower=p3_lower,
                     pc_lower=pc_lower,
                     mi_lower=mi_lower,
                     p1_upper=p1_upper,
                     p2_upper=p2_upper,
                     p3_upper=p3_upper,
                     pc_upper=pc_upper,
                     mi_upper=mi_upper,
                     p1_length=p1_length,
                     p2_length=p2_length,
                     p3_length=p3_length,
                     pc_length=pc_length,
                     mi_length=mi_length,
                     Xi_1=Xi_1,Xi_2=Xi_2,Xi_3=Xi_3,n=n,N=N,pars_0=pars_0)
save(boot_results,file="results/n250/1.RData")
#
stopCluster(cl)
#
#END
