#
# Example 2: Entry Game via SMC algorithm and likelihood criterion function
#
rm(list=ls())
#
n <- 500   # sample size
N <- 1000  # number of replications
#
source("fns/gen_data.R") 
source("fns/par_to_prob.R")
source("fns/mle.R")
source("fns/SMC.R")
source("fns/par_transform.R")
source("fns/likelihood_ratios.R")
source("fns/ID_set.R")
require("parallel")
require("pbivnorm")
#
set.seed(1988)
#
trans_pars <- TRUE    # whether to transform the parameters to have unbounded support
#
# parameters
Delta1 <- -0.5
Delta2 <- -0.5
beta1  <- 0.2
beta2  <- 0.2
rho    <- 0.5
s      <- 0.5
#
lower_bounds <- c(-2,-2,-1,-1,0,0)
upper_bounds <- c(0,0,2,2,1,1)
#
# true parameters
pars_0 <- c(Delta1,Delta2,beta1,beta2,rho,s)
#
initial_vals <- pars_0 # initial values for MLE, if needed
#
scalepar <- list(Nphi=200,Nmh=4,Nblocks=1,lambda=2,scale_cov=1)
keep <- 10000   # number of draws
#
alpha_vec <- c(0.9,0.95,0.99)  # significance levels
#
# pre-allocate for acceptance ratio
accept <- numeric(N)
#
# pre-allocate for coverage checking and CS length
p1_full_check <- matrix(0,nrow=N,ncol=length(alpha_vec))
p1_proj_check <- p2_check <- p3_check <- pc_check <- matrix(0,nrow=N,ncol=2*length(alpha_vec))
p1_lower <- p2_lower <- p3_lower <- pc_lower <- matrix(0,nrow=N,ncol=2*length(alpha_vec))
p1_upper <- p2_upper <- p3_upper <- pc_upper <- matrix(0,nrow=N,ncol=2*length(alpha_vec))
p1_length <- p2_length <- p3_length <- pc_length <- matrix(0,nrow=N,ncol=2*length(alpha_vec))
Xi_1 <- Xi_2_Delta_1 <- Xi_2_beta_1 <- Xi_3 <- matrix(0,nrow=N,ncol=length(alpha_vec))
#
# set up parallelization 
cl <- makeCluster(getOption("cl.cores", 20),type="PSOCK")
clusterExport(cl=cl, list("par_trans","par_inv_trans","par_to_prob","pbivnormu","pbivnorml","pbivnorm","PKernel_temp","loglike_short","loglike_short_true","find_Delta_1_lower","find_Delta_1_upper","find_beta_1_lower","find_beta_1_upper","kldist_fixed_raw","klobj","klobj_fixed","profile_LL_Delta_1","mle_Delta_1_fixed","loglike_Delta_1_fixed","profile_LL_beta_1","mle_beta_1_fixed","loglike_beta_1_fixed"))
#
M_I_Delta_1 <- seq(-1.42490,-0.000001,length.out=200)
M_I_beta_1  <- seq(-0.04309, 0.64224,length.out=200)
#
for(boot_iter in 1:N){
  #
  mdata <- gen_data(pars_0,n)
  mdata_short <- c(colSums(mdata),n)
  #
  pars_mle <- mle(initial_vals,mdata_short,lower_bounds,upper_bounds)
  #
  smc_run <- SMC(mdata_short,lower_bounds,upper_bounds,keep,scalepar,cl)
  #
  accept[boot_iter] <- smc_run$accept_ratio
  #
  # critical values
  # procedure 1: full set
  xi_1 <- p1_crit_weighted(pars_mle,mcmc_draws=smc_run$Draws,mdata_short,weights=smc_run$weights,alpha_vec,cl)
  p1_full_check[boot_iter,] <- full_check(pars_mle,pars_0,mdata_short,xi=xi_1)
  #
  # store critical value
  Xi_1[boot_iter,] <- xi_1
  #
  # procedure 2: profiling for mu
  xi_2_Delta_1 <- p2_crit_Delta_1_weighted(pars_mle,smc_run$Draws,mdata_short,weights=smc_run$weights,alpha_vec,cl)
  xi_2_beta_1 <- p2_crit_beta_1_weighted(pars_mle,smc_run$Draws,mdata_short,weights=smc_run$weights,alpha_vec,cl)
  #
  # store critical value
  Xi_2_Delta_1[boot_iter,] <- xi_2_Delta_1
  Xi_2_beta_1[boot_iter,] <- xi_2_beta_1
  #
  # procedure 3: chi-square 1 upper bound
  xi_3 <- qchisq(alpha_vec,df=1)/2
  #
  subvec_Delta_1_cvg_out <- subvec_Delta_1_cvg(pars_mle,M_I_Delta_1,mdata_short,lower_bounds,upper_bounds,cl,xi=c(xi_1,xi_2_Delta_1,xi_3))
  subvec_Delta_1_length_out <- subvec_Delta_1_length(pars_mle,mdata_short,lower_bounds,upper_bounds,cl,xi=c(xi_1,xi_2_Delta_1,xi_3))
  #
  subvec_beta_1_cvg_out <- subvec_beta_1_cvg(pars_mle,M_I_beta_1,mdata_short,lower_bounds,upper_bounds,cl,xi=c(xi_1,xi_2_beta_1,xi_3))
  subvec_beta_1_length_out <- subvec_beta_1_length(pars_mle,mdata_short,lower_bounds,upper_bounds,cl,xi=c(xi_1,xi_2_beta_1,xi_3))
  #
  # procedure 1: projection for mu
  p1_proj_check[boot_iter,] <- c(subvec_Delta_1_cvg_out[(1:3)],subvec_beta_1_cvg_out[(1:3)])
  p1_lower[boot_iter,] <- c(subvec_Delta_1_length_out[1,(1:3)],subvec_beta_1_length_out[1,(1:3)])
  p1_upper[boot_iter,] <- c(subvec_Delta_1_length_out[2,(1:3)],subvec_beta_1_length_out[2,(1:3)])
  p1_length[boot_iter,] <- c(subvec_Delta_1_length_out[3,(1:3)],subvec_beta_1_length_out[3,(1:3)])
  #
  # procedure 2: profiling for mu
  p2_check[boot_iter,] <- c(subvec_Delta_1_cvg_out[(4:6)],subvec_beta_1_cvg_out[(4:6)])
  p2_lower[boot_iter,] <- c(subvec_Delta_1_length_out[1,(4:6)],subvec_beta_1_length_out[1,(4:6)])
  p2_upper[boot_iter,] <- c(subvec_Delta_1_length_out[2,(4:6)],subvec_beta_1_length_out[2,(4:6)])
  p2_length[boot_iter,] <- c(subvec_Delta_1_length_out[3,(4:6)],subvec_beta_1_length_out[3,(4:6)])
  #
  # procedure 3: chi-square 1 upper bound
  p3_check[boot_iter,] <- c(subvec_Delta_1_cvg_out[(7:9)],subvec_beta_1_cvg_out[(7:9)])
  p3_lower[boot_iter,] <- c(subvec_Delta_1_length_out[1,(7:9)],subvec_beta_1_length_out[1,(7:9)])
  p3_upper[boot_iter,] <- c(subvec_Delta_1_length_out[2,(7:9)],subvec_beta_1_length_out[2,(7:9)])
  p3_length[boot_iter,] <- c(subvec_Delta_1_length_out[3,(7:9)],subvec_beta_1_length_out[3,(7:9)])
  #
  # percentile
  pc_Delta_1_out <- percentile_CS_weighted(smc_run$Draws[,1],weights=smc_run$weights,alpha_vec,M_I_Delta_1)
  pc_beta_1_out <- percentile_CS_weighted(smc_run$Draws[,3],weights=smc_run$weights,alpha_vec,M_I_beta_1)
  pc_check[boot_iter,] <- c(pc_Delta_1_out$cvg_check,pc_beta_1_out$cvg_check)
  pc_lower[boot_iter,] <- c(pc_Delta_1_out$lower,pc_beta_1_out$lower)
  pc_upper[boot_iter,] <- c(pc_Delta_1_out$upper,pc_beta_1_out$upper)
  pc_length[boot_iter,] <- c(pc_Delta_1_out$length,pc_beta_1_out$length)
  #
  # print outputs
  if(boot_iter%%1==0){
    cat("Iteration: ",boot_iter,". Time: ",as.character(Sys.time()),". Chain acceptance: ",sum(accept)/boot_iter,".\n",sep="")
  }
  if(boot_iter%%10==0){
    print_p1_full <- colSums(p1_full_check[1:boot_iter,])/boot_iter
    cat("Coverage p1 full: (",print_p1_full[1],",",print_p1_full[2],",",print_p1_full[3],").\n",sep="")
    print_p1_proj <- colSums(p1_proj_check[1:boot_iter,])/boot_iter
    cat("Coverage p1 proj: (",print_p1_proj[1],",",print_p1_proj[2],",",print_p1_proj[3],",",print_p1_proj[4],",",print_p1_proj[5],",",print_p1_proj[6],").\n",sep="")
    print_p2 <- colSums(p2_check[1:boot_iter,])/boot_iter
    cat("Coverage p2:      (",print_p2[1],",",print_p2[2],",",print_p2[3],",",print_p2[4],",",print_p2[5],",",print_p2[6],").\n",sep="")
    print_p3 <- colSums(p3_check[1:boot_iter,])/boot_iter
    cat("Coverage p3:      (",print_p3[1],",",print_p3[2],",",print_p3[3],",",print_p3[4],",",print_p3[5],",",print_p3[6],").\n",sep="")
    print_pc <- colSums(pc_check[1:boot_iter,])/boot_iter
    cat("Coverage PC:      (",print_pc[1],",",print_pc[2],",",print_pc[3],",",print_pc[4],",",print_pc[5],",",print_pc[6],").\n",sep="")
  }
}
boot_results <- list(p1_full_check=p1_full_check,
                     p1_proj_check=p1_proj_check,
                     p2_check=p2_check,
                     p3_check=p3_check,
                     pc_check=pc_check,
                     p1_lower=p1_lower,
                     p2_lower=p2_lower,
                     p3_lower=p3_lower,
                     pc_lower=pc_lower,
                     p1_upper=p1_upper,
                     p2_upper=p2_upper,
                     p3_upper=p3_upper,
                     pc_upper=pc_upper,
                     p1_length=p1_length,
                     p2_length=p2_length,
                     p3_length=p3_length,
                     pc_length=pc_length,
                     Xi_1=Xi_1,Xi_2_Delta_1=Xi_2_Delta_1,Xi_2_beta_1=Xi_2_beta_1,Xi_3=Xi_3,n=n,N=N,pars_0=pars_0)
save(boot_results,file="results/n500/1.RData")
#
stopCluster(cl)
#
#END