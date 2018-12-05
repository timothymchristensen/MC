#
# Example 1: Missing Data via SMC algorithm and CU-GMM criterion function
#
rm(list=ls())
#
n <- 500  # sample size
N <- 5000  # number of replications
#
source("fns/gen_data.R")
source("fns_gmm/moments.R")
source("fns_gmm/gmm_objfn.R")
source("fns_gmm/gmm_argmin.R")
source("fns_gmm/qlrs.R")
source("fns_gmm/SMC.R")
source("fns/ex1_gms.R")
require("parallel")
#
set.seed(1234567)
#
trans_pars <- TRUE    # whether to transform the parameters to have unbounded support
#
mu    <- 0.50    # Pr(Y_i = 1)
beta  <- 0.50    # Pr(Y_i = 1 | D_i = 0)
gamma <- 1-2/sqrt(n)    # Pr(D_i = 1)
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
scalepar <- list(Nphi=200,Nmh=1,lambda=2,scale_cov=1)
keep <- 10000   # number of draws
#
alpha_vec <- c(0.9,0.95,0.99)  # significance levels
n_alpha <- length(alpha_vec)
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
clusterExport(cl=cl, list("enforce_restrictions","gmm_objfn","PKernel_temp","gmm_mu_subvec","gmm_objfn","gmm_objfn_profile"))
#
for(boot_iter in 1:N){
  #
  mdata <- gen_data(pars_0,n)
  #
  # pre-compute some statistics for GMM criterion
  mdata_short <- data_trans(mdata)$mdata_short
  W_mat <- data_trans(mdata)$W_mat
  #
  initial_vals <- c(0.45,0.45,mdata_short[2]) 
  pars_hat <- gmm_argmin(initial_vals,mdata_short,W_mat)
  #
  smc_run <- SMC(mdata_short,W_mat,trans_pars,keep,scalepar,cl)
  #
  accept[boot_iter] <- smc_run$accept_ratio
  #
  # critical values
  # procedure 1: full set
  xi_1 <- p1_crit_weighted(pars_hat,mcmc_draws=smc_run$Draws,mdata_short,W_mat,weights=smc_run$weights,alpha_vec)
  Xi_1[boot_iter,] <- xi_1
  #
  p1_full_check[boot_iter,] <- full_check(pars_hat,pars_0,mdata_short,W_mat,xi=xi_1)
  #
  # procedure 2: profiling for mu
  xi_2 <- p2_crit_par_weighted(pars_hat,smc_run$Draws,mdata_short,W_mat,weights=smc_run$weights,alpha_vec,cl)
  Xi_2[boot_iter,] <- xi_2
  #
  # procedure 3: chi-square 1 upper bound
  xi_3 <- qchisq(alpha_vec,df=1)
  Xi_3[boot_iter,] <- xi_3
  #
  subvec_check <- subvec_mu_cvg(pars_hat,M_I,mdata_short,W_mat,xi=c(xi_1,xi_2,xi_3))
  subvec_out <- subvec_mu_length(pars_hat,mdata_short,W_mat,xi=c(xi_1,xi_2,xi_3))
  #
  p1_proj_check[boot_iter,] <- subvec_check[1:n_alpha]
  p1_lower[boot_iter,] <- subvec_out[1,1:n_alpha]
  p1_upper[boot_iter,] <- subvec_out[2,1:n_alpha]
  p1_length[boot_iter,] <- subvec_out[3,1:n_alpha]
  #
  p2_check[boot_iter,] <- subvec_check[(n_alpha+1):(2*n_alpha)]
  p2_lower[boot_iter,] <- subvec_out[1,(n_alpha+1):(2*n_alpha)]
  p2_upper[boot_iter,] <- subvec_out[2,(n_alpha+1):(2*n_alpha)]
  p2_length[boot_iter,] <- subvec_out[3,(n_alpha+1):(2*n_alpha)]
  #
  p3_check[boot_iter,] <- subvec_check[(2*n_alpha+1):(3*n_alpha)]
  p3_lower[boot_iter,] <- subvec_out[1,(2*n_alpha+1):(3*n_alpha)]
  p3_upper[boot_iter,] <- subvec_out[2,(2*n_alpha+1):(3*n_alpha)]
  p3_length[boot_iter,] <- subvec_out[3,(2*n_alpha+1):(3*n_alpha)]
  #
  # percentile
  pc_out <- percentile_CS_weighted(smc_run$Draws[,1],weights=smc_run$weights,alpha_vec,M_I)
  pc_check[boot_iter,] <- pc_out$cvg_check
  pc_lower[boot_iter,] <- pc_out$lower
  pc_upper[boot_iter,] <- pc_out$upper
  pc_length[boot_iter,] <- pc_out$length
  #
  # comparison with moment inequalities
  GMS_out <- GMS_mu(mdata[,c(2,1)],alpha_vec)
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
    print_pc <- colSums(pc_check[1:boot_iter,])/boot_iter
    cat("Coverage pc:      (",print_pc[1],",",print_pc[2],",",print_pc[3],").\n",sep="")
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
save(boot_results,file="results/n500/2.RData")
#
stopCluster(cl)
#
#END
