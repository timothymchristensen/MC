#
# Entry game application
#
rm(list=ls())
#
source("fns/get_data.R") 
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
alpha_vec <- c(0.9,0.95,0.99)
#
# same selection probabilties across regressors?
ssame <- TRUE
#
scalepar <- list(Nphi=200,Nmh=4-ssame*2,Nblocks=2,lambda=2,scale_cov=1)
keep <- 10000
#
if(ssame){
  lower_bounds <- c(-2,-2,rep(-1,3),rep(-1,3),0,0)
  upper_bounds <- c(0,0,rep(2,3),rep(2,3),1,1)
  path <- "game0.RData"
}else{
  lower_bounds <- c(-2,-2,rep(-1,3),rep(-1,3),0,rep(0,8))
  upper_bounds <- c(0,0,rep(2,3),rep(2,3),1,rep(1,8))
  path <- "game1.RData"
}
npar <- length(lower_bounds)
#
# pre-allocate for coverage checking and CS length
p1_lower <- p2_lower <- p3_lower <- pc_lower <- p1_upper <- p2_upper <- p3_upper <- pc_upper <- p1_length <- p2_length <- p3_length <- pc_length <- matrix(0,nrow=npar,ncol=length(alpha_vec))
Xi_2 <- matrix(NA,nrow=npar,ncol=length(alpha_vec))
#
data <- get_data()
mdata_short <- data$mdata_short
reg_short <- data$reg_short
#
# set up parallelization 
cl1  <- makeCluster(getOption("cl.cores", 1),type="FORK")
cl20 <- makeCluster(getOption("cl.cores", 20),type="FORK")
#
ptm <- proc.time()
smc_run <- SMC(mdata_short,lower_bounds,upper_bounds,ssame,keep,scalepar,cl1)
cat("smc time",as.numeric(proc.time()-ptm)[3],"\n")
#
pars_mle <- mle(colMeans(smc_run$Draws),mdata_short,lower_bounds,upper_bounds,ssame)
#
# critical values
# procedure 1: full set
xi_1 <- p1_crit_weighted(pars_mle,mcmc_draws=smc_run$Draws,mdata_short,ssame,weights=smc_run$weights,alpha_vec,cl20)
#
# procedure 2: profiling
for(fixed_ix in 1:npar){
  ptm <- proc.time()
  Xi_2[fixed_ix,] <- p2_crit_weighted(fixed_ix,pars_mle,smc_run$Draws,mdata_short,lower_bounds,upper_bounds,ssame,weights=smc_run$weights,alpha_vec,cl20)
  cat("p2 time for index",fixed_ix,": ",as.numeric(proc.time()-ptm)[3],"\n")
}
#
game_results <- list(p1_lower=p1_lower,
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
                     xi_1=xi_1,
                     Xi_2=Xi_2,
                     smc_run=smc_run)
save(game_results,file=path)
#
# procedure 3: chi-square 1 upper bound
xi_3 <- qchisq(alpha_vec,df=1)/2
#
# CSs
for(fixed_ix in 1:npar){
  ptm <- proc.time()
  subvec_length_out <- subvec_length(fixed_ix,pars_mle,mdata_short,lower_bounds,upper_bounds,ssame,cl20,xi=c(xi_1,Xi_2[fixed_ix,],xi_3))
  #
  p1_lower[fixed_ix,] <- subvec_length_out[1,(1:3)]
  p1_upper[fixed_ix,] <- subvec_length_out[2,(1:3)]
  p1_length[fixed_ix,] <- subvec_length_out[3,(1:3)]
  #
  p2_lower[fixed_ix,] <- subvec_length_out[1,(4:6)]
  p2_upper[fixed_ix,] <- subvec_length_out[2,(4:6)]
  p2_length[fixed_ix,] <- subvec_length_out[3,(4:6)]
  #
  p3_lower[fixed_ix,] <- subvec_length_out[1,(7:9)]
  p3_upper[fixed_ix,] <- subvec_length_out[2,(7:9)]
  p3_length[fixed_ix,] <- subvec_length_out[3,(7:9)]
  cat("CS time for index",fixed_ix,": ",as.numeric(proc.time()-ptm)[3],"\n")
}
#
# percentile
for(fixed_ix in 1:npar){
  pc_out <- percentile_CS_weighted(smc_run$Draws[,fixed_ix],weights=smc_run$weights,alpha_vec)
  pc_lower[fixed_ix,] <- pc_out$lower
  pc_upper[fixed_ix,] <- pc_out$upper
  pc_length[fixed_ix,] <- pc_out$length
}
#
game_results <- list(p1_lower=p1_lower,
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
                     xi_1=xi_1,
                     Xi_2=Xi_2,
                     smc_run=smc_run)
save(game_results,file=path)
#
stopCluster(cl1)
stopCluster(cl20)
#
#END