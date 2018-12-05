#
# First file for Helpman-Melitz-Rubinstein Trade Application
# - generates draws from posterior via SMC
# - computes t-statistic based CSs and percentile CSs
#
rm(list=ls())
#
library(foreign)
library(parallel)
#
source("fns/clean_data.R")
source("fns/trade_likelihood.R")
source("fns/SMC.R")
source("fns/two_step.R")
#
alpha_vec <- c(0.90,0.95,0.99)
#
# clean and import the data
# note: trade friction variables are stored as negatives, so the coefficients have the correct sign
data <- clean_data()$data_1986
#
# MLE 
initial_vals <- two_step(data,FALSE)
mle_out_0 <- mle(initial_vals,data,FALSE)
#
het <- TRUE
initial_vals <- two_step(data,het)
mle_out <- mle(initial_vals,data,het)
#
# parameters for SMC routine
scalepar <- list(Nphi=200,Nmh=8,Nblocks=6,lambda=2,scale_cov=1)
prior_par <- 100
keep <- 10000
#
# run SMC routine
cl <- makeCluster(getOption("cl.cores", 6),type="PSOCK")
clusterExport(cl=cl, list("loglike_trade","PKernel_trade","log_prior","log_initial","log_initial_pk","mle"))
smc_run <- SMC(data,het,prior_par,keep,scalepar,mle_out,cl)
stopCluster(cl)
#
results <- list(smc_run=smc_run,mle_out=mle_out,mle_out_0=mle_out_0,data=data,alpha_vec=alpha_vec,scaelpar=scalepar)
save(file="trade.RData",results)
#
# confidence sets assuming point identification
ml_lower_0 <- ml_upper_0 <- ml_lower <- ml_upper <- matrix(NA,9,length(alpha_vec))
zcrit <- qnorm(1-(1-alpha_vec)/2)
for(i in 2:10){
  ml_lower_0[i-1,] <- mle_out_0$mle[i] - zcrit*mle_out_0$ses[i]
  ml_lower[i-1,]   <- mle_out$mle[i] - zcrit*mle_out$ses[i]
  ml_upper_0[i-1,] <- mle_out_0$mle[i] + zcrit*mle_out_0$ses[i]
  ml_upper[i-1,]   <- mle_out$mle[i] + zcrit*mle_out$ses[i]
}
#
# percentile CS
pc_lower <- pc_upper <- matrix(NA,9,length(alpha_vec))
for(i in 2:10){
  pc_lower[i-1,] <- wquantile(smc_run$Draws[,i],smc_run$weights,(1-alpha_vec)/2)
  pc_upper[i-1,] <- wquantile(smc_run$Draws[,i],smc_run$weights,1-(1-alpha_vec)/2)
}
#
results <- c(results,list(ml_lower_0=ml_lower_0,ml_upper_0=ml_upper_0,ml_lower=ml_lower,ml_upper=ml_upper,pc_lower=pc_lower,pc_upper=pc_upper))
save(file="trade.RData",results)
#
# END