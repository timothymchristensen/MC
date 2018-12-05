#
# Third file for Helpman-Melitz-Rubinstein Trade Application
# - constructs CSs for friction variables via profiling
# - set up to run as an array job 
#
rm(list=ls())
#
library(foreign)
library(parallel)
#
source("fns/trade_likelihood.R")
source("fns/likelihood_ratios.R")
source("fns/SMC.R")
#
load("trade.RData")
#
alpha_vec <- results$alpha_vec
mle_out <- results$mle_out
smc_run <- results$smc_run
data <- results$data
het <- TRUE
#
xi_1 <- p1_crit_weighted(mle_out$mle,smc_run$Draws,data,het,smc_run$weights,alpha_vec)
#
Xi_2 <- results$Xi_2
#
xi_3 <- qchisq(df=1,alpha_vec)/2
#
# construct CSs
slurm_taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(slurm_taskid)
fixed_ix <- as.numeric(slurm_taskid)
#
subvec_length_out <- subvec_length(fixed_ix,mle_out,data,het,cl,xi=c(xi_1,Xi_2[fixed_ix-1,],xi_3))
#
path <- paste0("subvec_",fixed_ix,".RData")
save(file=path,subvec_length_out)
#
# END