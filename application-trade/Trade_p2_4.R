#
# Helpman-Melitz-Rubinstein Trade Application
# - generates procedure 2 critical values
#
# Expect command line args at the end
slurm_taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(slurm_taskid)
# coerce the value to an integer
num_ix <- as.numeric(slurm_taskid)
print(num_ix)
#
par_ix <- 4
#
source("../fns/trade_likelihood.R")
source("../fns/likelihood_ratios.R")
source("../fns/ID_set.R")
#
load("../trade.RData")
smc_run <- results$smc_run
mle_out <- results$mle_out
data <- results$data
het  <- TRUE
#
mle_LL <- -loglike_trade(mle_out$mle,data,het)
#
draws_ix <- seq((500*(num_ix-1)+1),500*num_ix)
head(draws_ix)
fixed_ix <- par_ix
#
pqlr <- numeric(500)
for(i in 1:500){
  pqlr[i] <- mle_LL - p2_crit_fn(smc_run$Draws[draws_ix[i],],fixed_ix,mle_out,data,het)
  cat("iteration: ",i," parameter: ",par_ix," block: ",num_ix,"\n")
}
#
path <- paste0("p2_",par_ix,"_",num_ix,".RData")
save(file=path,pqlr)
#
# END