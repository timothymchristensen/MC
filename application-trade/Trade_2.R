#
# Second file for Helpman-Melitz-Rubinstein Trade Application
# - imports procedure 2 critical values and saves
# - run Trade_p2_k.R before this for k = 2,...,10 to generate PQ(M(theta^b))
#
rm(list=ls())
#
load("trade.RData")
#
alpha_vec <- results$alpha_vec
#
Xi_2 <- matrix(NA,9,length(alpha_vec))
#
for(i in 2:10){
  pq <- numeric(0)
  for(j in 1:50){
    #
    path <- paste0("p2_",i,"_",j,".RData")
    if(file.exists(path)){
      load(path)
      pq <- c(pq,pqlr)
    }
  }
  Xi_2[i-1,] <- quantile(pq,alpha_vec)
}
#
results <- c(results,list(Xi_2=Xi_2))
save(file="trade.RData",results)
#
# END