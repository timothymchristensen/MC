#
# Final file for Helpman-Melitz-Rubinstein Trade Application
# - prints CSs 
# - saves all results
#
rm(list=ls())
#
load("trade.RData")
#
coef_of_interest <- c("Distance","Border","Island","Landlock","Legal","Language","Colonial","Currency","FTA")
alpha_ix <- 2 # for 95% CS
#
n_alpha <- length(results$alpha_vec)
p2_lower <- p2_upper <- p3_lower <- p3_upper <- matrix(NA,9,n_alpha)
#
for(i in 2:10){
  #
  mle0 <- paste0(round(results$mle_out_0$mle[i],3)," & ")
  mle1 <- paste0(round(results$mle_out$mle[i],3)," & ")
  ml0 <- paste0(" [",round(results$ml_lower_0[i-1,alpha_ix],3),",",round(results$ml_upper_0[i-1,alpha_ix],3),"] & ")
  ml1 <- paste0(" [",round(results$ml_lower[i-1,alpha_ix],3),",",round(results$ml_upper[i-1,alpha_ix],3),"] & ")
  pc1 <- paste0(" [",round(results$pc_lower[i-1,alpha_ix],3),",",round(results$pc_upper[i-1,alpha_ix],3),"]  ")
  #
  path <- paste0("subvec_",i,".RData")
  load(path)
  #
  p2_lower[i-1,] <- subvec_length_out$CS_length[1,(n_alpha+1):(2*n_alpha)]
  p2_upper[i-1,] <- subvec_length_out$CS_length[2,(n_alpha+1):(2*n_alpha)]
  #
  p3_lower[i-1,] <- subvec_length_out$CS_length[1,(2*n_alpha+1):(3*n_alpha)]
  p3_upper[i-1,] <- subvec_length_out$CS_length[2,(2*n_alpha+1):(3*n_alpha)]
  #
  p2 <- paste0(" [",round(subvec_length_out$CS_length[1,3+alpha_ix],3),",",round(subvec_length_out$CS_length[2,3+alpha_ix],3),"] & ")
  p3 <- paste0(" [",round(subvec_length_out$CS_length[1,6+alpha_ix],3),",",round(subvec_length_out$CS_length[2,6+alpha_ix],3),"] & ")
  # 
  cat(coef_of_interest[i-1]," & ",mle0,ml0,mle1,ml1,p2,p3,pc1,"\\\\ \n")
}
#
results <- c(results,list(p2_lower=p2_lower,p2_upper=p2_upper,p3_lower=p3_lower,p3_upper=p3_upper))
save(file="trade.RData",results)
#
# END