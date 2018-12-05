#
# populates tables for Example 2
#
cvg_full <- matrix(NA,4,3)
cvg_delta <- cvg_beta <- lb_delta <- ub_delta <- lb_beta <- ub_beta <- array(NA,c(4,3,4))
#
cvg_read <- function(path){
  #
  load(path)
  #
  cvg_beta_temp <- lb_beta_temp <- ub_beta_temp <- cvg_delta_temp <- lb_delta_temp <- ub_delta_temp <- array(NA,c(1,3,4))
  #
  cvg_temp <- colMeans(boot_results$p1_full_check)
  #
  cvg_delta_temp[1,,1] <- colMeans(boot_results$p2_check[,1:3])
  cvg_delta_temp[1,,2] <- colMeans(boot_results$p3_check[,1:3])
  cvg_delta_temp[1,,3] <- colMeans(boot_results$p1_proj_check[,1:3])
  cvg_delta_temp[1,,4] <- colMeans(boot_results$pc_check[,1:3])
  #
  cvg_beta_temp[1,,1] <- colMeans(boot_results$p2_check[,4:6])
  cvg_beta_temp[1,,2] <- colMeans(boot_results$p3_check[,4:6])
  cvg_beta_temp[1,,3] <- colMeans(boot_results$p1_proj_check[,4:6])
  cvg_beta_temp[1,,4] <- colMeans(boot_results$pc_check[,4:6])
  #
  lb_delta_temp[1,,1] <- colMeans(boot_results$p2_lower[,1:3])
  lb_delta_temp[1,,2] <- colMeans(boot_results$p3_lower[,1:3])
  lb_delta_temp[1,,3] <- colMeans(boot_results$p1_lower[,1:3])
  lb_delta_temp[1,,4] <- colMeans(boot_results$pc_lower[,1:3])
  #
  lb_beta_temp[1,,1] <- colMeans(boot_results$p2_lower[,4:6])
  lb_beta_temp[1,,2] <- colMeans(boot_results$p3_lower[,4:6])
  lb_beta_temp[1,,3] <- colMeans(boot_results$p1_lower[,4:6])
  lb_beta_temp[1,,4] <- colMeans(boot_results$pc_lower[,4:6])
  #
  ub_delta_temp[1,,1] <- colMeans(boot_results$p2_upper[,1:3])
  ub_delta_temp[1,,2] <- colMeans(boot_results$p3_upper[,1:3])
  ub_delta_temp[1,,3] <- colMeans(boot_results$p1_upper[,1:3])
  ub_delta_temp[1,,4] <- colMeans(boot_results$pc_upper[,1:3])
  #
  ub_beta_temp[1,,1] <- colMeans(boot_results$p2_upper[,4:6])
  ub_beta_temp[1,,2] <- colMeans(boot_results$p3_upper[,4:6])
  ub_beta_temp[1,,3] <- colMeans(boot_results$p1_upper[,4:6])
  ub_beta_temp[1,,4] <- colMeans(boot_results$pc_upper[,4:6])
  #
  return(list(cvg_temp=cvg_temp,cvg_delta_temp=cvg_delta_temp,cvg_beta_temp=cvg_beta_temp,lb_delta_temp=lb_delta_temp,ub_delta_temp=ub_delta_temp,lb_beta_temp=lb_beta_temp,ub_beta_temp=ub_beta_temp))
}
#
for(i in 1:4){
  #
  if(i == 1){
    pn <- "n100"
  }
  if(i == 2){
    pn <- "n250"
  }
  if(i == 3){
    pn <- "n500"
  }
  if(i == 4){
    pn <- "n1000"
  }
  #
  for(j in 1:3){
    #
    path <- paste0("./",pn,"/1.RData")
    #
    cvg_out <- cvg_read(path)
    #
    cvg_full[i,] <- cvg_out$cvg_temp
    cvg_delta[i,,] <- cvg_out$cvg_delta_temp
    cvg_beta[i,,] <- cvg_out$cvg_beta_temp
    lb_delta[i,,] <- cvg_out$lb_delta_temp
    lb_beta[i,,] <- cvg_out$lb_beta_temp
    ub_delta[i,,] <- cvg_out$ub_delta_temp
    ub_beta[i,,] <- cvg_out$ub_beta_temp
  }
}
#
stringtemp <- function(rvec){
  str <- paste()
  for(i in 1:length(rvec)){
    str <- paste(str," & ",numformat3(rvec[i])," & --- ")
  }
  return(str)
}
#
numformat <- function(val) { format(round(val,digits=2),nsmall=2) }
#
numformat3 <- function(val) { format(round(val,digits=3),nsmall=3) }
#
inttemp <- function(rvec,lvec,uvec){
  str <- paste()
  for(i in 1:length(lvec)){
    str <- paste(str," & ",numformat3(rvec[i])," & [$",numformat(lvec[i]),",\\!",numformat(uvec[i]),"$]")
  }
  return(str)
}
#
sink('tab_int.txt',append=FALSE)
cat("\n")
cat(" & & \\multicolumn{4}{c}{CSs for the identified set $\\Theta_I$} & \\\\ \n") 
cat(" & & \\multicolumn{4}{c}{$\\widehat{\\Theta}_{\\alpha}$ (Procedure 1)} & \\\\ \n")
cat("100  ",stringtemp(cvg_full[1,]),"\\\\ \n")
cat("250  ",stringtemp(cvg_full[2,]),"\\\\ \n")
cat("500  ",stringtemp(cvg_full[3,]),"\\\\ \n")
cat("1000 ",stringtemp(cvg_full[4,]),"\\\\ \\hline \n")
cat("& & \\multicolumn{4}{c}{CSs for the identified set for $\\Delta_1$} & \\\\ \n ")
cat("& & \\multicolumn{4}{c}{$\\widehat{M}_{\\alpha}$ (Procedure 2)} & \\\\ \n")
cat("100  ",inttemp(cvg_delta[1,,1],lb_delta[1,,1],ub_delta[1,,1]),"\\\\ \n")
cat("250  ",inttemp(cvg_delta[2,,1],lb_delta[2,,1],ub_delta[2,,1]),"\\\\ \n")
cat("500  ",inttemp(cvg_delta[3,,1],lb_delta[3,,1],ub_delta[3,,1]),"\\\\ \n")
cat("1000 ",inttemp(cvg_delta[4,,1],lb_delta[4,,1],ub_delta[4,,1]),"\\\\ \n")
cat("& & \\multicolumn{4}{c}{$\\widehat{M}_{\\alpha}^\\chi$ (Procedure 3)} & \\\\ \n") 
cat("100  ",inttemp(cvg_delta[1,,2],lb_delta[1,,2],ub_delta[1,,2]),"\\\\ \n")
cat("250  ",inttemp(cvg_delta[2,,2],lb_delta[2,,2],ub_delta[2,,2]),"\\\\ \n")
cat("500  ",inttemp(cvg_delta[3,,2],lb_delta[3,,2],ub_delta[3,,2]),"\\\\ \n")
cat("1000 ",inttemp(cvg_delta[4,,2],lb_delta[4,,2],ub_delta[4,,2]),"\\\\ \n")
cat("& & \\multicolumn{4}{c}{$\\widehat{M}_{\\alpha}^{proj}$ (Projection)} & \\\\ \n") 
cat("100  ",inttemp(cvg_delta[1,,3],lb_delta[1,,3],ub_delta[1,,3]),"\\\\ \n")
cat("250  ",inttemp(cvg_delta[2,,3],lb_delta[2,,3],ub_delta[2,,3]),"\\\\ \n")
cat("500  ",inttemp(cvg_delta[3,,3],lb_delta[3,,3],ub_delta[3,,3]),"\\\\ \n")
cat("1000 ",inttemp(cvg_delta[4,,3],lb_delta[4,,3],ub_delta[4,,3]),"\\\\ \n")
cat("& & \\multicolumn{4}{c}{$\\widehat{M}_{\\alpha}^{perc}$ (Percentiles)} & \\\\ \n")
cat("100  ",inttemp(cvg_delta[1,,4],lb_delta[1,,4],ub_delta[1,,4]),"\\\\ \n")
cat("250  ",inttemp(cvg_delta[2,,4],lb_delta[2,,4],ub_delta[2,,4]),"\\\\ \n")
cat("500  ",inttemp(cvg_delta[3,,4],lb_delta[3,,4],ub_delta[3,,4]),"\\\\ \n")
cat("1000 ",inttemp(cvg_delta[4,,4],lb_delta[4,,4],ub_delta[4,,4]),"\\\\ \\hline \n")
cat("& & \\multicolumn{4}{c}{CSs for the identified set for $\\beta_1$} & \\\\ \n ")
cat("& & \\multicolumn{4}{c}{$\\widehat{M}_{\\alpha}$ (Procedure 2)} & \\\\ \n")
cat("100  ",inttemp(cvg_beta[1,,1],lb_beta[1,,1],ub_beta[1,,1]),"\\\\ \n")
cat("250  ",inttemp(cvg_beta[2,,1],lb_beta[2,,1],ub_beta[2,,1]),"\\\\ \n")
cat("500  ",inttemp(cvg_beta[3,,1],lb_beta[3,,1],ub_beta[3,,1]),"\\\\ \n")
cat("1000 ",inttemp(cvg_beta[4,,1],lb_beta[4,,1],ub_beta[4,,1]),"\\\\ \n")
cat("& & \\multicolumn{4}{c}{$\\widehat{M}_{\\alpha}^\\chi$ (Procedure 3)} & \\\\ \n") 
cat("100  ",inttemp(cvg_beta[1,,2],lb_beta[1,,2],ub_beta[1,,2]),"\\\\ \n")
cat("250  ",inttemp(cvg_beta[2,,2],lb_beta[2,,2],ub_beta[2,,2]),"\\\\ \n")
cat("500  ",inttemp(cvg_beta[3,,2],lb_beta[3,,2],ub_beta[3,,2]),"\\\\ \n")
cat("1000 ",inttemp(cvg_beta[4,,2],lb_beta[4,,2],ub_beta[4,,2]),"\\\\ \n")
cat("& & \\multicolumn{4}{c}{$\\widehat{M}_{\\alpha}^{proj}$ (Projection)} & \\\\ \n") 
cat("100  ",inttemp(cvg_beta[1,,3],lb_beta[1,,3],ub_beta[1,,3]),"\\\\ \n")
cat("250  ",inttemp(cvg_beta[2,,3],lb_beta[2,,3],ub_beta[2,,3]),"\\\\ \n")
cat("500  ",inttemp(cvg_beta[3,,3],lb_beta[3,,3],ub_beta[3,,3]),"\\\\ \n")
cat("1000 ",inttemp(cvg_beta[4,,3],lb_beta[4,,3],ub_beta[4,,3]),"\\\\ \n")
cat("& & \\multicolumn{4}{c}{$\\widehat{M}_{\\alpha}^{perc}$ (Percentiles)} & \\\\ \n")
cat("100  ",inttemp(cvg_beta[1,,4],lb_beta[1,,4],ub_beta[1,,4]),"\\\\ \n")
cat("250  ",inttemp(cvg_beta[2,,4],lb_beta[2,,4],ub_beta[2,,4]),"\\\\ \n")
cat("500  ",inttemp(cvg_beta[3,,4],lb_beta[3,,4],ub_beta[3,,4]),"\\\\ \n")
cat("1000 ",inttemp(cvg_beta[4,,4],lb_beta[4,,4],ub_beta[4,,4]),"\\\\ \\hline \n")
cat("\n")
sink()



