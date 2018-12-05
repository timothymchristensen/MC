#
# populates tables for Example 1
#
cvg_full <- array(NA,c(4,9,6))
lb_full <- ub_full <- array(NA,c(4,9,5))
#
cvg_read <- function(path){
  #
  load(path)
  #
  cvg_temp <- array(NA,c(1,3,6))
  lb_temp <- ub_temp <- array(NA,c(1,3,5))
  #
  cvg_temp[1,(1:3),1] <- colMeans(boot_results$p1_full_check)
  cvg_temp[1,(1:3),2] <- colMeans(boot_results$p2_check)
  cvg_temp[1,(1:3),3] <- colMeans(boot_results$p3_check)
  cvg_temp[1,(1:3),4] <- colMeans(boot_results$p1_proj_check)
  cvg_temp[1,(1:3),5] <- colMeans(boot_results$pc_check)
  cvg_temp[1,(1:3),6] <- colMeans(boot_results$mi_check)
  #
  lb_temp[1,(1:3),1] <- colMeans(boot_results$p2_lower)
  lb_temp[1,(1:3),2] <- colMeans(boot_results$p3_lower)
  lb_temp[1,(1:3),3] <- colMeans(boot_results$p1_lower)
  lb_temp[1,(1:3),4] <- colMeans(boot_results$pc_lower)
  lb_temp[1,(1:3),5] <- colMeans(boot_results$mi_lower)
  #
  ub_temp[1,(1:3),1] <- colMeans(boot_results$p2_upper)
  ub_temp[1,(1:3),2] <- colMeans(boot_results$p3_upper)
  ub_temp[1,(1:3),3] <- colMeans(boot_results$p1_upper)
  ub_temp[1,(1:3),4] <- colMeans(boot_results$pc_upper)
  ub_temp[1,(1:3),5] <- colMeans(boot_results$mi_upper)
  #
  return(list(cvg_temp=cvg_temp,lb_temp=lb_temp,ub_temp=ub_temp))
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
    path <- paste0("./",pn,"/",3-j,".RData")
    #
    cvg_out <- cvg_read(path)
    #
    cvg_full[i,((3*(j-1)+1):(3*j)),] <- cvg_out$cvg_temp
    lb_full[i,((3*(j-1)+1):(3*j)),] <- cvg_out$lb_temp
    ub_full[i,((3*(j-1)+1):(3*j)),] <- cvg_out$ub_temp
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
numformat <- function(val) { sub("^(-?)0.", "\\1.", format(round(val,digits=2),nsmall=2)) }
#
numformat3 <- function(val) { sub("^(-?)0.", "\\1.", format(round(val,digits=3),nsmall=3)) }
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
cat("	& & \\multicolumn{16}{c}{$\\widehat{\\Theta}_{\\alpha}$ (Procedure 1) for $\\Theta_I$} & \\\\","\n")
cat("100  ",stringtemp(cvg_full[1,,1]),"\\\\ \n")
cat("250  ",stringtemp(cvg_full[2,,1]),"\\\\ \n")
cat("500  ",stringtemp(cvg_full[3,,1]),"\\\\ \n")
cat("1000 ",stringtemp(cvg_full[4,,1]),"\\\\ \n")
cat("	& & \\multicolumn{16}{c}{$\\widehat{M}_{\\alpha}$ (Procedure 2) for $M_I$} & \\\\","\n")
cat("100  ",inttemp(cvg_full[1,,2],lb_full[1,,1],ub_full[1,,1]),"\\\\ \n")
cat("250  ",inttemp(cvg_full[2,,2],lb_full[2,,1],ub_full[2,,1]),"\\\\ \n")
cat("500  ",inttemp(cvg_full[3,,2],lb_full[3,,1],ub_full[3,,1]),"\\\\ \n")
cat("1000 ",inttemp(cvg_full[4,,2],lb_full[4,,1],ub_full[4,,1]),"\\\\ \n")
cat("	& & \\multicolumn{16}{c}{$\\widehat{M}^{\\chi}_{\\alpha}$ (Procedure 3) for $M_I$} & \\\\","\n")
cat("100  ",inttemp(cvg_full[1,,3],lb_full[1,,2],ub_full[1,,2]),"\\\\ \n")
cat("250  ",inttemp(cvg_full[2,,3],lb_full[2,,2],ub_full[2,,2]),"\\\\ \n")
cat("500  ",inttemp(cvg_full[3,,3],lb_full[3,,2],ub_full[3,,2]),"\\\\ \n")
cat("1000 ",inttemp(cvg_full[4,,3],lb_full[4,,2],ub_full[4,,2]),"\\\\ \n")
cat("	& & \\multicolumn{16}{c}{$\\widehat{M}^{proj}_{\\alpha}$ (Projection) for $M_I$} & \\\\","\n")
cat("100  ",inttemp(cvg_full[1,,4],lb_full[1,,3],ub_full[1,,3]),"\\\\ \n")
cat("250  ",inttemp(cvg_full[2,,4],lb_full[2,,3],ub_full[2,,3]),"\\\\ \n")
cat("500  ",inttemp(cvg_full[3,,4],lb_full[3,,3],ub_full[3,,3]),"\\\\ \n")
cat("1000 ",inttemp(cvg_full[4,,4],lb_full[4,,3],ub_full[4,,3]),"\\\\ \n")
cat("	& & \\multicolumn{16}{c}{$\\widehat{M}^{perc}_{\\alpha}$ (Percentiles) for $M_I$} & \\\\","\n")
cat("100  ",inttemp(cvg_full[1,,5],lb_full[1,,4],ub_full[1,,4]),"\\\\ \n")
cat("250  ",inttemp(cvg_full[2,,5],lb_full[2,,4],ub_full[2,,4]),"\\\\ \n")
cat("500  ",inttemp(cvg_full[3,,5],lb_full[3,,4],ub_full[3,,4]),"\\\\ \n")
cat("1000 ",inttemp(cvg_full[4,,5],lb_full[4,,4],ub_full[4,,4]),"\\\\ \n")
cat("	& & \\multicolumn{16}{c}{GMS CSs for $\\mu$ via moment inequalities} & \\\\","\n")
cat("100  ",inttemp(cvg_full[1,,6],lb_full[1,,5],ub_full[1,,5]),"\\\\ \n")
cat("250  ",inttemp(cvg_full[2,,6],lb_full[2,,5],ub_full[2,,5]),"\\\\ \n")
cat("500  ",inttemp(cvg_full[3,,6],lb_full[3,,5],ub_full[3,,5]),"\\\\ \n")
cat("1000 ",inttemp(cvg_full[4,,6],lb_full[4,,5],ub_full[4,,5]),"\\\\ \n")
cat("\n")
sink()


