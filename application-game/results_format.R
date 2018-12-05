#
# Generate plots and table for game empirical application
#
rm(list=ls())
#
require("ggplot2")
#
ix <- 2 # 1 = 90%, 2 = 95%, 3 = 99%
#
# first make the plots
for(ssame in c(TRUE,FALSE)){
  #
  if(ssame){
    load("game0.RData")
    labels <- c(expression(Delta[OA]),
                expression(Delta[LC]),
                expression(beta[OA]^{0}),
                expression(beta[OA]^{MS}),
                expression(beta[OA]^{MP}),
                expression(beta[LC]^{0}),
                expression(beta[LC]^{MS}),
                expression(beta[LC]^{MP}),
                expression(rho),
                expression(s))
    l0 <- cbind(game_results$p2_lower[,ix],
                game_results$p3_lower[,ix],
                game_results$p1_lower[,ix],
                game_results$pc_lower[,ix])
    u0 <- cbind(game_results$p2_upper[,ix],
                game_results$p3_upper[,ix],
                game_results$p1_upper[,ix],
                game_results$pc_upper[,ix])
  }else{
    load("game1.RData")
    labels <- c(expression(Delta[OA]),
                expression(Delta[LC]),
                expression(beta[OA]^{0}),
                expression(beta[OA]^{MS}),
                expression(beta[OA]^{MP}),
                expression(beta[LC]^{0}),
                expression(beta[LC]^{MS}),
                expression(beta[LC]^{MP}),
                expression(rho),
                expression(s['000']),
                expression(s['001']),
                expression(s['010']),
                expression(s['100']),
                expression(s['011']),
                expression(s['101']),
                expression(s['110']),
                expression(s['111']))
    l1 <- cbind(game_results$p2_lower[,ix],
                game_results$p3_lower[,ix],
                game_results$p1_lower[,ix],
                game_results$pc_lower[,ix])
    u1 <- cbind(game_results$p2_upper[,ix],
                game_results$p3_upper[,ix],
                game_results$p1_upper[,ix],
                game_results$pc_upper[,ix])
  }
  #
  smc_run <- game_results$smc_run
  #
  npar <- ncol(smc_run$Draws)
  for(i in 1:npar){
    dfmu <- data.frame(smc_run$Draws[,i])
    colnames(dfmu) <- c("mu")
    saveG <- TRUE
    if(class(dev.list()) != "NULL"){dev.off()}
    path <- paste("plots/game","_",1-ssame*1,"_",i,".pdf",sep="")
    if(saveG==TRUE){cairo_pdf(file=path,height=9,width=13)}
    g <- ggplot(dfmu, aes(x = mu))+
      geom_histogram(bins=40,
                     fill=I("black"))+
      ylab("")+
      xlab(labels[i])+
      theme_minimal()+
      theme(legend.position = "right",
            legend.title=element_blank(),
            axis.text.x=element_text(size=30),
            axis.title.x=element_text(size=40),
            axis.text.y=element_text(size=30))
    print(g)
    if(saveG==TRUE){dev.off()}
  }
}
#
# populate table
#
numformat3 <- function(val) {  format(round(val,digits=3),nsmall=3) }
#
inttemp <- function(lvec,uvec){
  str <- paste()
  for(i in 1:length(lvec)){
    str <- paste(str," & [$",numformat3(lvec[i]),",\\!",numformat3(uvec[i]),"$]")
  }
  return(str)
}
#
sink('tab_int.txt',append=FALSE)
cat("\n")
cat(" &  Procedure 2 & Procedure 3 & Projection & Percentile &  Procedure 2 & Procedure 3 & Projection & Percentile \\\\ \n")
cat("$\\Delta_{OA}$     ",inttemp(cbind(l1[1,],l0[1,]),cbind(u1[1,],u0[1,])),"\\\\ \n")
cat("$\\Delta_{LC}$     ",inttemp(cbind(l1[2,],l0[2,]),cbind(u1[2,],u0[2,])),"\\\\ \n")
cat("$\\beta_{OA}^0$    ",inttemp(cbind(l1[3,],l0[3,]),cbind(u1[3,],u0[3,])),"\\\\ \n")
cat("$\\beta_{OA}^{MS}$ ",inttemp(cbind(l1[4,],l0[4,]),cbind(u1[4,],u0[4,])),"\\\\ \n")
cat("$\\beta_{OA}^{MP}$ ",inttemp(cbind(l1[5,],l0[5,]),cbind(u1[5,],u0[5,])),"\\\\ \n")
cat("$\\beta_{LC}^0$    ",inttemp(cbind(l1[6,],l0[6,]),cbind(u1[6,],u0[6,])),"\\\\ \n")
cat("$\\beta_{LC}^{MS}$ ",inttemp(cbind(l1[7,],l0[7,]),cbind(u1[7,],u0[7,])),"\\\\ \n")
cat("$\\beta_{LC}^{MP}$ ",inttemp(cbind(l1[8,],l0[8,]),cbind(u1[8,],u0[8,])),"\\\\ \n")
cat("$\\rho$            ",inttemp(cbind(l1[9,],l0[9,]),cbind(u1[9,],u0[9,])),"\\\\ \n")
cat("$s$                "," & --- & --- & --- & --- ",inttemp(l0[10,],u0[10,]),"\\\\ \n")
cat("$s_{000}$          ",inttemp(l1[10,],u1[10,])," & --- & --- & --- & --- ","\\\\ \n")
cat("$s_{001}$          ",inttemp(l1[11,],u1[11,])," & --- & --- & --- & --- ","\\\\ \n")
cat("$s_{010}$          ",inttemp(l1[12,],u1[12,])," & --- & --- & --- & --- ","\\\\ \n")
cat("$s_{100}$          ",inttemp(l1[13,],u1[13,])," & --- & --- & --- & --- ","\\\\ \n")
cat("$s_{011}$          ",inttemp(l1[14,],u1[14,])," & --- & --- & --- & --- ","\\\\ \n")
cat("$s_{101}$          ",inttemp(l1[15,],u1[15,])," & --- & --- & --- & --- ","\\\\ \n")
cat("$s_{110}$          ",inttemp(l1[16,],u1[16,])," & --- & --- & --- & --- ","\\\\ \n")
cat("$s_{111}$          ",inttemp(l1[17,],u1[17,])," & --- & --- & --- & --- ","\\\\ \n")
cat("\n")
sink()





#END