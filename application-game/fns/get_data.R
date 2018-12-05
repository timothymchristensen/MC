# 
# Imports data for entry game application
# 
get_data <- function(){
  #
  main_data <- as.matrix(read.table("airlinedata.dat", sep = ','))
  #
  colnames(main_data) <- c("mdist","mp_LC","mp_OA","msize","LC","OA")
  #
  for(i in c("mp_LC","mp_OA","msize")){
    main_data <- cbind(main_data,1*(median(main_data[,i]) < main_data[,i]))
    colnames(main_data)[ncol(main_data)] <- paste0(i,"_d")
  }
  #
  X_OA <- cbind(1, main_data[,c(9,8)]) # intercept, market size, marketpresence OA
  X_LC <- cbind(1, main_data[,c(9,7)]) # intercept, market size, marketpresence LC
  #
  # construct dependent variable
  D <- cbind(as.character(main_data[,6]),as.character(main_data[,5])) # OA, LC
  D <- matrix(apply(D,1,function(x) paste0(x[1],x[2])),ncol=1)
  D <- 1*t(apply(D,1,function(x) x==c("00","01","10","11"))) #(0,0), (0,1), (1,0), (1,1)
  colnames(D) <- c("00","01","10","11")
  #
  # construct regressor index for selection probability
  s_ind <- cbind(main_data[,9],main_data[,8],main_data[,7])
  s_ind <- matrix(apply(s_ind,1,function(x) paste0(x[1],x[2],x[3])),ncol=1)
  s_ind <- 1*t(apply(s_ind,1,function(x) x==c("000","001","010","100","011","101","110","111")))
  colnames(s_ind) <- c("000","001","010","100","011","101","110","111")
  #
  mdata <- list(X_LC=X_LC,X_OA=X_OA,D=D,s_ind=s_ind)
  #
  mdata_short <- reg_short <- matrix(NA,8,ncol(D))
  for(i in 1:8){
    mdata_short[i,] <- colSums(D[s_ind[,i]==1,])
    reg_short[i,] <- c(colMeans(X_OA[s_ind[,i]==1,2:3]),colMeans(X_LC[s_ind[,i]==1,2:3]))
  }
  mdata_short <- as.numeric(mdata_short)
  #
  return(list(mdata=mdata,mdata_short=mdata_short,reg_short=reg_short))
}
#
# END