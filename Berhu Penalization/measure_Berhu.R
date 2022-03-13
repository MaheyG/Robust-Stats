#Compute MSE or MAE or Deviance criterium for cross validation
measure.berhu <- function(y, yhat, args) {
  r <- y-yhat
  type.measure <- args$type.measure

  if (type.measure == "deviance") {
    method <- args$method
    if(method=="Huber"){
      sigma<-args$sigma
      M<-args$M
      #val <- hloss(r,M) #RECOMPUTE WITH THE SIGMA
      val <- hloss_sigma(r,M,sigma)}
    else {val<-r^2}
  } else if (type.measure == "mse") {
    val <- r^2
  } else {
    val <- abs(r)
  }
  val
}

hloss_sigma <- function(r,M,sigma) {
  rr <- abs(r)/sigma
  k<-dim(r)[2]
  for (i in 1:k){
    if(rr[i]<=M){
      rr[i]<-rr[i]^2/(2*M)
    }else{
      (rr[i]-M/2)
    }
  }
  rr<-rr*sigma
  return(rr)
}

hloss <- function(r, gamma) {
  rr <- abs(r)
  ifelse(rr <= gamma, rr^2/(2*gamma), rr-gamma/2)
}
