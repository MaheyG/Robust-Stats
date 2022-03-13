######################################
####Function Centrar########## OK
#Center the matrix X by means 0 (default) or by row vector mu

centrar<-function(X,mu){
  n=dim(X)[1]
  if(nargs()<2){
    mu=colMeans(X)
  }
  return(X-matrix(1,n,1)%x%t(mu)) #Kronecker product
}
