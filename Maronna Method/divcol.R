######################################
####Function Divcol########## 
#Divide columns of X by sig or by their standart deviation (default)

divcol<-function(X,sig){
  n=dim(X)[1]; p=dim(X)[2]
  if(missing(sig)){
    sig=matrix(apply(X,2,sd),1,p)
  }
  return(X/(matrix(1,n,1)%x%sig))
}
