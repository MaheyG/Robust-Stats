######################################
####Function Prepara########## 
#Centers y and the columnds of X to zero location and normalize X to unit scale
# If robu>0 (default) : with robust M-estimates; else means and SD
#Xnor=centered and normalized X ; ycen= centered y ; mux=location vector  of x
#muy  ;sigy=location and scale of y 

prepara<- function(X,y,robu=TRUE){
  p=dim(X)[2] #  n=dim(X)[1] ; 
  if(robu){
    mux=Mloca(X)
    Xnor=centrar(X,mux) 
    sigx=matrix(0,1,p)
    muy=Mloca(y)
    ycen=y-muy
    sigy=mscale_2(ycen)
    for(i in 1:p){sigx[i]=mscale_2(Xnor[,i])}
  }
  else{
    mux=mean(X)
    Xnor=centrar(Xn,mux)
    sigx=colSds(X)
    muy=mean(y)
    sigy=sd(y)
    ycen=y-muy
  }
  Xnor=divcol(Xnor,sigx)
  result<-list(Xnor=Xnor, ycen=ycen, mux=mux,sigx=sigx,muy=muy,sigy=sigy)
  return(result)
}
