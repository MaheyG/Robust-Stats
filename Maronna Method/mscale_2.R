######################################
####Function M-scale ########## 
#sigma solution of sumrho(x_i/sigma)=delta  where rho is bisquare 
#delta: optional, default=0.5 
#tole : optional, error tolerance,  default=1.e-5
#normz: optional; if >0 (default), normalize sig for consistency at the normal


mscale_2<-function(x,normz=1,delta=0.5,tole=1.e-5){

  x<-matrix(x)
  n=dim(x)[1]; p=dim(x)[2]
  if(p>n){x<-t(x)} 
  y=apply(abs(x),2,sort)

  n1=floor(n*(1-delta)) #floor : lower integer
  n2=ceiling(n*(1-delta)/(1-delta/2)) #ceiling : upper integer
  qq=cbind(y[n1],y[n2])

  u=rhoinv(delta/2)
  sigin=cbind(qq[1],qq[2]/u)

  if(qq[1]>=1){tolera=tole} #relative or absolute tolerance for sigma >or< 1
  else{tolera=tole*qq[1]}
  
  if(mean(x==0)>1-delta){sig=0}
  else{#MINIMIZATION of the function
    sig=uniroot(f=averho,interval=sigin,x=x,delta=delta,tol=tolera)$root
  }
  if(normz>0){sig=sig/1.56} #Normalize
  sig
} 

rhobisq<-function(x){1-(1-x^2)^3*(abs(x)<=1)} #Bisquare
averho<-function(x,sig,delta){mean(rhobisq(x/sig))-delta} #Function to put on 0
rhoinv<-function(u){return(sqrt(1-(1-u)^(1/3)))} #Inverse function of rho -> Ok for x in [0,1/2]
