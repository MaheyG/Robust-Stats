##########################
###function Mloca
#mu=Mloca(y,sig,efi)
#mu (column vector) = columnwise bisquare location M-estimator of matrix y. 
#sig (optional): scale, either a single number or a column vector.
#       Default=Normalized columnwise MAD
#efi (optional): desired efficiency, Choose 0.9 (default), 0.85 or 0.95

Mloca<-function(y,sig,efi=0.9){
  if(missing(sig)){sig=mad(y)/0.675}
  
  if(efi==0.85) {kefi=3.44}
  else if (efi==0.9){kefi=3.88}
  else if (efi==0.95){kefi=4.68}
  else {disp("wrong efi in Mloca")}
  q=dim(y)[1]; n=dim(y)[2]
  niter=10
  sig<-matrix(sig)
  nsig=dim(sig)[1]
  if(nsig==n){sigrep=t(sig)%x%matrix(1,q,1)}
  else {sigrep=sig%x%matrix(1,q,n)}
  if(q==1){mu=y}
  else {
    mume=apply(y,2,median) #initial
    mu=mume
    for (j in 1:niter){
      z=(y-t(mu%x%matrix(1,1,q)))/sigrep
      w=bisq(z/kefi)
      sw=colSums(w) ; nulos=(sw==0)
      if(sum(nulos)==0){mu=sum(y*w)/sw}
      else{nonul=c(nonul=!nulos) #to avoid division by 0
      mu[nonul]=sum(y[,nonul]*w[,nonul])/sw[nonul]
      mu[nulos]=mume[nulos]
      } 
    }
  }
  return(mu)
}

bisq<-function(z){#bisquare weight function
  t=z^2
  w=(t<=1)*(1-t)^2
  return(w)
}
