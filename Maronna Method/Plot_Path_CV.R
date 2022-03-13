#Plot Path and CV
#Compute and plot the path for the selecting lambda

plot_path<-function(X,y,lam,edf,ncore,lamin){
  n=dim(X)[2] ; p=length(lam)
  if(ncore==1){
    betaMM<-sapply(lam,RidSEMM,X,y)[2,]
  }

  else{
    clust <- makeCluster(ncore,type="FORK")
    betaMM<- parSapply(clust,lam,RidSEMM,X,y)[2,]
    stopCluster(clust)
  }
  betaMM=matrix(unlist(betaMM),n+1,p)
  lambda=log(matrix(lam))
  plot(0, 0, typ="n",main = "",
       xlab = expression(log(lambda)), ylab = expression(beta), ylim=range(betaMM),
       xlim = range(lambda))
  matlines(lambda, t(betaMM[0:n,]),lwd=1,type='l',lty=1)
  axis(3,at=lambda,labels=round(edf)) ; mtext("edf",side=3,line=2)
  
  #best lambda choosen
  abline(v = log(lamin),lty=2,lwd=0.7)
}

plot_cv<-function(historia,lamin){
  lambda=log(historia[,1])
  plot(0,0,xlab=expression(log(lambda)),ylab="",type="n",ylim=range(historia[,2:5]),xlim=range(lambda))
  matlines(lambda,historia[,2:5],lwd=2,type='l',lty=1)
  axis(3,at=lambda,labels=round(historia[,7])) ; mtext("edf",side=3,line=2)
  legend("topleft", legend = colnames(historia)[2:5], col = 1:4, lty = 1,cex=0.7)
  
  #Vertical line , standard error
  L=historia[,2]-historia[,6]
  U=historia[,2]+historia[,6]
  ylim=range(c(L,U))
  ind <- ((U-L)/diff(ylim) > 1e-3)
  arrows(x0=lambda[ind], x1=lambda[ind], y0=L[ind], y1=U[ind], code=3, angle=90, col="gray60", length=.03)
  
  #best lambda choosen
  abline(v = log(lamin),lty=2,lwd=0.7)
}
