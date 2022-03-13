##########################################################################################
#DESROBRID [beta res sig edf w crit] =desrobrid(betin,lam,x,y,delsca,niter,tolcrit,tolres)
#Descent for RR-SE starting from initial estimate betin
#x,y= data,  lam= penalty
#Intercept always at the end
#delsca= "delta" of M-scale
#niter: max. iterations, default=100
#tolcrit,tolres: tolerances for criterion and residuals,
#         defaults 1.e-3 and1.e-4
#w=weights; lala=lambda modified by w
#crit=criterion

desrobrid<-function(betin,lam,x,y,delsca,niter=100,tolcrit=1.e-3,tolres=1.e-4){
  n=dim(x)[1]; p=dim(x)[2]
  betinslo=betin[1:p] ; betinte=betin[p+1]
  res0=y-x%*%betinslo-betinte
  sig0=mscale_2(res0,0,delsca)
  iter=0; delta=Inf ; conve=Inf
  crit0=n*sig0^2+lam*t(betinslo)%*%betinslo
  binter=betinte
  while(iter<niter &&(delta>tolcrit || conve>tolres)){
    iter=iter+1
    tt=res0/sig0
    w=weights(tt) ; rw=sqrt(w)
    ycen=y-binter
    xw=x*(matrix(1,1,p)%x%rw)
    yw=ycen*rw
    lala=mean(w*tt^2)*lam 
    xau=rbind(xw,sqrt(lala)*diag(p))
    yau=rbind(yw,matrix(0,p,1))
    beta=matrix(lsfit(xau,yau,intercept = FALSE)$coeff) 
    resin=y-x%*%beta
    binter=sum(resin*w)/sum(w)
    res=resin-binter
    sig=mscale_2(res,0,delsca)
    crit=n*sig^2+lam*t(beta)%*%beta ; deltold=delta
    delta=1-crit/crit0
    conve=max(abs(res-res0))/sig #Measure convergence of residuals
    res0=res ; sig0=sig ; crit0=crit
  }
  beta=rbind(beta,binter)
  hmat=xau%*%lsfit(t(xau)%*%xau,t(xau),intercept = FALSE)$coeff
  h=diag(hmat) ; edf=sum(h[1:n])
  resp<-list(beta=beta,res=res,sig=sig,edf=edf,w=w,crit=crit)
  return(resp)
}

weights<-function(r){ #w=(psibis(r)/r)/2 to eliminate the 2 from 2*lam
  w=(abs(r)<=1)*(1-r^2)^2
}