######################################
####Function MMRid########## 

#MMRID   [beta res edf w mse]=MMRid(X,y,lam, betin,sigma,kefi,niter,tol)
#RR-MM descent starting from initial estimate "betin", with given scale
#"sigma" and penalty "lam"
#Minimizes criterion (k*sigma)^2*sum{rho(res/k*sigma)}+lamda*||beta1||^2
#Here rho''(0)=2
#niter=#(max. iterations),  default=50
#tol: tolerance for relative change of criterion and residuals, default=1.e-3;
#edf= equiv. deg fr.
#kefi= constant for efficiency, default =3.88 
#For the following efficiencies: %0.8   0.85   0.9   0.95,
# use kefi= 3.14   3.44   3.58   4.68
#res= residuals, w =weights.
#mse=estimated pred. error
#mse(j) for j=1:3 are based on: FPE, CV(n) and GCV(n) 


MMRid<-function(X,y,lam,betin,sigma,kefi=3.88,niter=50,tol=1.e-3){
  n=dim(X)[1]; p=dim(X)[2]
  kasig=kefi*sigma
  betinte=betin[p+1]
  betinslo=betin[1:p]
  res0=y-X%*%betinslo-betinte #Y-XB-B_0
  crit0=kasig^2*sum(rho(res0/kasig))+lam*sum(betinslo^2)
  
  #Iteration
  iter=0;delta=Inf;conve=Inf;
  binter=betinte #Beta of the iteration
  while((iter<niter) && ((delta>tol) ||(conve>tol))){ #atain the max of iteration && convergence criterium
    iter=iter+1
    tt=res0/kasig
    w=weights(tt) 
    rw=sqrt(w)
    ycen=y-binter
    Xw=X*(matrix(1,1,p)%x%rw) # X * Matrix of rw as column(Kronecker product)->normalize the column of X
    yw=ycen*rw
    Xau=rbind(Xw,sqrt(lam)*diag(p)) #Add raw-> Xau in M_(n+p)xp
    yau=rbind(yw,matrix(0,p,1)) #Add raw -> yau in R^(n+p)
    beta=matrix(lsfit(Xau,yau, intercept = FALSE)$coeff) #RESOLVE THE SYSTEM yau=Xaubeta in ls sense
    resin=y-X%*%beta
    if(sum(w)>0){binter=sum(resin*w)/sum(w)
    }
    else {binter=colMedians(resin)
    }
    res=resin-binter #centered residuals
    crit=kasig^2*sum(rho(res/kasig))+lam*sum(beta^2)
    deltold=delta
    delta=1-crit/crit0
    conve=max(abs(res-res0))/sigma #Measure convergence of residuals
    res0=res0;crit0=crit
  }
  beta=rbind(beta,binter)
  hmat=Xau%*%(lsfit(t(Xau)%*%Xau,t(Xau),intercept = FALSE)$coeff) 
  h=diag(hmat)
  edf=sum(h[1:n])
  result<- list(beta=beta,res=res,edf=edf,w=w,kasig=kasig)
  return(result)
}

rho<-function(r){ #To make rho''(0)=2
  r=(1-(1-r^2)^3*abs(r)<=1)/3
}
weights<-function(r){ #w=(psibis(r)/r)/2 to eliminate the 2 from 2*lam
  w=(abs(r)<=1)*(1-r^2)^2
}


