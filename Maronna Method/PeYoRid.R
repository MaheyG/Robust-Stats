library(MASS)
###########PeYoRid function##############
#[betamin resid sigma edf pesos]=PeYoRid(X,y,lam,deltaesc,nkeep,niter,tolcrit,tolres)
#Peña-Yohai for RR and M-scale S-estimator
#X and y assumed normalized. 
#Intercept added at the end
#niter= #final iterations; 
#tolcrit, tolres=tolerances for relative change in criterion and residuals
#pseos= final weights

PeYoRid<-function(X,y,lam,deltaesc,nkeep=5,niter=50,tolcrit=1.e-3,tolres=1.e-4){
  factr=facon(deltaesc) #corrects lambda to match RR-LS
  lamfac=factr*lam
  n=dim(X)[1]; p=dim(X)[2]
  Xuno=cbind(X,matrix(1,n,1))
  singu=0 #If too much colinearity choose singu=1
  
  temp=regrid(X,y,lam); betals=temp$beta; Xau=temp$Xau ; yau=temp$yau #RR-LS
  reslau=yau-Xau%*%betals
  temp=svdecon(Xau) ; A=temp$a ; sval=temp$s ; B=temp$b
  h=matrix(apply(A^2,1,sum)) 
  w=(reslau/(1-h))^2
  ab=A%*%t(B) ; Q=t(ab)%*%diag(c(w))%*%ab
  temp=eigen(Q) ; U=temp$vectors ; d=temp$values
  Z=A%*%t(B)%*%U; Z=Z[1:n,] #Take out added "observation
  
  sigls=mscale_2(reslau[1:n],0,deltaesc)
  critls=n*sigls^2+lamfac* sum(betals[1:p]^2)^2
  critkeep=critls ; betakeep=betals
  
  #Prosac=proportion of extrem residuals to be omitted
  prosac=deltaesc
  m=ceiling(prosac*n+1.e-6) #Strange
  n1=n-m ; cuales=0
  lam1=lam*n1/n #recall that there are only n1 "actual" observation
  
  #Loop which can be optimized
  #To process the Pena Yohai procedure
  for(j in 1:dim(Z)[2]){ 
    zj=Z[,j] ; zjord=sort(zj) ; ii=order(zj)
    Xj=X[ii,] ; yord=y[ii]
    
    Xjj=Xj[1:n1,] ; yj=matrix(yord[1:n1]) #Higher
    betaj=regrid(Xjj,yj,lam1)$beta ; resj=y-Xuno%*%betaj
    sigj=mscale_2(resj,0,deltaesc)
    critj=n*sigj^2+lamfac*sum(betaj[1:p]^2)^2
    critkeep=rbind(critkeep,critj) ; betakeep=cbind(betakeep,betaj)

    Xjj=Xj[(m+1):n,] ; yj=matrix(yord[(m+1):n]) #Lower
    beta=regrid(Xjj,yj,lam1)$beta ; resj=y-Xuno%*%betaj
    sigj=mscale_2(resj,0,deltaesc)
    critj=n*sigj^2+lamfac*sum(betaj[1:p]^2)^2
    critkeep=rbind(critkeep,critj) ; betakeep=cbind(betakeep,betaj)
    
    zjabord=sort(abs(zj)) ; ii=order(zj)#Higher absolute value
    yord=y[ii]
    Xj=X[ii,]
    Xjj=Xj[1:n1,] ; yj=matrix(yord[1:n1]) 
    betaj=regrid(Xjj,yj,lam1)$beta ; resj=y-Xuno%*%betaj
    sigj=mscale_2(resj,0,deltaesc)
    critj=n*sigj^2+lamfac*sum(beta[1:p]^2)^2
    critkeep=rbind(critkeep,critj) ; betakeep=cbind(betakeep,betaj)
    nk=dim(critkeep)[1]
    if( nk>nkeep){ #Omit unneeded candidates
      temp=unitol(critkeep); critor=temp$y ; ii=temp$ind #Omit repeated ones and sort
      betakeep=betakeep[,ii]
      nuni=min(nkeep,length(ii))
      critkeep=matrix(critor[1:nuni]) ; betakeep=betakeep[,1:nuni]
    }
  }
  #Descent for M-scale and S-estimator
  critmin=Inf
  for(k in 1:nkeep){
    betak=betakeep[,k]
    temp=desrobrid(betak,lamfac,X,y,deltaesc,niter,tolcrit,tolres)
    betak=temp$beta ; resk=temp$res ; sigmak=temp$sig ; edfk=temp$edf ; pesok=temp$w ; critk=temp$crit
    if(critk<critmin){
      critmin=critk ; sigma=sigmak ; betamin=betak ; resid=resk ; edf=edfk ; pesos=pesok
    }
  }
  
  deltult=0.5*(1-(edf+1)/n) #"delta" for MM
  deltult=max(deltult,0.25)
  #c0 =constant for consistency
  c0=7.8464-34.6565*deltult + 75.2573*deltult^2 -62.5880*deltult^3
  sigma=mscale_2(resid,0,deltult)/c0
  q=(edf+1)/n ; k1=1.29 ; k2=-6.02 ; corrige=1/(1-(k1+k2/n)*(q/n)) 
  sigma=sigma*corrige #biais correction
  result=list(betamin=betamin,resid=resid,sigma=sigma,edf=edf,pesos=pesos)
  return(result)
}

regresing<-function(X,y,singu=0){ #Don't called
  if(singu>0){beta=ginv(X)%*%y} #Moore-Penrose Pseudoinverse given by package MASS
  else {beta=solve(X,y)}
}

regrid<-function(X,y,lam){##RR-LS 
  n=dim(X)[1] ; p =dim(X)[2]
  Xau=rbind(cbind(X,matrix(1,n,1)),cbind(sqrt(lam)*diag(p),matrix(0,p,1)))
  yau=rbind(y,matrix(0,p,1)) 
  beta=matrix(lsfit(Xau,yau,intercept = FALSE)$coeff) #Compue Least Square soluton of the system
  result=list(beta=beta,Xau=Xau,yau=yau)
  return(result)
}

facon<-function(delta){#Factor that corrects lambda to match RR-LS
  ff=23.9716 -73.4391*delta+ 64.9480*delta^2
  return(ff)
}
