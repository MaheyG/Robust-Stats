######################################
####Function SPC########## 
#Compute the SPC of X, 
#cent is the centeration what dim? ? If cent>0 (default) x is centered if not "is used for RobRidge?"
# Given lambda robust eigen value (pseudo eigen value), b matrix of eigenvector, mu spatial median 
#and score projection of X (centered) on the eigen vectors

SPC<-function(X, cent=1){ 
  n=dim(X)[1];p=dim(X)[2]
  
  temp=spamed(X,cent) 
  mu=temp$mu ; w=temp$w
  
  xcen=centrar(X,mu) 
  
  y=xcen*(w*rep(1,p)) 
  
  temp=svdecon(y) #Economic SVD
  a=temp$a 
  s=temp$s 
  b=temp$b
  scores=xcen%*%b
  lambda=robsq(scores) #lambda=squared robust scale along principal directions
  
  
  #sort lambdas increasing
  ind=order(lambda) #give the position for each element to be an increasing vector
  lambda=sort(lambda)
  b=b[,ind] ; scores=scores[,ind]
  result<-list(lambda=lambda,b=b,mu=mu,scores=scores)
  return(result)
}


spamed<-function(x,cent){#Compute the weights and the spatial median
  # w=weights ; 
  #if cent>0; mu=spatial median else mu=0
  #w=1/||x-mu|| normalized 
  n=dim(x)[1];p=dim(x)[2]
  del0=1.e-5
  if (cent>0){
    niter=20 #Nb of iteration
    tol=1.e-5 
    mu0=median(x) 
    dife=Inf 
    ite=0 #Number of iteration done
    while(ite<niter && dife>tol*p){ #Compute the spatial median by iteration
      ite=ite+1
      xc=centrar(x,mu0) 
      w=sqrt(sum(xc^2)) 
      deldis=del0*median(w) 
      w=w.*(w>=deldis)+deldis*(w<deldis); #truncate just in case(divide by zero)-> coordinate more little than deldis are put on deldis
      w=1/w
      mu=sum((w*ones(1,p))*x)/sum(w) 
      dife=sqrt((mu0-mu)*(mu0-mu)) 
      mu0=mu
    }
  }
  else{
    mu=matrix(0,p,1)
    w=sqrt(apply(x^2,2,sum))
    deldis=del0*median(w)
    w=w*(w>=deldis)+deldis*(w<deldis) 
    w=1/w
  }
  w=w/sum(w);
  result<-list(mu=mu,w=w)
  return(result)
}

robsq<-function(z){#Squared Dispersion
  n=dim(z)[1];p=dim(z)[2]
  mu=apply(z,2,median)
  zcen=centrar(z,mu)
  v=matrix(0,p,1)
  for(k in 1:p){
   v[k]=mscale_2(zcen[,k])^2
    } 
  return(v)
}


