####################################################            
#ROBRIDGE [beta resid edf lambest]= RobRidge(X,y, numlam,cualcv, showhist,nkeep)
#Solves n*sig^2 *sum{rho(resid/sig)+lam*||beta1||^2} = min
# Required input: X, y data set
# Optional input:
#numlam: number of lambda values, default =min(n,p,20)
#cualcv: method for estimating prediction error. If cualcv=0 (default): approximate leave-one-out CV; 
#if >0: actual cualcv--fold CV ("N_{lambda}";
#if 0<cualcv<1: random test sample of size n*cualcv
#showhist: if >0, print edf and mse for each lambda (default=0)
#nkeep= number of candidates to be kept for full iteration in the Peña-Yohai procedure (default=5)
#Output
#beta= (p+1)-vector of regression parameters, %beta(p+1)=intercept
#resid= residual vector
#edf= final equivalent degrees of freedom ("p hat")
#lambest= optimal lambda

RobRidge<-function(X,y,numlam,cualcv=2,showhist=FALSE,nkeep=5,ncore=1,plot.path=FALSE,bestlambda=c("1se","min"),
	measure=c("robust mse","mae","trimmed mse","mse"),fold.id){ 
  n=dim(X)[1]; p=dim(X)[2]
  
  if(missing(numlam)){numlam=min(n,p,20)}
  
  #Normalize and center X and y
  temp=prepara(X,y)
  Xnor=temp$Xnor ; ynor=temp$ycen ; mux=temp$mux ; sigx=temp$sigx ; muy=temp$muy
  
  #Spherical Principal Components (no centering)
  #privar, Beig= vector of robust "eigenvalues" and matrix of eigenvectors
  #Xnor is now =PCA scores= "orthonormalized Xnor "
  temp=SPC(Xnor,0)
  privar=temp$lambda ; Beig=temp$b ; muspam=temp$mu ; Xnor=temp$scores
  n=dim(Xnor)[1]; p=dim(Xnor)[2] #p is now the 'actual' dimension
  privar=privar*n #Makes the robust eigenvalues of the same order as those of classifcal PCA used for LS
  #Creation of the lambda path
  nlam=min(p,numlam)
  pmax=min(p,n/2) #edf<=n/2 to keep BDP >=0.25
  pp=seq(1, pmax, length.out = nlam) #candidate edf'
  lamdas=findlam(privar,pp) #find lambdas corresponding to the edf
  lamdas[lamdas==0]=0.01
  
  deltas=0.5*(1-pp/n) #for the M-scale used with Pena Yohai

  #Cross Validation
  historia=NULL
  plot.mse<-NULL
  if(missing(fold.id)){fold.id<-ceiling(sample(1:n)/n*cualcv)}
  for(klam in 1:nlam){
    lam=lamdas[klam] ; deltaesc=deltas[klam]
    cat("Start Cross Validation for lambda =",lam,'#',klam,"in",nlam,": \n") 
    mse=CVRidRob(Xnor,ynor,cualcv,fold.id,lam,pp[klam],ncore)
    historia=rbind(historia,cbind(lam,mse))
  }
  
  if(match.arg(measure)=="mse"){r=2}
  if(match.arg(measure)=="mae"){r=3}
  if(match.arg(measure)=="trimmed mse"){r=4}
  if(match.arg(measure)=="robust mse"){r=5}
  historia=cbind(historia,pp) 
  
  #Select the well suited lambda with criterium match.arg(measure)
  i=which.min(historia[,r])
  mse.min=historia[i,r] ; lam.min=lamdas[i] ; delta.min=deltas[i]
  j = min(which(historia[,r] < mse.min+historia[i,6]))
  mse.lse=historia[j,r] ; lam.lse=lamdas[j] ; delta.lse=deltas[j]
  
  #Best lambda from CV
  if(match.arg(bestlambda)=="lse"){lambest=lam.lse ; delmin=delta.lse}
  else{lambest=lam.min ; delmin=delta.min} 
  
  #Plot the result of Cross Validation
  plot_cv(historia,lambest)
  
  #Plot Path
  if(plot.path){
    cat("Start Computing the path for",nlam, "lambda \n")
    plot_path(Xnor,ynor,lamdas,pp,ncore,lambest)
  }
  
  #RidSEMM for best lambda
  temp=RidSEMM(lambest,Xnor,ynor,delmin,nkeep)
  betaSE=temp$betaSE ; beta=temp$betaMM ; residSE=temp$residSE ; resid=temp$residMM
  edfSE=temp$edfSE ; edf=temp$edfMM
  
  #Denormalize beta
  betaslo=Beig%*%beta[1:p] ; bint=beta[p+1]
  beta=desprepa(betaslo,mux,sigx,muy+bint)
  result<-list(beta=beta,resid=resid,edf=edf,lambest=lambest,lamdas=lamdas)
  return(result)
}

  
