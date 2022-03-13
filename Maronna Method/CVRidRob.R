####################################################
###Function CVRidRob#########"
#mse=CVRidRob(XX,yy,nfold,orden,lam,gradlib)
#XX,yy= data
#performs nfold-CV 
#gradlib= degrees of freedom = edf
#orden= vector of random ordering
#resulting robust MSE or MSE or trimmed MSE or MAE

CVRidRob<-function(XX,yy,nfold,fold.id,lam,gradlib,ncore){
  
  n=dim(XX)[1]
  nestim=n*(1-1/nfold)
  deltaesc=0.5*(1-gradlib/nestim)
  inint=ceiling(seq(0,n,length.out=nfold+1))
  resid=matrix(0,n,1)
  
  if(ncore==1){ #No parallelism
    for(kk in 1:nfold){
      cat("folder #",kk," in ",nfold,"\n")
      resid<-cross_validation(kk,XX,yy,fold.id,lam,resid,deltaesc)
    }
  }
  else { #Parallelisme
    cat("Start Parallel Cross Validation for",nfold,"folds","\n")
    max.cores <- detectCores()
    if (ncore > max.cores) {
      cat("The number of cores specified (", ncores, ") is larger than 
          the number of avaiable cores (", max.cores, "), so", max.cores, "cores are used.", "\n")
      ncore = max.cores
    }
    cluster <- makeCluster(ncore,type="FORK")
    fold.results <- parSapply(cl = cluster, X = 1:nfold, FUN=cross_validation,
                              XX=XX,yy=yy,fold.id=fold.id,lam=lam,resid=resid,deltaesc=deltaesc)
    stopCluster(cluster)
    resid<-matrix(rowSums(fold.results))
  }
  return(compute_mse(resid))
  
}

cross_validation<-function(kk,XX,yy,fold.id,lam,resid,deltaesc){
  Xtrain<-XX[fold.id !=kk,,drop=FALSE]
  ytrain<-matrix(yy[fold.id !=kk])
  Xtest<-XX[fold.id == kk,,drop=FALSE]
  ytest<-matrix(yy[fold.id == kk])
  
  temp=RidSEMM(lam,Xtrain,ytrain,deltaesc)
  betaMM=matrix(temp$betaMM) ; r=dim(betaMM)[1]
  beta=matrix(betaMM[1:r-1]) ; bint=betaMM[r]
  fitk=Xtest%*%beta+bint
  resid[fold.id==kk,]=ytest-fitk
  return(resid)
}
compute_mse<-function(resid){
  n=dim(resid)[1]
  mse<-mean(resid^2)
  mae<-mean(abs(resid))
  deviance=tauscale(resid,ktau=5)^2
  colnames(deviance)=c("robust_mse")
  resid<-abs(resid)
  trimmed_resid=resid[resid<=quantile(resid,prob=c(0.9))]
  trimmed_mse=mean(trimmed_resid^2)
  cvse=sd(resid^2)/sqrt(n)
  return(cbind(mse,mae,trimmed_mse,deviance,cvse))
}
