######################################
####Function Despreparation########## A Vérifier avec données
#Inverse function of prepa return to original scale

desprepa<-function(beta0,mux,sigx,muy){
  beta<-beta0/t(sigx)
  bettint=muy-t(matrix(mux))%*%beta 
  beta=rbind(beta,bettint)
  return(beta)  
}