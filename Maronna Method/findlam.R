########################################################"
###function lamr= findlam(vals,r)###
#FINDLAM lamr= findlam(vals,r)   column vector
#Finds lamdas which yield edf=r
#Resolve f(lambda)=0  as in (22)
findlam<-function(vals,edf){
  p=length(vals); nr=length(edf)
  lamr=matrix(0,nr,1)
  lam1=0 ; lam2=max(vals)*(p/edf-1)
  for(i in 1:nr){
    lam0=cbind(lam1,lam2[i]+0.5)#the value 0.5 is for lam=0
    lamr[i]=uniroot(f=sumrat,interval = lam0, vals=vals,edf=edf[i])$root
  }
  return(lamr)
}
sumrat<-function(lam,vals,edf){
  sum(vals/(vals+lam))-edf 
}
