#######################################
######Function RidSEMM
#[betaSE betaMM residSE residMM edfSE edfMM mseMM]=RidSEMM(X,y,lam,deltaesc,nkeep,keff)
#Computes RR-SE and RR-MM for a given lambda=lam
#X,Y= data, X is assumed normalized
#deltaesc (optional)= "delta" for scale M-estimate (default: =0.5)
#nkep (optional)= # candidates for full iteration in Peña-YOhai (default=5)
# keff (optional): constant "c" for efficiency from lilst below:
  # efipos=[0.8 0.85 0.9 0.95]; %possible efficiencies
#   keff= [3.14 3.44 3.88 4.68]; %sus constantes
# default: efficiemcy = nominal 0.85 adjusted by edf
# betaSE, betaMM=estimates, with intercept at the end
# residSE, residMM= residuals; edfSE, edfMM= equivalent d.f.
# mseMM= estimated robust MSE

RidSEMM<-function(lam,X,y,deltaesc=0.5,nkeep=5,keff){
  
  #Compute Pena Yohai, S-estimator and M-scale estimator
  temp=PeYoRid(X,y,lam,deltaesc,nkeep)
  betaSE=temp$betamin ; residSE=temp$resid ; sigma=temp$sigma
  edfSE=temp$edf
  
  if(missing(keff)){#Correction
    psn=edfSE/length(y)
    if(psn<0.1){keff=3.5}
    else if(psn<0.2){keff=3.77}
    else if (psn<0.33){keff=4}
    else {keff=4}
  }
  
  #Compue MM estimator 
  temp=MMRid(X,y,lam,betaSE,sigma,keff)
  betaMM=temp$beta ; residMM=temp$res ; edfMM=temp$edf ; w=temp$w ; kasig=temp$kasig
  result<-list(betaSE=betaSE,betaMM=betaMM,residSE=residSE,residMM=residMM,edfSE=edfSE,
               edfMM=edfMM,kasig=kasig,sigma=sigma)
  return(result)
}
  
