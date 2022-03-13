####Function tauscale#############
# tau scales (row vector) of x for several constants ktau (row)
#compute the robuste mse based on a tau-scale estimator

tauscale<-function(x,ktau,delta=0.5){
  sigmas=NULL
  s0=mscale_2(x,0,delta)
  c0=7.8464-34.6565*delta + 75.2573*delta^2 -62.5880*delta^3
  s0=s0/c0
  for(k in ktau) {
    romed=mean(rho(x/(s0*k)))
    sig=k*s0*sqrt(romed)
    sigmas=cbind(sigmas,sig)
  }
  return(sigmas)
}

rho<-function(x){ #Bisquare for rho''(0)=2
  return((1-(1-x^2)^3*(abs(x)<=1))/3)
}
