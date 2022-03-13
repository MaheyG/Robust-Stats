library(CVXR)
library(hqreg)
library(scs)
library(parallel)
#library(care)  ; data("efron2004")# to get the data of Efron
#par(mfrow=c(1,2))

#Replace the ... by the way of your own computeur to charge functions in your environment
# source("~/Documents/M2 - Math En Action/Stage - Robust/TP Owen CVXR/Berhu.R")
# source("~/Documents/M2 - Math En Action/Stage - Robust/TP Owen CVXR/Berhu_Parallel.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP Owen CVXR/Cross_Validation_Berhu.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP Owen CVXR/measure_Berhu.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP Owen CVXR/predict_Berhu.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP Owen CVXR/Plot_Berhu.R")
# source("~/Documents/M2 - Math En Action/Stage - Robust/TP Owen CVXR/Adaptive_Berhu.R")
# source("~/Documents/M2 - Math En Action/Stage - Robust/TP Owen CVXR/Adaptive_Berhu_Parallel.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP Owen CVXR/Regression_ls_huber.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP Owen CVXR/Berhu_warmstart.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP Owen CVXR/Adaptive_Berhu_warmstart.R")

# 
#  n=442
#  p=10
#  X<-matrix(efron2004$x[1:n,],n,p)
#  y<-matrix(efron2004$y[1:n])




### Path with Berhu/Adaptive Berhu penalization with Huber/Least square regression ###

# path=Berhu_warmstart(X,y,pena="Berhu",loss="Huber",nlambda = 100,L=2,lambda_max =0.35 )
# plot.path.Berhu(path)
# 
# path=Berhu_warmstart(X,y,pena="Berhu",loss="Huber",nlambda = 100,L=2,lambda_max =0.35,max_iters=500L )
# plot.path.Berhu(path)
# 
# path=Adaptive_Berhu_warmstart(X,y,loss="Huber",beta_ini ="unpenalized",nlambda = 20,L=2)
# plot.path.Berhu(path)
# 
# 
# ### Cross Validation with Berhu/Adaptive Berhu penalization with Huber/Least square regression ###
# cv<-cv_Berhu(X,y,FUN="Berhu",nlambda=100 ,seed=123,nfolds = 5,ncore=5,type.measure="deviance",L=2)
# plot.cv.Berhu(cv)
# 
# cv<-cv_Berhu(X,y,FUN="Berhu",nlambda=100 ,nfolds=5, ncore = 5,type.measure="mse",seed=123)
# plot.cv.Berhu(cv)
# 
# cv<-cv_Berhu(X,y,FUN="Adaptive_Berhu",nlambda=100,nfolds=5,ncore=5,type.measure = "mse",seed=123,lambda_max=0.35)
# plot.cv.Berhu(cv)
# 
# #Path and Cross validation with hqreg package
# 
# fit= hqreg(X, y,nlambda = 100,alpha=0,method="huber")
# plot(fit) 
# 
# 
# cv=cv.hqreg(X,y,method="huber",type.measure="mse",seed=123,FUN="hqreg_raw", nfolds = 5,alpha=1)
# plot(cv)

