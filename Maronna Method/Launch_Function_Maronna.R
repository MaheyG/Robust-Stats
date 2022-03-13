library(matrixStats)
#library(pense)
#library(profvis)
library(parallel)
library(Matrix)

#Change the ... by the path of your computer to load the function on the environement
if(TRUE){
source("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/R Maronna/Maronna_Changed/RobRidge.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/R Maronna/Maronna_Changed/centrar.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/R Maronna/Maronna_Changed/CVRidRob.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/R Maronna/Maronna_Changed/desprepa.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/R Maronna/Maronna_Changed/Desrobrid.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/R Maronna/Maronna_Changed/divcol.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/R Maronna/Maronna_Changed/findlam.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/R Maronna/Maronna_Changed/Mloca.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/R Maronna/Maronna_Changed/MMRid.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/R Maronna/Maronna_Changed/mscale_2.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/R Maronna/Maronna_Changed/PeYoRid.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/R Maronna/Maronna_Changed/prepara.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/R Maronna/Maronna_Changed/RidSEMM.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/R Maronna/Maronna_Changed/SPC.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/R Maronna/Maronna_Changed/SVDEco.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/R Maronna/Maronna_Changed/tauscale.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/R Maronna/Maronna_Changed/unitol.R")
source("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/R Maronna/Maronna_Changed/Plot_Path_CV.R")  
}


 # #Data Effron
 # library(care)  ; data("efron2004")# to get the data of Efron
 # X<-matrix(efron2004$x[1:442,],442,10)
 # y<-matrix(efron2004$y[1:442])
# 
# 
# #Data Vessel
#  X<-read.delim("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/Vessel_X.txt",header=FALSE,sep=",",dec=".")
#  X<-matrix(unlist(X),180,301)
# Y<-read.delim("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/Vessel_Y.txt",header=FALSE,sep=",",dec=".")
#  Y<-matrix(unlist(Y),180,13)
#  y<-matrix(Y[,1])

# #Maronna's method
# a<-RobRidge(X_train,y_train,numlam = 5 ,cualcv = 5,ncore=5,plot.path=F,measure="robust mse")

# l=a$lamdas
# 
# #LS Ridge and  Huber Ridge
# library(hqreg)
# fit<-hqreg_raw(Xnor,ynor,alpha=0,method = "ls",lambda=l)
# plot(fit)
# cv=cv.hqreg(Xnor,ynor,FUN="hqreg_raw",type.measure="mse",nfolds = 5,method="ls",alpha=0,lambda=l)
# plot(cv)
# 
# #Pense Package (very slow)
# library(parallel)
# cv_results <- pense_cv(X, y, alpha = 0, cv_k = 5,ncores=5,nlambda=100)
# plot(cv_results)
# path<-pense(X,y,alpha=0,nlambda=100,ncores=5)
# plot(path)
