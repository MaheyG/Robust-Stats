
rm(list=ls())


#Data Vessel
X<-read.delim("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/Vessel_X.txt",header=FALSE,sep=",",dec=".")
X<-matrix(unlist(X),180,301)
y<-read.delim("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/Vessel_Y.txt",header=FALSE,sep=",",dec=".")
y<-matrix(unlist(y),180,13)
y<-matrix(y[,1]) # y[,i] for 1<i<13

X_train =X[0:149,]
X_test=X[150:180,]
y_train=matrix(y[0:149,])


# fixing the folds for CV - helpful in reproducing results
# CV is needed in the regularization methods
n <- length(y)
set.seed(123)
nb_folds=5
foldid = sample(rep(seq(nb_folds), length = n))


## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ELASTIC NET  %%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##

library (glmnet)
ptm<-proc.time()

enet.obj <- cv.glmnet(X_train, y_train,  
                      family="gaussian",
                      type="mse", 
                      foldid=foldid,
                      alpha=0.5,   
                      nlambda=100)

lambda.EN    <- enet.obj$lambda.min
beta.EN     <- coef(enet.obj,s=lambda.EN)[-1,1] # to exclude the intercept
intercept.EN <- coef(enet.obj,s=lambda.EN)[1,1]
time.EN=proc.time()-ptm
# predicting for the validation set with the lambda found previously
predict.EN <- predict(enet.obj, newx=X_test, s=lambda.EN)+intercept.EN
## %%%%%%%%%%%%%%%%%%%%%%%%% END ELASTIC NET  %%%%%%%%%%%%%%%%%%%%%%%%%%% ##



## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LASSO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##

library (glmnet)
ptm<-proc.time()
lasso.obj <- cv.glmnet(X_train, y_train,  
                       family="gaussian",
                       type="mse",  
                       foldid=foldid,
                       alpha=1,     
                       nlambda=100)

lambda.LASSO    <- lasso.obj$lambda.min
beta.Lasso      <- coef(lasso.obj,s=lambda.LASSO)[-1,1] # to exclude the intercept
intercept.LASSO <- coef(lasso.obj,s=lambda.LASSO)[1,1]
time.Lasso=proc.time()-ptm
# predicting for the validation set with the lambda found previously
predict.Lasso <- predict(lasso.obj, newx=X_test, s=lambda.LASSO)+intercept.LASSO
## %%%%%%%%%%%%%%%%%% END LASSO  %%%%%%%%%%%%%%%%%% ##   



## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%% RIDGE REGRESSION %%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##  
library (glmnet)
ptm<-proc.time()
rr.obj <- cv.glmnet(X_train, y_train,  
                    family="gaussian",
                    type="mse",  
                    foldid=foldid,
                    alpha=0,     
                    nlambda=100)

lambda.RR    <- rr.obj$lambda.min
beta.Ridge      <- coef(rr.obj,s=lambda.RR)[-1,1] # to exclude the intercept
intercept.RR <- coef(rr.obj,s=lambda.RR)[1,1]

# predicting for the validation set with the lambda found previously
predict.Ridge <- predict(rr.obj, newx=X_test, s=lambda.RR)+intercept.RR
time.Ridge=proc.time()-ptm
## %%%%%%%%%%%%%%%%%% END RIDGE REGRESSION  %%%%%%%%%%%%%%%%%% ## 


# unloading the library
detach("package:glmnet", unload=TRUE)  


## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%% RIDGE REGRESSION BLUP %%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 
# loading the library
library(rrBLUP)
ptm<-proc.time()
# fitting the model  
rrblup.obj <- mixed.solve(y_train,Z=X_train)

# the betas here are the blups
beta.RidgeBlup      <- rrblup.obj$u
intercept.RRblup <- rrblup.obj$beta[1]

# predicting for the validation set
coeffs_x         <- as.matrix(beta.RRblup)
Predict          <- X_test %*% coeffs_x   
predict.RidgeBlup   <- Predict + intercept.RRblup
time.RidgeBlup=-proc.time()-ptm
## %%%%%%%%%%%%%%%%%% END RR-BLUP %%%%%%%%%%%%%%%%%% ##


# unloading the librar
detach("package:rrBLUP", unload=TRUE)  



## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%% ROBUST RIDGE REGRESSION %%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 

####################################################
## Method 1 : Huber Loss, L2 Penalization (hqreg) ##
####################################################
librar (hqreg)
ptm<-proc.time()
enet.obj<-cv.hqreg(X_train,y_train,method="huber",
                   tpe.measure="mse",
                   fold.id = foldid,
                   alpha=0)

lambda.ENET    <- enet.obj$lambda.min
beta.RobRidgeHqreg      <- coef(enet.obj,s=lambda.ENET)[-1] # to exclude the intercept
intercept.ENET <- coef(enet.obj,s=lambda.ENET)[1]

# predicting for the validation set with the lambda found previousl
predict.RobRidgeHqreg <- predict(enet.obj, X=X_test, s=lambda.ENET)+intercept.ENET
time.RobRidgeHqreg=proc.time()-ptm
detach("package:hqreg", unload=TRUE)  

#####################################################
## Method 2 : Scaling Huber Loss , L2 Penalization ##
#####################################################

source("~/Documents/M2 - Math En Action/Stage - Robust/TP Owen CVXR/Launch_Function.R")
ptm<-proc.time()
cv<-cv_Berhu(X_train,y_train,pena="Elastic Net",loss="Huber",M=1.345,alpha=0.1, #It's not working if alpha=0
             FUN="Berhu",nlambda=100,
             fold.id=foldid,nfolds = nb_folds,ncore=5,type.measure="deviance")

lam=cv$lambda.min
beta.RobRidgeScs=cv$beta.min
intercept=cv$mu.min
predict.RobRidgeScs=X_test%*%beta.RobRidgeScs+intercept
time.RobRidgeScs=proc.time()-ptm
#FUN can be Berhu or Adaptive_Berhu for adaptive weights on the penalization
#type.measure can be deviance , mse , mae 
# lam can be lambda.min , lambda.1se


#################################
## Method 3 : Maronna's Method ##
#################################
source("~/Documents/M2 - Math En Action/Stage - Robust/TP MM Regression/R Maronna/Maronna_Changed/Launch_Function_Maronna.R")
ptm<-proc.time()
beta_Maronna<-RobRidge(X_train,y_train,numlam=100 ,cualcv =nb_folds,
                        ncore=7,measure="robust mse",plot.path = F,bestlambda = "1se")$beta
 
beta.RobRidgeMaronna<-beta_Maronna[1:p] beta_Maronna[1:p]
intercept<-beta_Maronna[1:p]
predict.RobRidgeMaronna <- X_test%*%beta.RobRidgeMaronna +intercept
time.RobRidgeMaronna=proc.time()-ptm


#Measure can be robust mse , mae , trimmed mse,mse

## IMPORTANT IT CANT NOT BE THE SAME FOLD (DIMENSION IS REDUCED DURING THE METHOD)

## %%%%%%%%%%%%%%%%%% END ROBUST RIDGE  %%%%%%%%%%%%%%%%%% ##

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%% ROBUST LASSO REGRESSION %%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 

####################################################
## Method 1 : Huber Loss, L1 Penalization (hqreg) ##
####################################################
library (hqreg)
ptm<-proc.time()
enet.obj<-cv.hqreg(X_train,y_train,method="huber",
                   type.measure="mse",
                   fold.id = foldid,
                   alpha=1)

lambda.ENET    <- enet.obj$lambda.min
beta.RobLassoHqreg      <- coef(enet.obj,s=lambda.ENET)[-1] # to exclude the intercept
intercept.ENET <- coef(enet.obj,s=lambda.ENET)[1]

# predicting for the validation set with the lambda found previously
predict.RobLassoHqreg <- predict(enet.obj, X=X_test, s=lambda.ENET)+intercept.ENET
time.RobLassoHqreg=proc.time()-ptm

detach("package:hqreg", unload=TRUE)  
####################################################
## Method 2 : Scaling Huber Loss, L1 Penalization ##
####################################################

source("~/Documents/M2 - Math En Action/Stage - Robust/TP Owen CVXR/Launch_Function.R")
ptm<-proc.time()
cv<-cv_Berhu(X_train,y_train,pena="Elastic Net",loss="Huber",M=1.345,alpha=1,
             FUN="Berhu",nlambda=100,
             fold.id=foldid,nfolds = nb_folds,
             ncore=5,type.measure="deviance")

lam=cv$lambda.min
beta.RobLassoScs=cv$beta.min
intercept=cv$mu.min
predict.RobLassoScs=X_test%*%beta.RobLassoScs+intercept
time.RobLassoScs=proc.time()-ptm
#FUN can be Berhu or Adaptive_Berhu for adaptive weights on the penalization
#type.measure can be deviance , mse , mae 
# lam can be lambda.min , lambda.1se


## %%%%%%%%%%%%%%%%%% END ROBUST LASSO  %%%%%%%%%%%%%%%%%% ##

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%% ROBUST BERHU REGRESSION %%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 

#######################################################
## Method 1 : Scaling Huber loss, Berhu Penalization ##
#######################################################

source("~/Documents/M2 - Math En Action/Stage - Robust/TP Owen CVXR/Launch_Function.R")

ptm<-proc.time()
cv<-cv_Berhu(X_train,y_train,pena="Berhu",loss="Huber",M=1.345,L=1.345,
             FUN="Berhu",nlambda=100 ,
             fold.id=foldid,nfolds = nb_folds,ncore=5,type.measure="deviance")


lam=cv$lambda.min
beta.RobBerhu=cv$beta.min
intercept=cv$mu.min

predict.RobBerhu=X_test%*%beta.RobBerhu+intercept
time.RobBerhu=proc.time()-ptm
#FUN can be Berhu or Adaptive_Berhu for adaptive weights on the penalization
#type.measure can be deviance , mse , mae 
# lam can be lambda.min , lambda.1se

## %%%%%%%%%%%%%%%%%% END ROBUST BERHU  %%%%%%%%%%%%%%%%%% ##



# %%%%%%%%%%%%%%%%%%%%%%   TABLES WITH RESULTS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
methods.names<-c("EN","Lasso","Ridge","RidgeBlup",
                 "RobRidgeHqreg","RobRidgeScs",#"RobRidgeMaronna",
                 "RobLassoHqreg","RobLassoScs","RobBerhu")

# keeping the estimates
estim.beta.coefs<-cbind(beta.EN,beta.Lasso,beta.Ridge,beta.RidgeBlup,
                        beta.RobRidgeHqreg,beta.RobRidgeScs,#beta.RobRidgeMaronna,
                        beta.RobLassoHqreg,beta.RobLassoScs,beta.RobBerhu)
colnames(estim.beta.coefs)<-methods.names

# keeping the predictions
predictions <- cbind(predict.EN,predict.Lasso,predict.Ridge,predict.RidgeBlup,
                    predict.RobRidgeHqreg,predict.RobRidgeScs,#predict.RobRidgeMaronna,
                    predict.RobLassoHqreg,predict.RobLassoScs,predict.RobBerhu)
colnames(predictions)<-methods.names

time<-cbind(time.EN,time.Lasso,time.Ridge,time.RidgeBlup,
            time.RobRidgeHqreg,time.RobRidgeScs,#time.RobRidgeMaronna,
            time.RobLassoHqreg,time.RobLassoScs,time.RobBerhu)
# %%%%%%%%%%%%%%%%%%%%%%   SAVING ALL RESULTS TO AN RDATA FILE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# saving the workspace and thus all the results    
save.image("Results_full_data.Rdata")

