rm(list=ls())
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulated Data  %%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
library(MASS)
set.seed(123)

n <- 100
p <- 10

rho_V=0.9
SNR=10 #Signal to noise ratio
klev=20
kslo=10
eps=0.2 #Contanimation rate

## Generate problem data
noise <- matrix(rnorm(n))
V=matrix(rep(rho_V,p^2),p,p) ; diag(V)<-1 #The matrix of covariance for X_i
X<-mvrnorm(n=n,mu=rep(0,p),Sigma=V) #Create n multivariate normal distribution
temp<-matrix(rnorm(p)) ; b=temp/norm(temp,"F") #uniform distribution on the sphere (more or less)
beta_true=sqrt(p*SNR)*b

#Add Spasity
sparse <- matrix(1,nrow=p)
r=0 #100% means full sparsity and 0% means no sparsity
zero_index <- sample(p)[1:r*p/100]
sparse[zero_index] <- 0

beta_true=beta_true*sparse
y_true <- X %*% beta_true +noise

##Contamination
m=floor(n*eps)
y<-y_true
for(i in 0:m){
  a=runif(p,-1,1) ; a<-matrix(a-sum(a)/p)
  x_0=as.numeric(klev/(t(a)%*%(solve(V)%*%a))) ; x_0<-x_0*a
  X[i,]<-x_0
  y[i]<-X[i,] %*% beta_true*(1+kslo) + noise[i]
}
X_train=X
y_train=y

X_train =X[0:(2*n/3),] #2/3 of data are for training
X_test=X[(2*n/3):n,] #1/3 are fo testing
y_train=matrix(y[0:(2*n/3)])
y_test=matrix(y_true[(2*n/3):n]) #We test with regard with y_true (not contaminates data)



## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Vessel Data  %%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
set.seed(123)

X<-read.delim("~/Documents/Cours/M2 - Math En Action/Robust-Stats-main/Vessel Data/Vessel_X.txt",header=FALSE,sep=",",dec=".")
X<-matrix(unlist(X),180,301)
y<-read.delim("~/Documents/Cours/M2 - Math En Action/Robust-Stats-main/Vessel Data/Vessel_Y.txt",header=FALSE,sep=",",dec=".")
y<-matrix(unlist(y),180,13)
y<-matrix(y[,1]) # y[,i] for 1<i<13

X_train =X[0:149,]
X_test=X[150:180,]
y_train=matrix(y[0:149,])

##################################
## Cross Validation preparation ##
##################################
n <- length(y_train)
nb_folds=5
foldid = sample(rep(seq(nb_folds), length = n))
ncore=10


## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ELASTIC NET  %%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
# loading the library
library (glmnet)
ptm<-proc.time()
EN.obj <- cv.glmnet(X_train, y_train,  
                    family="gaussian",
                    type="mse",  
                    foldid=foldid,
                    alpha=.5,     
                    nlambda=100)

lambda.EN    <- EN.obj$lambda.min
beta.EN      <- coef(EN.obj,s=lambda.EN)[-1,1] # to exclude the intercept
intercept.EN <- coef(EN.obj,s=lambda.EN)[1,1]

# predicting for the validation set with the lambda found previously
predict.EN <- predict(EN.obj, newx=X_test, s=lambda.EN)+intercept.EN
time.EN=proc.time()-ptm
## %%%%%%%%%%%%%%%%%%%%%%%%% END ELASTIC NET  %%%%%%%%%%%%%%%%%%%%%%%%%%% ##



## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LASSO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
# loading the library
library (glmnet)
ptm<-proc.time()
Lasso.obj <- cv.glmnet(X_train, y_train,  
                    family="gaussian",
                    type="mse",  
                    foldid=foldid,
                    alpha=1,     
                    nlambda=100)

lambda.Lasso    <- Lasso.obj$lambda.min
beta.Lasso      <- coef(Lasso.obj,s=lambda.Lasso)[-1,1] # to exclude the intercept
intercept.Lasso <- coef(Lasso.obj,s=lambda.Lasso)[1,1]

# predicting for the validation set with the lambda found previously
predict.lasso <- predict(Lasso.obj, newx=X_test, s=lambda.Lasso)+intercept.Lasso
time.Lasso=proc.time()-ptm
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
beta.RRblup      <- rrblup.obj$u
intercept.RRblup <- rrblup.obj$beta[1]

# predicting for the validation set
coeffs_x         <- as.matrix(beta.RRblup)
Predict          <- X_test %*% coeffs_x   
predict.RidgeBlup   <- Predict + intercept.RRblup
time.RidgeBlup=-proc.time()-ptm
## %%%%%%%%%%%%%%%%%% END RR-BLUP %%%%%%%%%%%%%%%%%% ##


# unloading the library
detach("package:rrBLUP", unload=TRUE)  



## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%% ROBUST RIDGE REGRESSION %%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ## 

####################################################
## Method 1 : Huber Loss, L2 Penalization (hqreg) ##
####################################################
library(hqreg)
ptm<-proc.time()
enet.obj<-cv.hqreg(X_train,y_train,method="huber",
                   type.measure="mse",
                   fold.id = foldid,
                   alpha=0)

lambda.ENET    <- enet.obj$lambda.min
beta.RobRidgeHqreg      <- coef(enet.obj,s=lambda.ENET)[-1] # to exclude the intercept
intercept.ENET <- coef(enet.obj,s=lambda.ENET)[1]

# predicting for the validation set with the lambda found previousl
predict.RobRidgeHqreg <- predict(enet.obj, X=X_test, s=lambda.ENET)+intercept.ENET
time.RobRidgeHqreg=proc.time()-ptm
detach("package:hqreg", unload=TRUE)  



#################################
## Method 3 : Maronna's Method ##
#################################
source("~/Documents/Cours/M2 - Math En Action/Robust-Stats-main/Maronna Method/Launch_Function_Maronna.R")
ptm<-proc.time()
beta_Maronna<-RobRidge(X_train,y_train,numlam=100 ,cualcv =nb_folds,
                       ncore=ncore,measure="robust mse",plot.path = F,
                       bestlambda = "1se",plot.cv = F)$beta

beta.RobRidgeMaronna<-beta_Maronna[-length(beta_Maronna)]
intercept<-beta_Maronna[length(beta_Maronna)]
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


## %%%%%%%%%%%%%%%%%% END ROBUST LASSO  %%%%%%%%%%%%%%%%%% ##




# %%%%%%%%%%%%%%%%%%%%%%   TABLES WITH RESULTS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
methods.names<-c("EN","Lasso","Ridge","RidgeBlup",
                 "RobRidgeHqreg","RobRidgeMaronna",
                 "RobLassoHqreg")

# keeping the estimates
estim.beta.coefs<-cbind(beta.EN,beta.Lasso,beta.Ridge,beta.RRblup,
                        beta.RobRidgeHqreg,beta.RobRidgeMaronna,
                        beta.RobLassoHqreg)
colnames(estim.beta.coefs)<-methods.names

# keeping the predictions
predictions <- cbind(predict.EN,predict.Lasso,predict.Ridge,predict.RidgeBlup,
                     predict.RobRidgeHqreg,predict.RobRidgeMaronna,
                     predict.RobLassoHqreg)
colnames(predictions)<-methods.names

time<-cbind(time.EN,time.Lasso,time.Ridge,time.RidgeBlup,
            time.RobRidgeHqreg,time.RobRidgeMaronna,
            time.RobLassoHqreg)
# %%%%%%%%%%%%%%%%%%%%%%   SAVING ALL RESULTS TO AN RDATA FILE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# saving the workspace and thus all the results    
save.image("Results_full_data.Rdata")


