library(MASS)

####################
###### Data 1 ######
####################
set.seed(123)

n <- 50#80
p <- 10

rho_V=0.9
SNR=10 #Signal to noise ratio
klev=50
kslo=50 #n=80>=50 n=40> =10
eps=0.1 #Contanimation rate
eps_2=1 #Contamination on noise (noise=eps_2*normal)

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
y_true <- X %*% beta_true +eps_2*noise

##Contamination
m=floor(n*eps)
y<-y_true
for(i in 0:m){
  a=runif(p,-1,1) ; a<-matrix(a-sum(a)/p)
  x_0=as.numeric(klev/(t(a)%*%(solve(V)%*%a))) ; x_0<-x_0*a
  X[i,]<-x_0
  y[i]<-X[i,] %*% beta_true*(1+kslo) + noise[i]
}

X=X[sample(n),]
y=y[sample(n),]

#X_train =X[0:(2*n/3),] #2/3 of data are for training
#X_test=X[(2*n/3):n,] #1/3 are for testing
#y_train=matrix(y[0:(2*n/3)])
#y_test=matrix(y_true[(2*n/3):n]) #We test with regard with y_true (not contaminates data)

X_train=X
y_train=matrix(y)

X_test=X
y_test=matrix(y_true)

##################################
## Cross Validation preparation ##
##################################
n <- length(y_train)
nb_folds=5#floor(2*n)/3
foldid = sample(rep(seq(nb_folds), length = n))
ncore=10


## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%% ACCURACY MEASURE  %%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##

accuracy <- function(y_test,predict,method=c("trimmed mse","mse","mae","cov")){

  method<-match.arg(method)
  n=dim(y_train)[1]
  resid=y_test-predict

  if (method == "mse"){
    acc=mean(resid^2)}

  if(method== "cov"){
    acc=cor(y_test,predict)}

  if(method =="trimmed mse"){
    resid<-abs(resid)
    trimmed_resid=resid[resid<=quantile(resid,prob=c(0.9))] #
    acc=mean(trimmed_resid^2)}

  if(method == "mae"){
    acc=mean(abs(resid))}
  result<-list(Accuracy=acc,method=method)
}


## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ELASTIC NET  %%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##


###################
### I - LS Loss ###
###################
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
predict.EN <- predict(EN.obj, newx=X_test, s=lambda.EN)
time.EN=(proc.time()-ptm)[3]
acc.EN=accuracy(y_test,predict.EN)


#######################
### II - Huber Loss ###
#######################
library(hqreg)
ptm<-proc.time()
EN_huber.obj<-cv.hqreg(X_train,y_train,method="huber",
                   type.measure="mse",
                   fold.id = foldid,
                   alpha=.5)

lambda.EN_huber<- EN_huber.obj$lambda.min
beta.EN_huber<- coef(EN_huber.obj,s=lambda.EN_huber)[-1] # to exclude the intercept
intercept.EN_huber <- coef(EN_huber.obj,s=lambda.EN_huber)[1]

# predicting for the validation set with the lambda found previousl
predict.EN_huber <- predict(EN_huber.obj, X=X_test, s=lambda.EN_huber)
time.EN_huber=(proc.time()-ptm)[3]
acc.EN_huber=accuracy(y_test,predict.EN_huber)



######################
### III - LAD Loss ###
######################
library(hqreg)
ptm<-proc.time()
EN_lad.obj<-cv.hqreg(X_train,y_train,method="quantile",
                       type.measure="mse",
                       fold.id = foldid,
                       alpha=.5)

lambda.EN_lad<- EN_lad.obj$lambda.min
beta.EN_lad<- coef(EN_lad.obj,s=lambda.EN_lad)[-1] # to exclude the intercept
intercept.EN_lad <- coef(EN_lad.obj,s=lambda.EN_lad)[1]

# predicting for the validation set with the lambda found previousl
predict.EN_lad <- predict(EN_lad.obj, X=X_test, s=lambda.EN_lad)
time.EN_lad=(proc.time()-ptm)[3]
acc.EN_lad=accuracy(y_test,predict.EN_lad)



##########################
### IV - Pense Package ###
##########################
library(pense)
library(parallel)
ptm<-proc.time()
cluster <- makeCluster(ncore)
fit <- pense_cv(X_train, y_train, alpha = 0.5, cv_k = nb_folds, cl = cluster)
stopCluster(cluster)
beta.EN_pense<-fit$estimates[[1]]$std_beta
intercept.EN_pense<-fit$estimates[[1]]$std_intercept
predict.EN_pense<- X_test%*%beta.EN_pense +intercept.EN_pense
time.EN_pense=(proc.time()-ptm)[3]
acc.EN_pense=accuracy(y_test,predict.EN_pense)

## %%%%%%%%%%%%%%%%%%%%%%%%% END ELASTIC NET  %%%%%%%%%%%%%%%%%%%%%%%%%%% ##



## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LASSO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##

###################
### I - LS Loss ###
###################
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
predict.Lasso <- predict(Lasso.obj, newx=X_test, s=lambda.Lasso)
acc.Lasso=accuracy(y_test,predict.Lasso)
time.Lasso=(proc.time()-ptm)[3]


#######################
### II - Huber Loss ###
#######################
library(hqreg)
ptm<-proc.time()
Lasso_huber.obj<-cv.hqreg(X_train,y_train,method="huber",
                       type.measure="mse",
                       fold.id = foldid,
                       alpha=1)

lambda.Lasso_huber<- Lasso_huber.obj$lambda.min
beta.Lasso_huber<- coef(Lasso_huber.obj,s=lambda.Lasso_huber)[-1] # to exclude the intercept
intercept.Lasso_huber <- coef(Lasso_huber.obj,s=lambda.Lasso_huber)[1]

# predicting for the validation set with the lambda found previousl
predict.Lasso_huber <- predict(Lasso_huber.obj, X=X_test, s=lambda.Lasso_huber)
time.Lasso_huber=(proc.time()-ptm)[3]
acc.Lasso_huber=accuracy(y_test,predict.Lasso_huber)



######################
### III - LAD Loss ###
######################
library(hqreg)
ptm<-proc.time()
Lasso_lad.obj<-cv.hqreg(X_train,y_train,method="quantile",
                     type.measure="mse",
                     fold.id = foldid,
                     alpha=1)

lambda.Lasso_lad<- Lasso_lad.obj$lambda.min
beta.Lasso_lad<- coef(Lasso_lad.obj,s=lambda.Lasso_lad)[-1] # to exclude the intercept
intercept.Lasso_lad <- coef(Lasso_lad.obj,s=lambda.Lasso_lad)[1]

# predicting for the validation set with the lambda found previousl
predict.Lasso_lad <- predict(Lasso_lad.obj, X=X_test, s=lambda.Lasso_lad)
time.Lasso_lad=(proc.time()-ptm)[3]
acc.Lasso_lad=accuracy(y_test,predict.Lasso_lad)



##########################
### IV - Pense Package ###
##########################
library(pense)
library(parallel)

ptm<-proc.time()
cluster <- makeCluster(ncore)
fit <- pense_cv(X_train, y_train, alpha = 1, cv_k = nb_folds, cl = cluster)
stopCluster(cluster)
beta.Lasso_pense<-fit$estimates[[1]]$std_beta
intercept.Lasso_pense<-fit$estimates[[1]]$std_intercept
predict.Lasso_pense<- X_test%*%beta.Lasso_pense +intercept.Lasso_pense
time.Lasso_pense=(proc.time()-ptm)[3]
acc.Lasso_pense=accuracy(y_test,predict.Lasso_pense)
## %%%%%%%%%%%%%%%%%% END LASSO  %%%%%%%%%%%%%%%%%% ##

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%% RIDGE REGRESSION %%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##

###################
### I - LS Loss ###
###################
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
predict.Ridge <- predict(rr.obj, newx=X_test, s=lambda.RR)
time.Ridge=(proc.time()-ptm)[3]
acc.Ridge=accuracy(y_test,predict.Ridge)
detach("package:glmnet", unload=TRUE)


#################
### II - BLUP ###
#################
library(rrBLUP)
ptm<-proc.time()
# fitting the model
rrblup.obj <- mixed.solve(y_train,Z=X_train)

# the betas here are the blups
beta.Ridge_blup      <- rrblup.obj$u
intercept.Ridge_blup <- rrblup.obj$beta[1]

# predicting for the validation set
coeffs_x         <- as.matrix(beta.Ridge_blup)
Predict          <- X_test %*% coeffs_x
predict.Ridge_blup   <- Predict + intercept.Ridge_blup
time.Ridge_blup=(proc.time()-ptm)[3]
acc.Ridge_blup=accuracy(y_test,predict.Ridge_blup)

# unloading the library
detach("package:rrBLUP", unload=TRUE)

########################
### III - Huber Loss ###
########################
library(hqreg)
ptm<-proc.time()
Ridge_huber.obj<-cv.hqreg(X_train,y_train,method="huber",
                          type.measure="mse",
                          fold.id = foldid,
                          alpha=0)

lambda.Ridge_huber<- Ridge_huber.obj$lambda.min
beta.Ridge_huber<- coef(Ridge_huber.obj,s=lambda.Ridge_huber)[-1] # to exclude the intercept
intercept.Ridge_huber <- coef(Ridge_huber.obj,s=lambda.Ridge_huber)[1]

# predicting for the validation set with the lambda found previousl
predict.Ridge_huber <- predict(Ridge_huber.obj, X=X_test, s=lambda.Ridge_huber)
time.Ridge_huber=(proc.time()-ptm)[3]
acc.Ridge_huber=accuracy(y_test,predict.Ridge_huber)



#####################
### IV - LAD Loss ###
#####################
library(hqreg)
ptm<-proc.time()
Ridge_lad.obj<-cv.hqreg(X_train,y_train,method="quantile",
                     type.measure="mse",
                     fold.id = foldid,
                     alpha=0)

lambda.Ridge_lad<- Ridge_lad.obj$lambda.min
beta.Ridge_lad<- coef(Ridge_lad.obj,s=lambda.Ridge_lad)[-1] # to exclude the intercept
intercept.Ridge_lad <- coef(Ridge_lad.obj,s=lambda.Ridge_lad)[1]

# predicting for the validation set with the lambda found previousl
predict.Ridge_lad <- predict(Ridge_lad.obj, X=X_test, s=lambda.Ridge_lad)
time.Ridge_lad=(proc.time()-ptm)[3]
acc.Ridge_lad=accuracy(y_test,predict.Ridge_lad)
detach("package:hqreg", unload=TRUE)


#########################
### V - Pense Package ###
#########################
library(pense)
library(parallel)

ptm<-proc.time()
cluster <- makeCluster(ncore)
fit <- pense_cv(X_train, y_train, alpha = 0, cv_k = nb_folds, cl = cluster)
stopCluster(cluster)
beta.Ridge_pense<-fit$estimates[[1]]$std_beta
intercept.Ridge_pense<-fit$estimates[[1]]$std_intercept
predict.Ridge_pense<- X_test%*%beta.Ridge_pense +intercept.Ridge_pense
time.Ridge_pense=(proc.time()-ptm)[3]
acc.Ridge_pense=accuracy(y_test,predict.Ridge_pense)
detach("package:pense", unload=TRUE)

#############################
### VI - Maronna's Method ###
#############################

#install.packages("MaronnaRidge_0.1.0.tar.gz",repos=NULL, type="source")

library(MaronnaRidge)
ptm<-proc.time()
beta_Maronna<-RobRidge(X_train,y_train,numlam=100 ,cualcv =nb_folds,
                       ncore=ncore,measure="robust mse",plot.path = F,
                       bestlambda = "1se",plot.cv = F)$beta

beta.Ridge_Maronna<-beta_Maronna[-length(beta_Maronna)]
intercept<-beta_Maronna[length(beta_Maronna)]
predict.Ridge_Maronna <- X_test%*%beta.Ridge_Maronna +intercept
time.Ridge_Maronna=(proc.time()-ptm)[3]
acc.Ridge_Maronna=accuracy(y_test,predict.Ridge_Maronna)
detach("package:MaronnaRidge", unload = TRUE)
## %%%%%%%%%%%%%%%%%% END RIDGE REGRESSION  %%%%%%%%%%%%%%%%%% ##


## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%  BERHU REGRESSION %%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##

###################
### I - LS Loss ###
###################

#install.packages("BerhuPenalization_0.1.0.tar.gz",repos=NULL, type="source")

library(BerhuPenalization)
ptm<-proc.time()
cv<-cv_Berhu(X_train,y_train,pena="Berhu",loss="Least Square",M=1.345,
             FUN="Berhu",nlambda=100,
             fold.id=foldid,nfolds = nb_folds,
             ncore=ncore,type.measure="deviance",max_iters=500L)

#cv<-cv_Berhu(X_train,y_train,pena="Berhu",loss="Least Square",M=1.345,
#             FUN="Adaptive_Berhu",nlambda=100,
#             fold.id=foldid,nfolds = nb_folds,
#             ncore=ncore,type.measure="deviance",max_iters=5000L)

lam=cv$lambda.min
beta.Berhu=cv$beta.min
intercept.Berhu=cv$mu.min
predict.Berhu=X_test%*%beta.Berhu+intercept.Berhu
time.Berhu=(proc.time()-ptm)[3]
acc.Berhu=accuracy(y_test,predict.Berhu)

#FUN can be Berhu or Adaptive_Berhu for adaptive weights on the penalization
#type.measure can be deviance , mse , mae
# lam can be lambda.min , lambda.1se


#######################
### II - Huber Loss ###
#######################


ptm<-proc.time()
cv<-cv_Berhu(X_train,y_train,pena="Berhu",loss="Huber",M=1.345,
             FUN="Berhu",nlambda=100,
             fold.id=foldid,nfolds = nb_folds,
             ncore=ncore,type.measure="deviance",max_iters=5000L)

lam=cv$lambda.min
beta.Berhu_huber=cv$beta.min
intercept.Berhu_huber=cv$mu.min
predict.Berhu_huber=X_test%*%beta.Berhu_huber+intercept.Berhu_huber
time.Berhu_huber=(proc.time()-ptm)[3]
acc.Berhu_huber=accuracy(y_test,predict.Berhu_huber)
detach("package:BerhuPenalization", unload = TRUE)
#FUN can be Berhu or Adaptive_Berhu for adaptive weights on the penalization
#type.measure can be deviance , mse , mae
# lam can be lambda.min , lambda.1se

## %%%%%%%%%%%%%%%%%% END ROBUST BERHU  %%%%%%%%%%%%%%%%%% ##



# %%%%%%%%%%%%%%%%%%%%%%   TABLES WITH RESULTS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
methods.names<-c("EN","EN Huber","EN LAD","EN Pense",
                 "Lasso","Lasso Huber","Lasso LAD","Lasso Pense",
                 "Ridge","RidgeBlup","Ridge Huber","Ridge LAD","Ridge Pense","Ridge Maronna",
                 "Berhu","Berhu Huber")

# keeping the estimates
estim.beta.coefs<-cbind(beta.EN,beta.EN_huber,beta.EN_lad,beta.EN_pense,
                        beta.Lasso,beta.Lasso_huber,beta.Lasso_lad,beta.Lasso_pense,
                        beta.Ridge,beta.Ridge_blup ,beta.Ridge_huber,beta.Ridge_lad,beta.Ridge_pense,beta.Ridge_Maronna,
                        beta.Berhu,beta.Berhu_huber)
colnames(estim.beta.coefs)<-methods.names

# keeping the predictions
predictions <-cbind(predict.EN,predict.EN_huber,predict.EN_lad,predict.EN_pense,
                    predict.Lasso,predict.Lasso_huber,predict.Lasso_lad,predict.Lasso_pense,
                    predict.Ridge,predict.Ridge_blup,predict.Ridge_huber,predict.Ridge_lad,predict.Ridge_pense,predict.Ridge_Maronna,
                    predict.Berhu,predict.Berhu_huber)
colnames(predictions)<-methods.names

time<-cbind(time.EN,time.EN_huber,time.EN_lad,time.EN_pense,
            time.Lasso,time.Lasso_huber,time.Lasso_lad,time.Lasso_pense,
            time.Ridge,time.Ridge_blup,time.Ridge_huber,time.Ridge_lad,time.Ridge_pense,time.Ridge_Maronna,
            time.Berhu,time.Berhu_huber)
colnames(time)<-methods.names

Acc<-cbind(acc.EN,acc.EN_huber,acc.EN_lad,acc.EN_pense,
            acc.Lasso,acc.Lasso_huber,acc.Lasso_lad,acc.Lasso_pense,
            acc.Ridge,acc.Ridge_blup,acc.Ridge_huber,acc.Ridge_lad,acc.Ridge_pense,acc.Ridge_Maronna,
            acc.Berhu,acc.Berhu_huber)
colnames(Acc)<-methods.names
# %%%%%%%%%%%%%%%%%%%%%%   SAVING ALL RESULTS TO AN RDATA FILE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# saving the workspace and thus all the results
#save.image("Results_full_data.Rdata")


