rm(list=ls())

# loading the training dataset 3001x10008
# which has columns individual_ID Sire Dam Sex Generation    Trait1   Trait2   Trait3 m1 m2...m10000
  dat3 <- read.delim(file="~/Documents/Cours/M2 - Math En Action/Package Robust/Data/QTLMAS_training_2012.txt", header=TRUE,sep="\t" )
  attach(dat3)

# getting the SNP matrix 3001x10000 from dat3; because the last line has 8364 obs as NAs out of the 10000 and
# entry=NA the last line is removed and one now has a 3000x10000 matrix
  SNIPS <- data.frame(dat3[,substr(names(dat3),1,1)=="m"])
  x=as.matrix(SNIPS)
  x=x[-3001,]
  rownames(x)<-dat3[-3001,1]

# getting SNP markers with variance!=0; and variance==0
# I have discovered that some SNPs have NAs, in which case na.rm needs to
# be part of the code below;
# resulting matrices are of 3000x9969 and 3000x31
  x2     <-x[,which(apply(x,2,var,na.rm=T)!=0)]
  x2_var0<-x[,which(apply(x,2,var,na.rm=T)==0)]
  x2_var1<-colnames(x2_var0)  # names SNPs with 0 variance

# loading the validation/prediction/test dataset 1021x10008
# header goes like this: individual_ID Sire  Dam Sex Generation Trait1 Trait2 Trait3 m1 m2 ...
  Testdata <- read.delim(file="~/Documents/Cours/M2 - Math En Action/Package Robust/Data/QTLMAS_prediction_2012.txt")
  Testdata <- Testdata[-1021,] # last position is entry 1;
  # does not have trait value and is thus excluded

# extracting the matrix of SNPs for prediction 1020x10000
  SNIPSTEST<- data.frame(Testdata[,substr(names(Testdata),1,1)=="m"])
  rownames(SNIPSTEST)<-Testdata[,1]

# cleaning SNIPTEST from the SNPs with variance = 0; 1020 x 9880
  x2_val   <-SNIPSTEST[,which(apply(SNIPSTEST,2,var,na.rm=T)!=0)]

# Remove markers with zero variance in training data set from the test data set
# we get a matrix of 1020x9969
  Test_nonzero <-SNIPSTEST[ , !(names(SNIPSTEST) %in% x2_var1)]
  n_val<-dim(Test_nonzero)[1]
  print(dim(Test_nonzero))
  p_val<-dim(Test_nonzero)[2]

# variables for number of genotypes/individuals and markers
  n<-dim(x2)[1]
  p<-dim(x2)[2]

# Convert SNP predictor data to data matrix; matrix of 3000x9969
  X=as.matrix(x2[1:n,1:p])

# Convert SNP validation data to data matrix; matrix of 1021x9969
  X_val=as.matrix(Test_nonzero[1:n_val,1:p_val])

# read in response data (size 3001) and convert to a vector 1x3000
  T1=as.matrix(Trait1)
  y <- T1[1:n, ,drop=FALSE]
  y <- as.vector(y)

# read in the validation yield (traits 1, 2 & 3); 1020x4 table
# note that entry 1 is missing from this validation dataset
# testtraitdata starts with entry 3081 as does Testdata from where the
# test SNP data (X_val) is obtained
  TestTraitdata <- read.delim(file="~/Documents/Cours/M2 - Math En Action/Package Robust/Data/ValidationTraitValues.txt")
  T1_val=as.matrix(TestTraitdata[,2])
  y_val <- T1_val[1:n_val, ,drop=FALSE]
  y_val <- as.vector(T1_val)

# removing all the unnecessary things from memory
#  rm(dat3,SNIPS,SNIPSTEST,x2,x2_var0,x2_var1,x2_val,Test_nonzero,Testdata,TestTraitdata,T1,T1_val,x)

# randomly selecting 2000 SNPs out of the 9969 SNPs
  set.seed(123)
  SNPindex <- sample(1:9969,2000)
  X     <- X[,SNPindex]
  X_val <- X_val[,SNPindex]
  p <- p_val <- dim(X_val)[2]

  X_train=X
  y_train=y
  X_test=X_val
  y_test=y_val

  save.image("~/Documents/Cours/M2 - Math En Action/Package Robust/Data/Data2000SNPs.RData")
