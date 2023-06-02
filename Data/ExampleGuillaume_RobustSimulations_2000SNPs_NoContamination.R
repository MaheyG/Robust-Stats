
# -------------------------------------------------------------------------------- #
#   loading the animal QTLMAS data referring to trait T1  subsetted for 2000 SNPs  #
# -------------------------------------------------------------------------------- #
  source("~/Documents/Cours/M2 - Math En Action/Package Robust/Data/DataPrep2000SNPs.R")

# EXAMPLE OF ONE ROBUST MODEL FIT
# WE WANT TO KEEP PA pred.mspe.hENET pred.mape.hENET trim.mspe.hENET trim.mape.hENET
  # we need to fix the folds for 5-fold CV if we wish to reproduce results
    set.seed(123)
    foldid = sample(rep(seq(5), length = n))

  # setting the percentage of trimming for the trimmed prediciton errors
    n.trim.obs <- length(y_val)*0.1 # number of obs to trim for the trimmed MSPE and MAPE

  # loading the robust library
    library(hqreg)     # for robust regression with huber and LAD losses
                       # although it says one can pass fold.id, the cv.hqreg does not recognize this
                       # and uses the default 10-fold CV even when fold.id is passed
                       # we will therefore set nfolds = 5 in the model fits

# ----------------------------------------------------------------
#    Robust regression with Huber loss penalized by elastic net
# ----------------------------------------------------------------
  startTime.hENET <- Sys.time()
    set.seed(501)
    hENET.obj        <- cv.hqreg(X_train,y_train,method="huber",
                            type.measure="mse",
                            nfolds = 5,
                            alpha=.5,seed = 123)
    pred.hENET       <- predict(hENET.obj,X_test,s=hENET.obj$lambda.min)
    cor.pred.hENET   <- cor(pred.hENET,y_test)
    # computing the MSPE and MAPE
      pred.mspe.hENET <- mean((pred.hENET-y_test)^2)
      pred.mape.hENET <- mean(abs(pred.hENET-y_test))
      trim.mspe.hENET <- mean(sort((pred.hENET-y_test)^2)[1:(n_val-n.trim.obs)])
      trim.mape.hENET <- mean(sort(abs(pred.hENET-y_test))[1:(n_val-n.trim.obs)])
      res.hENET       <- cbind(cor.pred.hENET,pred.mspe.hENET,pred.mape.hENET,trim.mspe.hENET,trim.mape.hENET)
      rm(cor.pred.hENET,pred.mspe.hENET,pred.mape.hENET,trim.mspe.hENET,trim.mape.hENET)
      colnames(res.hENET)[1] <- "PA"
  endTime.hENET  <- Sys.time()
  save.image("ExampleNoCont.RData")
  res.hENET
    #        PA pred.mspe.hENET pred.mape.hENET trim.mspe.hENET trim.mape.hENET
    # 0.6156213        8992.502        76.70195        5657.035        63.58787
  difftime(Sys.time(), startTime.hENET, units = "hours")-difftime(Sys.time(), endTime.hENET, units = "hours")
    # Time difference of 0.007820854 hours


