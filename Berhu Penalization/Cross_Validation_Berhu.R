#Compute Cross Validation for Berhu or Adaptive Berhu penalization regression Huber or Least Square
#nfolds : Number of cross validation folds
#fold.id : A vector of value betwenn 1 and nfold indicating which fold each observation belongs to
#         By defalut the observation are randomly assigned
#type.measure : Mean Square Error, Mean Absolute Error, deviance (uses the choosen loss function of the model)
#seed : Seed for the random number generator in order to obtain reproducible results



cv_Berhu <- function(X, y, ..., FUN = c("Berhu","Adaptive_Berhu"), 
                    ncore = 1, nfolds=10, fold.id,type.measure = c("deviance", "mse", "mae"), seed) {
  
  n <- dim(X)[1] ; p<-dim(X)[2]
  
  
  FUN <- match.arg(FUN) 
  type.measure <- match.arg(type.measure) 
  
  if (!missing(seed)) set.seed(seed)
  if(missing(fold.id)) fold.id <- ceiling(sample(1:n)/n*nfolds) #Same fold.id for each lambda

  #Run one time to have parameters (can be adjusted)
  if(FUN=="Berhu"){
    fit <- Berhu_warmstart(X, y, ...)
    FUN="Berhu_warmstart"}
  if(FUN=="Adaptive_Berhu"){
    fit <- Adaptive_Berhu_warmstart(X, y, ...) 
    FUN="Adaptive_Berhu_warmstart"}
  cv.args <- list(...)
  cv.args$lambda <- fit$lambda

  measure.args <- list(type.measure = type.measure,sigma=fit$sigma, M=fit$M,method=fit$loss)
  E <- matrix(NA, nrow = n, ncol = length(cv.args$lambda))

  # Parallelisation
  parallel <- FALSE
  if (ncore > 1) {
    max.cores <- detectCores()
    if (ncore > max.cores) {
      cat("The number of cores specified (", ncore, ") is larger than 
          the number of avaiable cores (", max.cores, "), so", max.cores, "cores are used.", "\n")
      ncore = max.cores
    }
    cluster <- makeCluster(ncore) #Can add : type="FORK"
    if (!("cluster" %in% class(cluster))) stop("Cluster is not of class 'cluster'; see ?makeCluster")
    parallel <- TRUE
    cat("Start parallel computing for cross-validation...")
    clusterExport(cluster, c("fold.id", "X", "y", "cv.args", "measure.args","FUN","predict.berhu",
                             "measure.berhu","hloss","hloss_sigma","regression_ls_huber","Berhu_warmstart",
                             "Adaptive_Berhu_warmstart","SCS.dims_to_solver_dict","scs_control","scs"),envir = environment())
    clusterEvalQ(cluster,library(CVXR))
    fold.results <- parLapply(cl = cluster, X = 1:nfolds, fun = cvf.exp, XX = X, y = y, 
                              fold.id = fold.id, cv.args = cv.args, measure.args = measure.args,func=FUN)
    stopCluster(cluster)
  }
  
  for (i in 1:nfolds) { 
    if (parallel) {
      fit.i <- fold.results[[i]]
    } else {
      fit.i <- cvf.exp(i, X, y, fold.id, cv.args, measure.args, func=FUN)
      cat("CV fold #",i," in ",nfolds,sep="","\n")
    }
    E[fold.id == i, 1:fit.i$nl] <- fit.i$pe
  }
  
  #Eliminate saturated lambda values (delete infinite value)
  ind <- which(apply(is.finite(E), 2, all)) 
  E <- E[,ind]
  lambda <- cv.args$lambda[ind]
  
  # Results
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, sd) / sqrt(n)
  index.min <- which.min(cve)
  
  # adjust the selection using 1-SD method
  index.1se <- min(which(cve < cve[index.min]+cvse[index.min]))
  val <- list(cve = cve, cvse = cvse, type.measure = type.measure, fit = fit, 
              lambda.1se = lambda[index.1se], lambda.min = lambda[index.min],
              beta.min=fit$beta[,index.min],beta.1se=fit$beta[,index.1se],
              mu.min=fit$mu[index.min],mu.1se=fit$mu[index.1se])
  structure(val, class="cv.hqreg") 
  #class "cv.hqreg" because it comes from package hqreg (to change)
}

cvf.exp <- function(i, XX, y, fold.id, cv.args, measure.args, func) {
  #Training data
  cv.args$X <- XX[fold.id != i,,drop = FALSE] 
  cv.args$y <- y[fold.id != i]
  #Testing data
  X2 <- XX[fold.id == i,,drop = FALSE]
  y2 <- y[fold.id == i]
  fit.i <- do.call(func, cv.args) #Train the data
  yhat <- matrix(predict.berhu(fit.i, X2), length(y2)) #Compare trained data with test data
  list(pe = measure.berhu(y2, yhat, measure.args), nl = length(fit.i$lambda))
}
