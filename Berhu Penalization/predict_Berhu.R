#Compute the difference for the trained data and test data
predict.berhu <- function(object, X, lambda, type=c("response","coefficients","nvars"), exact = FALSE, ...) {
  #NEVER USE THE LAMBDA COEFF??
  type=match.arg(type)
  if (missing(X) && type == "response") stop("Need to supply 'X'")
  beta <- object$beta # coef.berhu(object, lambda = lambda, exact = exact)
  num_coef <- dim(beta)[1]
  intercept_added <- (num_coef > length(object$penalty.factor))
  if (type == "coefficients") return(beta)
  if (intercept_added) {
      b0 <-object$mu
      b <- beta
    if (type == "nvars") {
      if (is.matrix(b)) return(apply(b!=0, 2, sum))
      else return(sum(b!=0))
    }     
    if (type == "response") return(sweep(X %*% b, 2, b0, "+"))
  }
  NULL #WHY?
}

coef.berhu <- function(object, lambda, exact = FALSE, ...) {
  if (missing(lambda)) { #NEVER USE THE LAMBDA COEFF
    beta <- object$beta
  } else if (exact) {
    # augment the lambda sequence with the new values, and refit the model
    ls <- object$lambda
    ind <- match(lambda, ls, 0)
    if (any(ind == 0)) {
      ls <- unique(rev(sort(c(lambda,ls))))
      object <- update(object, lambda=ls)
      ind <- match(lambda, ls)
    }
    beta <- object$beta[, ind]
  } else {
    # use linear interpolation to estimate coefficients for supplied lambda
    ls <- object$lambda
    lambda[lambda>max(ls)] <- max(ls)
    lambda[lambda<min(ls)] <- min(ls)
    ind <- approx(ls, seq(ls), lambda)$y
    left <- floor(ind)
    right <- ceiling(ind)
    weight <- ind %% 1
    beta <- (1-weight)*object$beta[,left] + weight*object$beta[,right]
    if (length(lambda) > 1) colnames(beta) <- round(lambda, 4)
  }
  beta
}
