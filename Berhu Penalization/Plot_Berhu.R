#Plot the Path or the Cross validation for the Berhu/Adaptive Berhu penalization

plot.path.Berhu<-function(X){
  lambda=X$lambda
  beta=X$beta
  if(X$loss=="Least Square"){
    title="Lambda Path for Least Square Regression Penalization Berhu"
    }else{
    title="Lambda Path for Huber Regression Penalization Berhu"
    }
    
  lambda=X$lambda ; beta_vals=X$beta ; 
  plot(0, 0, typ="n",main = title,
       xlab = expression(log(lambda)), ylab = expression(beta), ylim=range(beta),
       xlim = c(max(log(lambda)), min(log(lambda))))
  matlines(log(lambda), t(beta),lwd=2,type='l',lty=1)

  #eps=max(max(abs(range(beta)))/100,0.1) #Under which epsilon beta_i is considerer null
  eps=10e-5
  nvars=colSums(abs(beta)>eps) #number of null variable
  axis(3,at=log(lambda),labels=nvars,tick=FALSE, line=-0.8)
}


plot.cv.Berhu <- function(x, log.l = TRUE, nvars = TRUE, ...){
  l <- x$fit$lambda
  if (log.l) {
    l <- log(l)
    xlab <- expression(log(lambda))
  } else xlab <- expression(lambda)

  ## Calculate y
  L.cve <- x$cve - x$cvse
  U.cve <- x$cve + x$cvse
  y <- x$cve
  L <- L.cve
  U <- U.cve
  ylab = switch(x$type.measure, deviance = "Deviance", mae = "Mean Absolute Error", mse = "Mean Squared Error")
  
  ylim <- range(c(L, U))
  ind <- ((U-L)/diff(ylim) > 1e-3)
  plot.args = list(x=l, y=y, ylim=ylim, xlab=xlab, ylab=ylab, type="n", xlim=rev(range(l)), las=1)
  new.args = list(...)
  #if (length(new.args)) plot.args[names(new.args)] = new.args
  do.call("plot", plot.args)
  suppressWarnings(arrows(x0=l[ind], x1=l[ind], y0=L[ind], y1=U[ind], code=3, angle=90, col="gray60", length=.03))
  points(l, y, col="red", pch=19, cex=.5)
  lines(l, y, col="red", pch=19, cex=.5)
  if (nvars) {
    x$fit<-structure(x$fit,class="hqreg")
    nv <- predict(x$fit, lambda = x$lambda, type = "nvars")
    
  }
}


##To plot the frontier run the following code in "Lauch_function" with Path the result of a Berhu or Adaptive_Berhu
#par(new=TRUE)
#plot(log(path$lambda), L*path$tau, typ="l",lty=6,lwd=4, ylim=range(path$beta),
#     xlim = c(max(log(path$lambda)), min(log(path$lambda))),xlab = "",ylab="")
#par(new=TRUE)
#plot(log(path$lambda), -L*path$tau, typ="l",lty=6,lwd=4, ylim=range(path$beta),
#     xlim = c(max(log(path$lambda)), min(log(path$lambda))),xlab = "",ylab="")
