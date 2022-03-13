#Compute the Regression for Huber or Least square without penalization
regression_ls_huber<-function(X,y,loss,M=1.345){
  
  n=dim(X)[1] ; p=dim(X)[2]
  
  #Error Checking
  if(M<1)stop("M must be greeter than 1")
  
  #Choose the method
  if(loss=="Least Square"){
    beta<-(lsfit(X,y,intercept = TRUE)$coeff)[1:p+1]
  }
  if(loss=="Huber"){
  #Variable to Minimize
  beta<-Variable(p) ; sigma<-Variable(1); v<-Variable(n) ; mu<-Variable(1) 
  
  #Loss
  Huber_Owen<-function(beta,v,mu, sigma){(quad_over_lin(y-X%*%beta-v-mu,sigma)+2*M*p_norm(v,1)+n*sigma)/(2*n)}
  
  #Constrain
  constr<-list(sigma>=0)
  
  obj <- Huber_Owen(beta,v,mu,sigma)
  prob <- Problem(Minimize(obj),constr)  
  result <- solve(prob)
  beta<-c(result$getValue(beta))
  }
  return(beta)
}

