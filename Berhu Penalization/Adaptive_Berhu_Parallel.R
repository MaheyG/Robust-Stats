#Compute the adaptive Berhu Penalization with Least square of Huber Loss
#beta_ini is either "unpenalized" or a vector of weight of R^p

Adaptive_Berhu_Parallel<-function(X,y,M=1.345,L=1.345,lambda,
                              beta_ini="unpenalized", loss=c("Huber","Least Square"),
                              nlambda=10,lambda_min=0.003,lambda_max,ncore=1){
  
  n=dim(X)[1] ; p=dim(X)[2]
  
  #Security Parallel Package
  if (ncore > 1) {
    max.cores <- detectCores()
    if (ncore > max.cores) {
      cat("The number of cores specified (", ncore, ") is larger than 
          the number of avaiable cores (", max.cores, "), so", max.cores, "cores are used.", "\n")
      ncore = max.cores}
  }
  #Error Checking
  if(M<1)stop("M must be greeter than 1")
  if(L<=0)stop("L must be greeter than 0")
  
  #Choose the method
  loss<-match.arg(loss)
  Huber<-1
  if(loss=="Least Square"){
    Huber<-0}
  
  #Compute the Lambda sequence
  if(missing(lambda)){
    if(missing(lambda_max)){
      lambda_max=0
      for (i in 1:p){
        if(lambda_max<abs(X[,i]%*%y)){
          lambda_max<-abs(X[,i]%*%y)}
      }
      lambda_max<-as.numeric(lambda_max/n) #TODO CHANGE in function of L
    }
    lambda <- exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))
  } 
  
  #Computation of the weights
  gamma<-1
  if(beta_ini[1]=="unpenalized"){
    weight=1/abs(regression_ls_huber(X,y,loss,M))^gamma #Compute the unpenalized regression
  }else{
    weight=1/abs(beta_ini)
  }
  
  #Variable to Minimize
  beta<-Variable(p) ; sigma<-Variable(1); tau<-Variable(1) ; mu<-Variable(1) # mu is the intercept
  v<-Variable(n) ;  w<-Variable(p) # v and w are used for the minimization of Huber and Berhu functions 
  
  #Loss
  Huber_Owen<-function(beta,v,mu,sigma){(quad_over_lin(y-X%*%beta-v-mu,sigma)+2*M*p_norm(v,1)+n*sigma)/(2*n)}
  loss_ls <-function(beta){sum((y - X %*% beta)^2) /(2*n)}
  
  #Penalization
  pen_Berhu<-function(lambda,beta,w,tau){lambda*(sum(tau/weight)+quad_over_lin(sqrt(weight)*w,2*L*tau)-sum(w*weight)+p_norm(weight*beta,1)+L*tau*sum(weight)*0.5)}


  #Constraints
  constr<-list(sigma>=0,tau>=0,w>=abs(beta),w>=L*tau) 
  
  beta_sigma_tau <- matrix(0, nrow = p+3, ncol = nlambda)
    
  # Concomitant Regression
  clust <- makeCluster(ncore,type="FORK") #Prepare for parallelism computation "FORK" only on LINUX
  clusterExport(clust, c('X','y') ,envir=environment())
  clusterCall(clust, function() require(CVXR))
  beta_sigma_tau <- parSapply(clust,lambda,
                              function (lambda) {
                                obj <- Huber_Owen(beta,v,mu,sigma)*Huber + loss_ls(beta)*(1-Huber) + pen_Berhu(lambda,beta,w,tau)
                                prob <- Problem(Minimize(obj),constr)  
                                result <- solve(prob,verbose=T)
                                list<-rbind(result$getValue(beta),result$getValue(mu),result$getValue(sigma),result$getValue(tau))
                                return(list)
                              })
  stopCluster(clust)
  
  beta_val=beta_sigma_tau[1:p,]
  mu=beta_sigma_tau[p+1,]#Intercept 
  sigma=beta_sigma_tau[p+2,]
  tau=beta_sigma_tau[p+3,]
  
  result<-list(beta=beta_val,mu=mu,lambda=lambda, sigma=sigma, tau=tau, M=M,L=L,pena="Berhu",loss=loss,weight=weight)
}

