#Compute the Berhu Penalization with Least square of Huber Loss


Berhu<-function(X,y,M=1.345,L=1.345,lambda,nlambda=10,lambda_min=0.003,lambda_max,
                loss=c("Huber","Least Square"),pena=c("Berhu","Elastic Net"),alpha=1){
  
  n=dim(X)[1] ; p=dim(X)[2]
  
  #Error Checking
  if(M<1)stop("M must be greeter than 1")
  if(L<=0)stop("L must be greeter than 0")
  if(alpha<0||alpha>1)stop("alpha must be in 0 and 1")
  
  #Choose the method
  loss<-match.arg(loss)
  H<-1
  if(loss=="Least Square"){H<-0}
  pena<-match.arg(pena)
  B<-1
  if(pena=="Elastic Net"){B<-0}
  
  #Compute the Lambda sequence as in the glmnet package
  if(missing(lambda)){
    if(missing(lambda_max)){
      lambda_max=0
      for (i in 1:p){
        if(lambda_max<abs(X[,i]%*%y)){
          lambda_max<-abs(X[,i]%*%y)
          }
      }
      lambda_max<-as.numeric(lambda_max/(n*alpha)) #TODO Change in function of L
    }
   lambda <- exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))
  }  

  #Variable to Minimize
  beta<-Variable(p) ; sigma<-Variable(1); tau<-Variable(1) ; mu<-Variable(1) # mu is the intercept
  v<-Variable(n) ;  w<-Variable(p) # v and w are used for the minimization of Huber and Berhu functions
  
  #Loss
  Huber<-function(beta,v,mu, sigma){(quad_over_lin(y-X%*%beta-v-mu,sigma)+2*M*p_norm(v,1)+n*sigma)/(2*n)}
  ls <-function(beta,mu){sum((y - X %*% beta-mu)^2) /(2*n)}

  #Penalization
  Berhu<-function(lambda,beta,w,tau){lambda*(p*tau+quad_over_lin(w,2*L*tau)-sum(w)+p_norm(beta,1)+L*p*tau/2)}
  EN<-function(lambda,alpha,beta){lambda*(alpha*p_norm(beta,1)+(1-alpha)*p_norm(beta,2)/2)}
  
  #Constraints
  constr<-list(sigma>=0,tau>=0,w>=abs(beta),w>=L*tau) 

  
  beta_sigma_tau <- matrix(0, nrow = p+3, ncol = nlambda)

  # Concomitant Regression using CVXR package
  beta_sigma_tau <- sapply(lambda,
                      function (lambda) {
                        obj <- Huber(beta,v,mu,sigma)*H + ls(beta,mu)*(1-H) + Berhu(lambda,beta,w,tau)*B + EN(lambda,alpha,beta)*(1-B) 
                        prob <- Problem(Minimize(obj),constr)  
                        result <- solve(prob)
                        list<-rbind(result$getValue(beta),result$getValue(mu),result$getValue(sigma),result$getValue(tau))
                        return(list)
                      })

  beta_val=beta_sigma_tau[1:p,]
  mu=beta_sigma_tau[p+1,]#Intercept 
  sigma=beta_sigma_tau[p+2,]
  tau=beta_sigma_tau[p+3,]
  
  result<-list(beta=beta_val,mu=mu,lambda=lambda, sigma=sigma, tau=tau,M=M,L=L,pena=pena,loss=loss)
}

