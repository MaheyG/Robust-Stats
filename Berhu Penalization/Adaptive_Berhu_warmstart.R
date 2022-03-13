# CVXR_Huber function but with the possibility of parallel compuation

SCS.dims_to_solver_dict <- function(cone_dims) {
  cones <- list(f = as.integer(cone_dims@zero),
                l = as.integer(cone_dims@nonpos),
                q = sapply(cone_dims@soc, as.integer),
                ep = as.integer(cone_dims@exp),
                s = sapply(cone_dims@psd, as.integer))
  return(cones)
}



Adaptive_Berhu_warmstart<-function(X,y,M=1.345,L=1.345,lambda,
                         beta_ini="unpenalized", loss=c("Huber","Least Square"),
                         nlambda=10,lambda_min=0.003,lambda_max,max_iters=1000L,verbose=F){
  
  #Error Checking
  if(M<1)stop("M must be greeter than 1")
  if(L<=0)stop("L must be greeter than 0")
  n=dim(X)[1] ; p=dim(X)[2]
  
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
  
  #Initialization of the weights
  gamma<-1
  if(beta_ini[1]=="unpenalized"){
    weight=1/abs(regression_ls_huber(X,y,loss,M))^gamma
  }else{
    weight=1/abs(beta_ini)
  }
  #Variable to Minimize
  beta<-Variable(p) ; sigma<-Variable(1); tau<-Variable(1) ; mu<-Variable(1) # mu is the intercept
  v<-Variable(n) ;  w<-Variable(p) # v and w are used for the minimization of Huber and Berhu functions
  
  #Loss
  Huber<-function(beta,v,mu, sigma){(quad_over_lin(y-X%*%beta-v-mu,sigma)+2*M*p_norm(v,1)+n*sigma)/(2*n)}
  ls <-function(beta,mu){sum((y - X %*% beta-mu)^2) /(2*n)}
  
  #Penalization
  Berhu<-function(lambda,beta,w,tau){lambda*(sum(tau/weight)+quad_over_lin(sqrt(weight)*w,2*L*tau)-sum(w*weight)+p_norm(weight*beta,1)+L*tau*sum(weight)*0.5)}
  EN<-function(lambda,alpha,beta){lambda*(alpha*p_norm(beta,1)+(1-alpha)*p_norm(beta,2)/2)}
  
  #Constraints
  constr<-list(sigma>=0,tau>=0,w>=abs(beta),w>=L*tau)
  
  beta_sigma_tau <- matrix(0, nrow = p+3, ncol = nlambda)
  
  
  #Initialisation with lambda_max
  lambda_max=lambda[1]
  obj <- Huber(beta,v,mu,sigma) + Berhu(lambda_max,beta,w,tau)
  prob <- Problem(Minimize(obj),constr)
  argu<-get_problem_data(prob,solver="SCS")
  A<-argu$data$A
  b<-argu$data$b
  object=argu$data$c
  cone=SCS.dims_to_solver_dict(argu$data$dims) #a CVXR function easy to transplant
  
  beta_mu_sigma_tau <- matrix(0, nrow = p+3, ncol = nlambda)
  #x<-0 ; y<-0 #Choose the dimension of x and y #CARE value for mu,v,w,tau,sigma
  
  #Solve the Problem  with warm_start 0
  control = scs_control(max_iters = max_iters,verbose = verbose) #warm_start=TRUE x,y
  resu=scs(A = A,b=b,obj=object,cone = cone,control = control)
  
  #Update the warm_start
  x<-resu$x ; y<-resu$y ;s<-resu$s
  
  #Collect the result
  direct_soln <- unpack_results(prob, resu, argu$chain, argu$inverse_data) #TAKE A LOOK
  beta_mu_sigma_tau[1:p,1]=direct_soln$getValue(beta)
  beta_mu_sigma_tau[p+1,1]=direct_soln$getValue(mu)
  beta_mu_sigma_tau[p+2,1]=direct_soln$getValue(sigma)
  beta_mu_sigma_tau[p+3,1]=direct_soln$getValue(tau)
  
  for (i in 1:(nlambda-1)){ #Maybe change for lambda\lambda_min
    #Update with new lambda
    l<-lambda[-1][i] #Without lambda_max
    object[(n+3)]=p*l*(1+L/2)
    object[(n+4)]=l
    object[(n+5):(n+4+p)]=-l
    object[(n+5+p):(n+4+2*p)]=l
    
    #Solve with new warm_start
    control = scs_control(max_iters = max_iters,
                          normalize = T,
                          verbose = verbose)#, 
    #eps = 1e-5, 
    #alpha=1.5) 
    initial=list(x=x,y=y,s=s)
    resu=scs(A = A,b=b,obj=object,cone = cone,control = control,initial=initial)
    
    #Update warm start 
    x<-resu$x ; y<-resu$y ;s=resu$s
    
    #Collect result
    direct_soln <- unpack_results(prob, resu, argu$chain, argu$inverse_data) #TAKE A LOOK
    beta_mu_sigma_tau[1:p,i+1]=direct_soln$getValue(beta)
    beta_mu_sigma_tau[p+1,i+1]=direct_soln$getValue(mu)
    beta_mu_sigma_tau[p+2,i+1]=direct_soln$getValue(sigma)
    beta_mu_sigma_tau[p+3,i+1]=direct_soln$getValue(tau)
  }
  
  beta=beta_mu_sigma_tau[1:p,]
  mu=beta_mu_sigma_tau[p+1,]#Intercept 
  sigma=beta_mu_sigma_tau[p+2,]
  tau=beta_mu_sigma_tau[p+3,]
  result<-list(beta=beta,mu=mu,lambda=lambda, sigma=sigma, tau=tau,M=M,L=L,loss=loss,pena="Berhu",wieght=weight)
}
