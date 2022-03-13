#Compute the economic SVD 
svdecon<-function(x){
  n=dim(x)[1];p=dim(x)[2]
  epsilon=1.e-8
  q=max(p,n)
  
  R=svd(x)
  s=R$d ; a=R$u ; b=R$v
  if(p>n){b=b[,1:n]}
  
  sm=q*s[1]*epsilon
  a=a[,s>sm] ; b=b[,s>sm]
  s=s[s>sm]
  
  result<-list(a=a,s=diag(s),b=b)
  return(result)
  
}
