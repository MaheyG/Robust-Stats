####Function unitol#############
#Returns the same data as in A (with tolerance eps), but with no repetitions and in sorted order.
#I used the package bazar to construct the function
unitol<- function(x, eps=1.e-10){#Work only for tolerance<1
  n=dim(x)[1]
  if(eps==0){result<-list(y=sort(x),ind=order(x))}
  else{
    y <- floor(x/eps)
    d <- duplicated(y)
    x<-matrix(x[!d])
    
    z<-matrix(1:n) ; z<-matrix(z[!d])
    y=sort(x)  ; ind=z[order(x)]
    result<-list(y=y,ind=ind)
  }
  return(result)
}
