# Original functions of Soni, Dewan and Jain (2012), verbatim (triangular kernel).
soni <- new.env()
local(envir = soni, {


dentquant <- function(z,n,h,...)
{ 
  u <- seq(1/(n+1),n/(n+1),1/(n+1))
  t <- seq(1/n,1,1/n)
  k <- matrix(0,nrow=n,ncol=n)
  k1 <- matrix(0,nrow=n,ncol=n)
  Q <- sort(z)
  for(i in 1:n)
  { 
    
    k[,i] <- (1-abs((Q[i]-z)/h))
    k[,i][k[,i] < 0] <- 0    
  }
  s1 <- apply(k,2,sum)
  f1 <- s1/(n*h)
  q2 <- 1/f1
  for(i in 1:n)
  { 
    
    k1[,i] <- (1-abs((t-u[i])/h))/f1
    k1[,i][k1[,i]<0] <- 0  
  }
  s2 <- apply(k1,2,sum)
  q1 <- s2/(n*h)
  return(q1)
}

## jones1


jon1 <- function(z,n,h)
  { 
  kj <- matrix(0,nrow=n,ncol=n) 
     Q <- sort(z)
    for(i in 1:n)
    {
  kj[,i] <- (1-abs((Q[i]-z)/h))
  kj[,i][kj[,i] < 0] <- 0    
}
s1 <- apply(kj,2,sum)
f1 <- s1/(n*h)
qj1 <- 1/f1
return(qj1)
}



## ## jones 2nd estimator
jon2 <- function(z,n,h)
{
f  <- matrix(0,nrow=n,ncol=n-1)
k2 <- matrix(0,nrow=n,ncol=n-1)
k1  <- matrix(0,nrow=n,ncol=n-1)
k  <- matrix(0,nrow=n,ncol=n-1)
q1 <- matrix(0,nrow=n,ncol=n-1)
Q2 <- numeric(length(n-1))
qj2 <- numeric(length(n))
u <- seq(1/(n+1),n/(n+1),1/(n+1))
for(i in 1:n-1)
{
  Q1 <- sort(z)
  Q2[i] <- Q1[i+1]
  k1[,i] <- (1-abs(u-(i/n))/h)/h
  k1[,i][k1[,i] < 0] <- 0 
  k2[,i] <- (1-abs(u-(i+1)/n)/h)/h
  k2[,i][k2[,i] < 0] <- 0 
  k[,i] <- k1[,i]-k2[,i]
  q1 [,i]<- k[,i]*Q2[i]
}
for(i in 1:n)
{   
  qj2[i] <- sum(q1[i,])
}
return(qj2)
}

n= 50
})
