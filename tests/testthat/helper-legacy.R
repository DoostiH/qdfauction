# Original functions from qdf_b_s.R (Doosti, Dewan & Talebian 2025), kept
# verbatim as the reference implementation, except: `n` is assigned before
# use in dqkc, and the stale-index `u[i]` in RLCVST_N/RLCVIP_N is replaced by `u`.
legacy <- new.env()
local(envir = legacy, {
ker=function(x){
  return(dnorm(x))
}

ker1=function(x,h,u){
  g<-ker((x-u)/h)/h
  return(g)
}

Ker1=function(u,h){
  g<-ker(u/h)/h
  return(g)
}



dqk=function(u,X,h){
  l=length(u);
  n=length(X);
  n1=1:n;
  n2=n1-1;
  dqkf=rep(0,l)
  for(i in 1:l){
    dqkf[i]=sum(X*(Ker1(n2/n-u[i],h)-Ker1(n1/n-u[i],h)))
  }
  return(dqkf)
}

BCVfK=function(m,X,Xsd,ni,n){
  h=1/m
  qi=rep(0,n);
  qiSQ=1/(dqk(ni,X,h))^2;
  for(i in 1:n){
    qi[i]=1/dqk(ni[i],X[-i],h)
  }
  return(sum(Xsd*qiSQ)-2*mean(qi))
}

lstar<-function(x,a){
  if(x>=a) value<-log(x)
  else value<-(x/a)-1+log(a)
  return(value)
}


a_n<-function(X){
  #Sample based th
  #eshold for Robust lcv
  s<-sd(X);
  n<-length(X)
  value<-(log(n)/2)^0.5/(s*n)
  return(value)
}

lsv<-function(x,a){
  lstar<-function(sx,a){
    if(sx>=a) value<-log(sx)
    else
      value<-log(a)-1+(sx/a)
    return(value)}
  return(sapply(x,lstar,a=a))}


RLCVK<-function(m,y){
  ###Cross validation for dqp smooth parameter
  h=1/m
  n=length(y)
  y<-sort(y)
  u<-c(1:n)/n
  lsv<-Vectorize(lstar,"x")
  bc<-numeric(0)
  for(i in 1:n) bc[i]<-dqk(u[i],y[-i],h)
  b<-mean(lsv(bc,a_n(y)))
  return(b)}

RLCVK_N<-function(m,y){
  ###Cross validation for indirect smooth Poisson estimator
  h=1/m
  n=length(y)
  y<-sort(y)
  u<-c(1:n)/n
  #lsv<-Vectorize(lstar,"x")
  bc<-numeric(0)
  for(i in 1:n) bc[i]<-dqk(u[i],y[-i],h)
  an=a_n(y)
  a=-mean(lsv(1/bc,an))
  qi=1/(dqk(u,y,h))
  qiSQ=qi^2;
  Xs<-c(0,y);
  Xsd<-Xs[2:(n+1)]-Xs[1:n]
  xs1<-Xsd*qi
  xs2<-Xsd*qiSQ
  b1<-sum(xs1[qi>=an])
  b2<-(1/(2*an))*sum(xs2[qi<an])
  value<-a+b1+b2
  return(value)}

#################################################Kernel correction#################################
dqkc<-function(u,X,h){
  l<-length(u);
  n<-length(X); L=1:n;
  n<-length(X);
  dqkc<-rep(0,l);
  S<-seq(0,1,by=1/n);
  #Wi<-rep(0,n);
  for(i in 1:l){
    Wi=pnorm((S[2:(n+1)]-u[i])/h)-pnorm((S[1:n]-u[i])/h)
    dqkc[i]<-(sum((X[2:n]-X[1:(n-1)])*Ker1(L[1:n-1]/n-u[i],h))-X[n]*Ker1(1-u[i],h)+X[1]*Ker1(-u[i],h)
              +(Ker1(1-u[i],h)-Ker1(-u[i],h))*sum(X*Wi)/sum(Wi))/sum(Wi)
  }
  return(dqkc)
}

BCVfKC=function(m,X,Xsd,ni,n){
  h=1/m
  qi=rep(0,n);
  qiSQ=1/(dqkc(ni,X,h))^2;
  for(i in 1:n){
    qi[i]=1/dqkc(ni[i],X[-i],h)
  }
  return(sum(Xsd*qiSQ)-2*mean(qi))
}

RLCVKC<-function(m,y){
  ###Cross validation for dqp smooth parameter
  h=1/m
  n=length(y)
  y<-sort(y)
  u<-c(1:n)/n
  lsv<-Vectorize(lstar,"x")
  bc<-numeric(0)
  for(i in 1:n) bc[i]<-dqkc(u[i],y[-i],h)
  b<-mean(lsv(bc,a_n(y)))
  return(b)}

RLCVKC_N<-function(m,y){
  ###Cross validation for indirect smooth Poisson estimator
  h=1/m
  n=length(y)
  y<-sort(y)
  u<-c(1:n)/n
  #lsv<-Vectorize(lstar,"x")
  bc<-numeric(0)
  for(i in 1:n) bc[i]<-dqkc(u[i],y[-i],h)
  an=a_n(y)
  a=-mean(lsv(1/bc,an))
  qi=1/(dqkc(u,y,h))
  qiSQ=qi^2;
  Xs<-c(0,y);
  Xsd<-Xs[2:(n+1)]-Xs[1:n]
  xs1<-Xsd*qi
  xs2<-Xsd*qiSQ
  b1<-sum(xs1[qi>=an])
  b2<-(1/(2*an))*sum(xs2[qi<an])
  value<-a+b1+b2
  return(value)}

################################################################
####################Poisson
####################

Wn=function(a,n){
  w=rep(0,n)
  for (i in 1:n){
    w[i]=dpois(i-1,a)
  }
  return(w)
}



EQ<-function(x,y){
  Y<-sort(y);
  l<-length(x);
  n<-length(y);
  eq<-rep(0,l);
  for(i in 1:l){
    if(x[i]==1)eq[i]<-Y[n]
    else if(x[i]==0) eq[i]<-0
    else eq[i]<-Y[floor(n*x[i])+1]
  }
  return(eq)
}



dqp<-function(u,X,h){
  step<-1/h;
  l<-floor(step);
  k<-seq(0,l,by=1)
  n<-length(u);
  dqp<-rep(0,n);
  eq=EQ(k/step,X)
  for(i in 1:n){
    W=Wn(u[i]*step,l+1);
    wpp=c(-W[1],W[1:l]-W[2:(l+1)]);
    Pn=sum(W);
    wp=(wpp+W*W[l+1]/Pn)/Pn
    dqp[i]=step*sum(wp*eq)
  }
  return(dqp)
}

BCVfP=function(m,X,Xsd,ni,n){
  h=1/m
  qi=rep(0,n);
  qiSQ=1/(dqp(ni,X,h))^2;
  for(i in 1:n){
    qi[i]=1/dqp(ni[i],X[-i],h)
  }
  return(sum(Xsd*qiSQ)-2*mean(qi))
}

RLCVP<-function(m,y){
  ###Cross validation for dqp smooth parameter
  h=1/m
  n=length(y)
  y<-sort(y)
  u<-c(1:n)/n
  lsv<-Vectorize(lstar,"x")
  bc<-numeric(0)
  for(i in 1:n) bc[i]<-dqp(u[i],y[-i],h)
  b<-mean(lsv(bc,a_n(y)))
  return(b)}

RLCVP_N<-function(m,y){
  ###Cross validation for indirect smooth Poisson estimator
  h=1/m
  n=length(y)
  y<-sort(y)
  u<-c(1:n)/n
  #lsv<-Vectorize(lstar,"x")
  bc<-numeric(0)
  for(i in 1:n) bc[i]<-dqp(u[i],y[-i],h)
  an=a_n(y)
  a=-mean(lsv(1/bc,an))
  qi=1/(dqp(u,y,h))
  qiSQ=qi^2;
  Xs<-c(0,y);
  Xsd<-Xs[2:(n+1)]-Xs[1:n]
  xs1<-Xsd*qi
  xs2<-Xsd*qiSQ
  b1<-sum(xs1[qi>=an])
  b2<-(1/(2*an))*sum(xs2[qi<an])
  value<-a+b1+b2
  return(value)}
#############################################################
#############Bernstein
#######################################################
dqb<-function(u,X,h){
  step<-1/h;
  K<-floor(step);
  k<-seq(0,K,by=1)
  K1=K-1;
  k1=seq(0,K1,by=1)
  n<-length(u);
  dqb<-rep(0,n);
  eq=EQ(k/step,X);
  deq=eq[2:(K+1)]-eq[1:K];
  for(i in 1:n){
    B=dbinom(k1,K1,u[i]);
    dqb[i]=step*sum(B*deq)
  }
  return(dqb)
}

BCVfB=function(m,X,Xsd,ni,n){
  h=1/m
  qi=rep(0,n);
  qiSQ=1/(dqb(ni,X,h))^2;
  for(i in 1:n){
    qi[i]=1/dqb(ni[i],X[-i],h)
  }
  return(sum(Xsd*qiSQ)-2*mean(qi))
}

RLCVB<-function(m,y){
  ###Cross validation for dqp smooth parameter
  h=1/m
  n=length(y)
  y<-sort(y)
  u<-c(1:n)/n
  lsv<-Vectorize(lstar,"x")
  bc<-numeric(0)
  for(i in 1:n) bc[i]<-dqb(u[i],y[-i],h)
  b<-mean(lsv(bc,a_n(y)))
  return(b)}

RLCVB_N<-function(m,y){
  ###Cross validation for indirect smooth Poisson estimator
  h=1/m
  n=length(y)
  y<-sort(y)
  u<-c(1:n)/n
  #lsv<-Vectorize(lstar,"x")
  bc<-numeric(0)
  for(i in 1:n) bc[i]<-dqb(u[i],y[-i],h)
  an=a_n(y)
  a=-mean(lsv(1/bc,an))
  qi=1/(dqb(u,y,h))
  qiSQ=qi^2;
  Xs<-c(0,y);
  Xsd<-Xs[2:(n+1)]-Xs[1:n]
  xs1<-Xsd*qi
  xs2<-Xsd*qiSQ
  b1<-sum(xs1[qi>=an])
  b2<-(1/(2*an))*sum(xs2[qi<an])
  value<-a+b1+b2
  return(value)}

##################################################################
#########Indirect Kernel
##################################################################
######################bandwidth selection for kernel ###################
idqfj=function(u,X,h){
  l<-length(u);
  n<-length(X);
  dpf<-rep(0,l);
  S<-seq(0,1,by=1/n);
  Wi<-rep(0,n);
  for(i in 1:l){
    dpf[i]<-1/(mean(dnorm((X-EQ(u[i],X))/h))/h)
  }
  return(dpf)
}


BCVfJ=function(m,X,Xsd,ni,n){
  h=1/m
  qi=rep(0,n);
  qiSQ=1/(idqfj(ni,X,h))^2;
  for(i in 1:n){
    qi[i]=1/idqfj(ni[i],X[-i],h)
  }
  return(sum(Xsd*qiSQ)-2*mean(qi))
}

RLCVJ<-function(m,y){
  ###Cross validation for dqp smooth parameter
  h=1/m
  n=length(y)
  y<-sort(y)
  u<-c(1:n)/n
  lsv<-Vectorize(lstar,"x")
  bc<-numeric(0)
  for(i in 1:n) bc[i]<-idqfj(u[i],y[-i],h)
  b<-mean(lsv(bc,a_n(y)))
  return(b)}

RLCVJ_N<-function(m,y){
  ###Cross validation for indirect smooth Poisson estimator
  h=1/m
  n=length(y)
  y<-sort(y)
  u<-c(1:n)/n
  #lsv<-Vectorize(lstar,"x")
  bc<-numeric(0)
  for(i in 1:n) bc[i]<-idqfj(u[i],y[-i],h)
  an=a_n(y)
  a=-mean(lsv(1/bc,an))
  qi=1/(idqfj(u,y,h))
  qiSQ=qi^2;
  Xs<-c(0,y);
  Xsd<-Xs[2:(n+1)]-Xs[1:n]
  xs1<-Xsd*qi
  xs2<-Xsd*qiSQ
  b1<-sum(xs1[qi>=an])
  b2<-(1/(2*an))*sum(xs2[qi<an])
  value<-a+b1+b2
  return(value)}
###############################################################
###########SDJ
###############################################################
idqfST=function(u,X,h,H){
  l<-length(u);
  n<-length(X);
  dpf<-rep(0,l);
  S<-seq(0,1,by=1/n);
  Wi<-rep(0,n);
  fn<-rep(0,n);
  for(i in 1:l){
    for(j in 1:n){
      Wi[j]<-integrate(ker1,lower=S[j],upper=S[j+1],u=u[i],h=H,abs.tol=0.1^100)$value
      fn[j]<-mean(dnorm((X-X[j])/h))/h
    }
    dpf[i]<-sum(Wi/fn)
  }
  return(dpf)
}

BCVfST=function(m,H,X,Xsd,ni,n){
  h=1/m
  qi=rep(0,n);
  qiSQ=1/(idqfST(ni,X,h,H))^2;
  for(i in 1:n){
    qi[i]=1/idqfST(ni[i],X[-i],h,H)
  }
  return(sum(Xsd*qiSQ)-2*mean(qi))
}

RLCVST<-function(m,H,y){
  ###Cross validation for dqp smooth parameter
  h=1/m
  n=length(y)
  y<-sort(y)
  u<-c(1:n)/n
  lsv<-Vectorize(lstar,"x")
  bc<-numeric(0)
  for(i in 1:n) bc[i]<-idqfST(u[i],y[-i],h,H)
  b<-mean(lsv(bc,a_n(y)))
  return(b)}




RLCVST_N<-function(m,H,y){
  ###Cross validation for indirect smooth Poisson estimator
  h=1/m
  n=length(y)
  y<-sort(y)
  u<-c(1:n)/n
  #lsv<-Vectorize(lstar,"x")
  bc<-numeric(0)
  for(i in 1:n) bc[i]<-idqfST(u[i],y[-i],h,H)
  an=a_n(y)
  a=-mean(lsv(1/bc,an))
  qi=1/(idqfST(u,y,h,H))
  qiSQ=qi^2;
  Xs<-c(0,y);
  Xsd<-Xs[2:(n+1)]-Xs[1:n]
  xs1<-Xsd*qi
  xs2<-Xsd*qiSQ
  b1<-sum(xs1[qi>=an])
  b2<-(1/(2*an))*sum(xs2[qi<an])
  value<-a+b1+b2
  return(value)}



#######################################################
#######Indirect Poisson Quantile 
#######################################################
SPDF<-function(x,X,H,u){
  step<-1/H;
  l<-length(x);
  smf<-rep(0,l);
  M=ceiling(max(X)*step);
  k1=0:M;
  G=stepfuncS(k1/step,X);
  for(i in 1:l){
    smf[i]<-1-sum(Wn(x[i]*step,M+1)*G)-u;
  }
  return(smf)
}   

IvSPDF<-function(u,X,H){
  l<-length(u);
  smfq<-rep(0,l)
  for(i in 1:l){
    smfq[i]=uniroot(SPDF,c(0,max(X)),extendInt="upX",X=X,H=H,u=u[i])$root
  }
  return(smfq)
}


stepfuncS=function(x,y){
  n=length(x)
  u=rep(0,n)
  for (i in 1:n){
    z<-y>x[i]
    u[i]=mean(z)
  }
  return(u)
}

stepfunc=function(x,y){
  n=length(x)
  u=rep(0,n)
  for (i in 1:n){
    z<-y<=x[i]
    u[i]=mean(z)
  }
  return(u)
}

idqfp=function(u,X,h){
  step=1/h;
  l=length(u);
  M=ceiling(max(X)*step)
  k1=c(1:M)
  k2=k1-1
  l1<-floor(step);
  k<-seq(0,l1,by=1)
  G1=stepfunc(k1*h,X)
  G2=stepfunc(k2*h,X)
  w=G1-G2
  dpf=rep(0,l)
  for (i in 1:l){
    dpf[i]=sum(Wn(IvSPDF(u[i],X,h)*step,M)*w)*step
  }
  return(1/dpf)
}

BCVfIP=function(m,X,Xsd,ni,n){
  h=1/m
  qi=rep(0,n);
  qiSQ=1/(idqfp(ni,X,h))^2;
  for(i in 1:n){
    qi[i]=1/idqfp(ni[i],X[-i],h)
  }
  return(sum(Xsd*qiSQ)-2*mean(qi))
}

idqfpT=function(u,X,h,H){
  step=1/h;
  l=length(u);
  M=ceiling(max(X)*step)
  k1=c(1:M)
  k2=k1-1
  l1<-floor(step);
  k<-seq(0,l1,by=1)
  G1=stepfunc(k1*h,X)
  G2=stepfunc(k2*h,X)
  w=G1-G2
  dpf=rep(0,l)
  for (i in 1:l){
    dpf[i]=sum(Wn(IvSPDF(u[i],X,H)*step,M)*w)*step
  }
  return(1/dpf)
}




BCVfIPT=function(h,H,X,Xsd,ni,n){
  qi=rep(0,n);
  qiSQ=1/(idqfpT(ni,X,h,H))^2;
  for(i in 1:n){
    qi[i]=1/idqfpT(ni[i],X[-i],h,H)
  }
  return(sum(Xsd*qiSQ)-2*mean(qi))
}

RLCVIP<-function(m,y){
  ###Cross validation for dqp smooth parameter
  h=1/m
  n=length(y)
  y<-sort(y)
  u<-c(1:n)/n
  lsv<-Vectorize(lstar,"x")
  bc<-numeric(0)
  for(i in 1:n) bc[i]<-idqfp(u[i],y[-i],h)
  b<-mean(lsv(bc,a_n(y)))
  return(b)}

RLCVIP_N<-function(m,y){
  ###Cross validation for indirect smooth Poisson estimator
  h=1/m
  n=length(y)
  y<-sort(y)
  u<-c(1:n)/n
  #lsv<-Vectorize(lstar,"x")
  bc<-numeric(0)
  for(i in 1:n) bc[i]<-idqfp(u[i],y[-i],h)
  an=a_n(y)
  a=-mean(lsv(1/bc,an))
  qi=1/(idqfp(u,y,h))
  qiSQ=qi^2;
  Xs<-c(0,y);
  Xsd<-Xs[2:(n+1)]-Xs[1:n]
  xs1<-Xsd*qi
  xs2<-Xsd*qiSQ
  b1<-sum(xs1[qi>=an])
  b2<-(1/(2*an))*sum(xs2[qi<an])
  value<-a+b1+b2
  return(value)}


RLCVIPT_N<-function(m,H,y){
  ###Cross validation for indirect smooth Poisson estimator
  h=1/m
  n=length(y)
  y<-sort(y)
  u<-c(1:n)/n
  #lsv<-Vectorize(lstar,"x")
  bc<-numeric(0)
  for(i in 1:n) bc[i]<-idqfpT(ni[i],X[-i],h,H)
  an=a_n(y)
  a=-mean(lsv(1/bc,an))
  qi=1/(idqfpT(ni[i],X[-i],h,H))
  qiSQ=qi^2;
  Xs<-c(0,y);
  Xsd<-Xs[2:(n+1)]-Xs[1:n]
  xs1<-Xsd*qi
  xs2<-Xsd*qiSQ
  b1<-sum(xs1[qi>=an])
  b2<-(1/(2*an))*sum(xs2[qi<an])
  value<-a+b1+b2
  return(value)}



  
  
tdqf=function(u,a,b){
  l=length(u);
  dqf=rep(0,l);
  for(i in 1:l){
    dqf[i]=1/dgamma(qgamma(u[i],a,b),a,b);
  }
  return(dqf)
} 
})
