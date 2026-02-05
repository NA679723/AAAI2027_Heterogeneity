library(readr)
library(tictoc)

pi<-3.141592653589

##true value #C=0.5

Na<-1000
C<-0.5
UX<-runif(Na,0,1)
UY<-runif(Na,0,1)
X<-C+UX
Y<-(0.5*X+0.1*C+1)*(-1*UY+0.5)

PI<-runif(Na,0,0.1)
Y0<-(0.5*PI+0.1*C+1)*(-1*UY+0.5)
Y2<-(0.5*(PI+1.9)+0.1*C+1)*(-1*UY+0.5)
## P-CACE
mean((Y2-Y0)*as.numeric(Y2-Y0>0))
## N-CACE
mean((Y0-Y2)*as.numeric(Y0-Y2>0))



h<-1
kernel<-function(x){
  (1/sqrt(2*pi))*exp(-1*(x^2/(2*h^2)))  
}



#mean(subset(Y,(X==x1&&C==c)))-mean(subset(Y,(X==x0&&C==c)))

dd<-1
Ite<-10
delta<-0.1
BOO<-100
BOOT<-100
CACE<-numeric(BOOT)
PCACE<-numeric(BOOT)
NCACE<-numeric(BOOT)
MSEs<-numeric(BOOT)

N<-10000
Sto<-10

c<-0.5

tic()
for(boot in 1:BOOT){
  for(sto in 1:Sto){
    
    x0<-runif(1,0,0.1)
    x1<-x0+1.9
    
    C<-runif(N,0,1)
    UX<-runif(N,0,1)
    UY<-runif(N,0,1)
    UY2<-runif(N,0,1)
    X<-C+UX
    Y<-(0.5*X+0.1*C+1)*(-1*UY+0.5)+UY2
    
    Y0<-(0.5*x0+0.1*C+1)*(-1*UY+0.5)
    Y2<-(0.5*x1+0.1*C+1)*(-1*UY+0.5)
    
    CACE1<-0
    PCACE1<-0
    NCACE1<-0
    
    MSE0<-numeric(BOO)
    
    ### monte Carlos Integration
    for(boo in 1:BOO){
      y<-runif(1,min(Y),max(Y))
      #y<-runif(1,quantile(Y,0.25),quantile(Y,0.75))
      #y<-rnorm(mean(Y),var(Y))
      
      
      ### leaning P(Y<y|X=x0,C=c)
      
      W0<-matrix(0,N,N)
     for(di in 1:N){
      W0[di,di]<-kernel(((X[di]-x0)^2+(C[di]-c)^2)^(1/2))
     }      
      XC0<-cbind(1,X-x0,C-c)
      Beta0<-solve(t(XC0)%*%W0%*%XC0)%*%(t(XC0)%*%W0%*%as.numeric(Y<y))
      alpha0<-Beta0[1]
      
      #MSE0[boo]<-alpha0-alpha
      
      ### leaning P(Y<y|X=x1,C=c)
      W1<-matrix(0,N,N)
      for(di in 1:N){
        W1[di,di]<-kernel(((X[di]-x1)^2+(C[di]-c)^2)^(1/2))
      }      
      XC1<-cbind(1,X-x1,C-c)
      Beta1<-solve(t(XC1)%*%W1%*%XC1)%*%(t(XC1)%*%W1%*%as.numeric(Y<y))
      alpha1<-Beta1[1]
      
      CACE1<-CACE1+(alpha1-alpha0)*(max(Y)-min(Y))*(1/0.1)/(BOO*Sto)
      #PCACE1<-PCACE1+max(alpha1-alpha0,0)*(max(Y)-min(Y))*(1/0.1)/(BOO*Sto)
      #NCACE1<-NCACE1+max(alpha0-alpha1,0)*(max(Y)-min(Y))*(1/0.1)/(BOO*Sto)
      
      #CACE1<-CACE1+(alpha1-alpha0)*(max(Y)-min(Y))*(1/0.1)/(BOO*Sto)
      PCACE1<-PCACE1+max(1-alpha0,alpha1)*(max(Y)-min(Y))*(1/0.1)/(BOO*Sto)
      NCACE1<-NCACE1+max(1-alpha1,alpha0)*(max(Y)-min(Y))*(1/0.1)/(BOO*Sto)
      
      #CACE1<-CACE1+(alpha1-alpha0)*(quantile(Y,0.75)-quantile(Y,0.25))/BOO
      #PCACE1<-PCACE1+max(alpha1-alpha0,0)*(quantile(Y,0.75)-quantile(Y,0.25))/BOO
      #NCACE1<-NCACE1+min(alpha1-alpha0,0)*(quantile(Y,0.75)-quantile(Y,0.25))/BOO
    }
  }
  
  
  CACE[boot]<-CACE1
  PCACE[boot]<-PCACE1
  NCACE[boot]<-NCACE1
  
  
  ##test (choosing h)
  
  Ct<-runif(N,0,1)
  UX<-runif(N,0,1)
  UY<-runif(N,0,1)
  Xt<-Ct+UX
  Yt<-(Xt+Ct+1)*(-1*UY+0.5)
  
  MSE<-0
  for(boo in 1:BOO){
    y<-runif(1,min(Y),max(Y))
    SK<-0
    for(i in 1:length(X)){
      SK<-SK+kernel(Xt[i]-x0)*kernel(Ct[i]-c)
      SK<-SK+kernel(Xt[i]-x1)*kernel(Ct[i]-c)
    }
    for(i in 1:length(X)){
      MSE<-MSE+(as.numeric(Yt[i]<y)-alpha0-beta0*(Xt[i]-x0)-gamma0*(Ct[i]-c))^2*kernel(Xt[i]-x0)*kernel(Ct[i]-c)/SK
    }
    for(i in 1:length(X)){
      MSE<-MSE+(as.numeric(Yt[i]<y)-alpha1-beta1*(Xt[i]-x1)-gamma1*(Ct[i]-c))^2*kernel(Xt[i]-x1)*kernel(Ct[i]-c)/SK
    }
  }
  MSEs[boot]<-sum(abs(MSE))
  
  
  
}

toc()
summary(CACE)
summary(PCACE)
summary(NCACE)

quantile(CACE,0.025)
mean(CACE)
quantile(CACE,0.975)

quantile(PCACE,0.025)
mean(PCACE)
quantile(PCACE,0.975)

quantile(NCACE,0.025)
mean(NCACE)
quantile(NCACE,0.975)

### choosing h
mean(abs(MSEs))
