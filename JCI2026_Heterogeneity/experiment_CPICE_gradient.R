library(readr)
library(tictoc)

pi<-3.141592653589

##true value #C=0.5

Na<-100
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

N<-100
Sto<-10

c<-0.5
un<-0.1


tic()
for(boot in 1:BOOT){
for(sto in 1:Sto){
    
  x0<-runif(1,0,un)
  x1<-x0+1.9
  
  C<-runif(N,0,1)
  UX<-runif(N,0,1)
  UY<-runif(N,0,1)
  X<-C+UX
  Y<-(0.5*X+0.1*C+1)*(-1*UY+0.5)
  
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
    
    alpha0<-mean(as.numeric(Y0<y))
    beta0<-0
    gamma0<-0
    alpha1<-mean(as.numeric(Y2<y))
    beta1<-0
    gamma1<-0
    
    ### leaning P(Y<y|X=x0,C=c)
    for(ite in 1:Ite){
      MSE<-0
      for(i in 1:length(X)){
        MSE<-MSE+(as.numeric(Y[i]<y)-alpha0-beta0*(X[i]-x0)-gamma0*(C[i]-c))^2*kernel(X[i]-x0)*kernel(C[i]-c)
      }
      MSEa<-0
      for(i in 1:length(X)){
        MSEa<-MSEa+(as.numeric(Y[i]<y)-(alpha0+delta)-beta0*(X[i]-x0)-gamma0*(C[i]-c))^2*kernel(X[i]-x0)*kernel(C[i]-c)
      }
      MSEb<-0
      for(i in 1:length(X)){
        MSEb<-MSEb+(as.numeric(Y[i]<y)-alpha0-(beta0+delta)*(X[i]-x0)-gamma0*(C[i]-c))^2*kernel(X[i]-x0)*kernel(C[i]-c)
      }
      MSEc<-0
      for(i in 1:length(X)){
        MSEc<-MSEc+(as.numeric(Y[i]<y)-alpha0-beta0*(X[i]-x0)-(gamma0+delta)*(C[i]-c))^2*kernel(X[i]-x0)*kernel(C[i]-c)
      }
      alpha0<-alpha0-dd*(MSEa-MSE)#/delta
      beta0<-beta0-dd*(MSEb-MSE)#/delta
      gamma0<-gamma0-dd*(MSEc-MSE)#/delta
    }
    
    #MSE0[boo]<-alpha0-alpha
    
    ### leaning P(Y<y|X=x1,C=c)
    for(ite in 1:Ite){
      MSE<-0
      for(i in 1:length(X)){
        MSE<-MSE+(as.numeric(Y[i]<y)-alpha1-beta1*(X[i]-x1)-gamma1*(C[i]-c))^2*kernel(X[i]-x1)*kernel(C[i]-c)
      }
      MSEa<-0
      for(i in 1:length(X)){
        MSEa<-MSEa+(as.numeric(Y[i]<y)-(alpha1+delta)-beta1*(X[i]-x1)-gamma1*(C[i]-c))^2*kernel(X[i]-x1)*kernel(C[i]-c)
      }
      MSEb<-0
      for(i in 1:length(X)){
        MSEb<-MSEb+(as.numeric(Y[i]<y)-alpha1-(beta1+delta)*(X[i]-x1)-gamma1*(C[i]-c))^2*kernel(X[i]-x1)*kernel(C[i]-c)
      }
      MSEc<-0
      for(i in 1:length(X)){
        MSEc<-MSEc+(as.numeric(Y[i]<y)-alpha1-beta1*(X[i]-x1)-(gamma1+delta)*(C[i]-c))^2*kernel(X[i]-x1)*kernel(C[i]-c)
      }
      alpha1<-alpha1-dd*(MSEa-MSE)#/delta
      beta1<-beta1-dd*(MSEb-MSE)#/delta
      gamma1<-gamma1-dd*(MSEc-MSE)#/delta
    }
    
    CACE1<-CACE1+(alpha1-alpha0)*(max(Y)-min(Y))*(1/un)/(BOO*Sto)
    PCACE1<-PCACE1+max(alpha1-alpha0,0)*(max(Y)-min(Y))*(1/un)/(BOO*Sto)
    NCACE1<-NCACE1+max(alpha0-alpha1,0)*(max(Y)-min(Y))*(1/un)/(BOO*Sto)
    
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
