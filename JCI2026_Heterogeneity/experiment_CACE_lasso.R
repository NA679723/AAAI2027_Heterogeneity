

library(readr)
library(tictoc)
library(Rfast)
library(glmnet)
library(MultiRNG)
pi<-3.141592653589
##true value #C=0.5
Na<-1000000
#C<-0.5
UX<-runif(Na,0,1)
UY<-runif(Na,0,1)
D<-100

mu<-rep(1,D)
sigma<-matrix(0,D,D)
for(i in 1:D){
  for(j in 1:D){
    if(i==j){
      sigma[i,j]<-1
    }else{
      sigma[i,j]<-0.9
    }
  }
}
c<-rep(1,D)
beta<-rep(1,D)



c<-rep(0.5,D)

C<-as.numeric(t(beta)%*%(c))
X<-C+UX
Y<-(0.5*X+0.1*C+1)*(-1*UY+0.5)
PI<-runif(Na,0,0.1)
Y0<-(0.5*0+0.1*C+1)*(-1*UY+0.5)
Y2<-(0.5*2+0.1*C+1)*(-1*UY+0.5)
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

#c<-0.5
c<-rep(0.5,D)
tic()



for(boot in 1:BOOT){
  #  for(sto in 1:Sto){
  x0<-0
  x1<-2
  
  
  #C<-runif(N,0,1)
  Cova<-matrix(0,N,D)
  Cova<-draw.d.variate.uniform(N,D,sigma)
  #for(i in 1:D){
  #  for(j in 1:N){
  #  Cova[j,i]<-runif(1,0,1)
  #}}
  #rmvnorm(N,mu,sigma)
  beta<-rep(1,D)
  C<-beta%*%t(Cova)
  UX<-runif(N,0,1)
  UY<-runif(N,0,1)
  X<-2*(C+UX-1.5)/(D-2)
  Y<-(0.5*X+0.1*C+1)*(-1*UY+0.5)
  #Y0<-(0.5*x0+0.1*C+1)*(-1*UY+0.5)
  #Y2<-(0.5*x1+0.1*C+1)*(-1*UY+0.5)
  CACE1<-0
  PCACE1<-0
  NCACE1<-0
  MSE0<-numeric(BOO)
  
  
  ### monte Carlos Integration
  for(boo in 1:BOO){
    y<-runif(1,quantile(Y,0.1),quantile(Y,0.9))
    #y<-runif(1,quantile(Y,0.25),quantile(Y,0.75))
    #y<-rnorm(mean(Y),var(Y))
    
    
    
    ### leaning P(Y<y|X=x0,C=c)
    W0<-matrix(0,N,N)
    for(di in 1:N){
      #W0[di,di]<-kernel(((X[di]-x0)^2+(Cova[di]-c)^2)^(1/2))
      W0[di,di]<-kernel(((X[di]-x0)^2+norm(Cova[di]-c, type="2")^2)^(1/2))
    }
    #XC0<-cbind(1,X-x0,C-c)
    XC0<-cbind(1,t(X-x0),Cova-c*matrix(1,N,D))
    
    #lm0<-.lm.fit(t((as.numeric(Y<y))%*%sqrt(W0)),sqrt(W0)%*%(XC0))
    #lm0<-lm(t((as.numeric(Y<y))%*%sqrt(W0))~sqrt(W0)%*%(XC0))
    #Beta0<-solve(t(XC0)%*%W0%*%XC0)%*%(t(XC0)%*%W0%*%as.numeric(Y<y))
    #alpha0<-Beta0[1]
    #coef(dia.ridge.cv, s="lambda.min")
    #alpha0<-lm0$coefficients[2]
    
    
    #lm0<-lm(t((as.numeric(Y<y))%*%sqrt(W0))~sqrt(W0)%*%(XC0))
    #alpha0<-lm0$coefficients[2]
    
    #dia.ridge.cv <- cv.glmnet(y=as.matrix(t((as.numeric(Y<y))%*%sqrt(W0))), x=as.matrix(sqrt(W0)%*%(XC0)), family="gaussian",standardize = FALSE, alpha=1, intercept=FALSE)
    dia.ridge.cv <- cv.glmnet(y=as.matrix(t((as.numeric(Y<y))%*%sqrt(W0))), x=as.matrix(sqrt(W0)%*%(XC0)), family="gaussian",standardize = FALSE, alpha=0, intercept=FALSE)
    alpha0<-coef(dia.ridge.cv, s="lambda.min")[2]
    
    
    
    #MSE0[boo]<-alpha0-alpha
    ### leaning P(Y<y|X=x1,C=c)
    W1<-matrix(0,N,N)
    for(di in 1:N){
      #W1[di,di]<-kernel(((X[di]-x1)^2+(C[di]-c)^2)^(1/2))
      W1[di,di]<-kernel(((X[di]-x1)^2+norm(Cova[di]-c, type="2")^2)^(1/2))
    }
    #XC1<-cbind(1,X-x1,C-c)
    XC1<-cbind(1,t(X-x1),Cova-c*matrix(1,N,D))
    
    #lm1<-.lm.fit(t((as.numeric(Y<y))%*%sqrt(W1)),sqrt(W1)%*%(XC1))
    #lm1<-lm(t((as.numeric(Y<y))%*%sqrt(W1))~sqrt(W1)%*%(XC1))
    #Beta1<-solve(t(XC1)%*%W1%*%XC1)%*%(t(XC1)%*%W1%*%as.numeric(Y<y))
    #alpha1<-Beta1[1]
    #alpha1<-lm1$coefficients[1]
    #coef(dia.ridge.cv, s="lambda.min")
    #alpha0<-lm0$coefficients[2]
    
    
    #lm1<-lm(t((as.numeric(Y<y))%*%sqrt(W1))~sqrt(W1)%*%(XC1))
    #alpha1<-lm1$coefficients[2]
    
    #dia.ridge.cv <- cv.glmnet(y=as.matrix(t((as.numeric(Y<y))%*%sqrt(W1))), x=as.matrix(sqrt(W1)%*%(XC1)), family="gaussian",standardize = FALSE, alpha=1, intercept=FALSE)
    dia.ridge.cv <- cv.glmnet(y=as.matrix(t((as.numeric(Y<y))%*%sqrt(W1))), x=as.matrix(sqrt(W1)%*%(XC1)), family="gaussian",standardize = FALSE, alpha=0, intercept=FALSE)
    alpha1<-coef(dia.ridge.cv, s="lambda.min")[2]
    
    
    
    
    CACE1<-CACE1+(alpha1-alpha0)*(max(Y)-min(Y))/(BOO)
    PCACE1<-PCACE1+max(alpha1-alpha0,0)*(max(Y)-min(Y))/(BOO)
    NCACE1<-NCACE1+max(alpha0-alpha1,0)*(max(Y)-min(Y))/(BOO)
    #}
  }
  CACE[boot]<-CACE1
  PCACE[boot]<-PCACE1
  NCACE[boot]<-NCACE1
  ##test (choosing h)
  #Ct<-runif(N,0,1)
  #UX<-runif(N,0,1)
  #UY<-runif(N,0,1)
  #Xt<-Ct+UX
  #Yt<-(Xt+Ct+1)*(-1*UY+0.5)
  #MSE<-0
  #for(boo in 1:BOO){
  #  y<-runif(1,min(Y),max(Y))
  #  SK<-0
  #  for(i in 1:length(X)){
  #    SK<-SK+kernel(Xt[i]-x0)*kernel(Ct[i]-c)
  #    SK<-SK+kernel(Xt[i]-x1)*kernel(Ct[i]-c)
  #  }
  # for(i in 1:length(X)){
  #  MSE<-MSE+(as.numeric(Yt[i]<y)-alpha0-beta0*(Xt[i]-x0)-gamma0*(Ct[i]-c))^2*kernel(Xt[i]-x0)*kernel(Ct[i]-c)/SK
  #}
  #for(i in 1:length(X)){
  #  MSE<-MSE+(as.numeric(Yt[i]<y)-alpha1-beta1*(Xt[i]-x1)-gamma1*(Ct[i]-c))^2*kernel(Xt[i]-x1)*kernel(Ct[i]-c)/SK
  #}
  #}
  #  MSEs[boot]<-sum(abs(MSE))
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
PCACE