library(readr)

data<-read_csv("insurance.csv")

data<-subset(data,data$sex=="male")
data<-subset(data,data$smoker=="no")
data<-subset(data,data$children==1)
#data<-subset(data,data$region=="southeast")


Xa<-data$bmi
Ca<-data$age
Ya<-data$charges
plot(Ya~Xa)


h<-1
kernel<-function(x){
  (1/sqrt(2*pi))*exp(-1*(x^2/(2*h^2)))  
}

##### whole data
alpha<-0
beta<-0
gamma<-0
alpha<-0
beta<-0
gamma<-0

### leaning P(Y<y|X=x0,C=c)
for(ite in 1:Ite){
  MSE<-0
  for(i in 1:length(X)){
    MSE<-MSE+(as.numeric(Y[i]<y)-alpha-beta*(X[i]-x0)-gamma*(C[i]-c))^2*kernel(X[i]-x0)*kernel(C[i]-c)
  }
  MSEa<-0
  for(i in 1:length(X)){
    MSEa<-MSEa+(as.numeric(Y[i]<y)-(alpha+delta)-beta*(X[i]-x0)-gamma*(C[i]-c))^2*kernel(X[i]-x0)*kernel(C[i]-c)
  }
  MSEb<-0
  for(i in 1:length(X)){
    MSEb<-MSEb+(as.numeric(Y[i]<y)-alpha-(beta+delta)*(X[i]-x0)-gamma*(C[i]-c))^2*kernel(X[i]-x0)*kernel(C[i]-c)
  }
  MSEc<-0
  for(i in 1:length(X)){
    MSEc<-MSEc+(as.numeric(Y[i]<y)-alpha-beta*(X[i]-x0)-(gamma0+delta)*(C[i]-c))^2*kernel(X[i]-x0)*kernel(C[i]-c)
  }
  alpha<-alpha-dd*(MSEa-MSE)#/delta
  beta<-beta-dd*(MSEb-MSE)#/delta
  gamma<-gamma-dd*(MSEc-MSE)#/delta
}
################

#h: 10(4*10^21),1(14),0.1

x0<-20
x1<-40
c<-30

#mean(subset(Y,(X==x1&&C==c)))-mean(subset(Y,(X==x0&&C==c)))

dd<-0.1
Ite<-10
delta<-0.1
BOO<-100
BOOT<-100
CACE<-numeric(BOOT)
PCACE<-numeric(BOOT)
NCACE<-numeric(BOOT)
MSEs<-numeric(BOOT)

for(boot in 1:BOOT){
  
  ## subsampling
  Nd<-length(data$age)
  ID_s<-sample(1:Nd, Nd, replace = TRUE)
  
  X<-Xa[ID_s]
  C<-Ca[ID_s]
  Y<-Ya[ID_s]
  
  CACE1<-0
  PCACE1<-0
  NCACE1<-0
  
  MSE0<-numeric(BOO)
  
  ### monte Carlos Integration
  for(boo in 1:BOO){
    y<-runif(1,min(Y),max(Y))
    #y<-runif(1,quantile(Y,0.25),quantile(Y,0.75))
    #y<-rnorm(mean(Y),var(Y))
    
    alpha0<-0
    beta0<-0
    gamma0<-0
    alpha1<-0
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
    
    MSE0[boo]<-alpha0-alpha
    
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
    
    CACE1<-CACE1+(alpha1-alpha0)*(max(Y)-min(Y))/BOO
    PCACE1<-PCACE1+max(alpha1-alpha0,0)*(max(Y)-min(Y))/BOO
    NCACE1<-NCACE1+min(alpha1-alpha0,0)*(max(Y)-min(Y))/BOO
    
    #CACE1<-CACE1+(alpha1-alpha0)*(quantile(Y,0.75)-quantile(Y,0.25))/BOO
    #PCACE1<-PCACE1+max(alpha1-alpha0,0)*(quantile(Y,0.75)-quantile(Y,0.25))/BOO
    #NCACE1<-NCACE1+min(alpha1-alpha0,0)*(quantile(Y,0.75)-quantile(Y,0.25))/BOO
  }
  
  CACE[boot]<-CACE1
  PCACE[boot]<-PCACE1
  NCACE[boot]<-NCACE1
  
  MSEs[boot]<-sum(abs(MSE))
}

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
