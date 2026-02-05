library(readr)

data<-read_csv("insurance.csv")

data<-subset(data,data$sex=="male")
data<-subset(data,data$smoker=="no")
data<-subset(data,data$children==1)
#data<-subset(data,data$region=="southeast")


Xa<-data$bmi
Ca<-data$age
Ya<-data$charges

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
  
  lm<-lm(Y~X+C)

  coeff<-lm$coefficients
  
  tbrc<-()
  
  
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