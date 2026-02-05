library(tikzDevice)

M<-numeric(100)
CAPCE<-numeric(100)
PCAPCE<-numeric(100)
NCAPCE<-numeric(100)

for(m in 1:100){
w0<-1
w1<-(0.1*m)
M[m]<-(0.1*m)

Tf<-function(x){
  pnorm(x,0,w0)-pnorm(x,w1,w1)
}
Pf<-function(x){
max(pnorm(x,0,w0)-pnorm(x,w1,w1),0)
}
Nf<-function(x){
  max(pnorm(x,w1,w1)-pnorm(x,0,w0),0)
}

X<-runif(1000000,-50,50)
CAPCE[m]<-mean(sapply(X,Tf))*40
PCAPCE[m]<-mean(sapply(X,Pf))*100
NCAPCE[m]<-mean(sapply(X,Nf))*100
}

CAPCE
PCAPCE
NCAPCE

plot(M,CAPCE,ylim=c(-1,5),type="l",xlab="",ylab="",lty=1)
par(new=T)
plot(M,PCAPCE,ylim=c(-1,5),type="l",xlab="",ylab="",lty=2)
par(new=T)
plot(M,NCAPCE,ylim=c(-1,5),type="l",xlab="",ylab="",lty=3)



tikz('tikz2.tex',width = 3.5, height = 3)

plot(M,CAPCE,ylim=c(-1,5),type="l",xlab="",ylab="",lty=1)
par(new=T)
plot(M,PCAPCE,ylim=c(-1,5),type="l",xlab="",ylab="",lty=2)
par(new=T)
plot(M,NCAPCE,ylim=c(-1,5),type="l",xlab="",ylab="",lty=3)


dev.off() 
