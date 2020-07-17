
library(bnlearn)
library(RSQLite)
library(gRim)
library(survey)
library(ggplot2)
library(ggm)
library(abind)
library(Rgraphviz)
library(RBGL)
library(gRain)
library(pcalg)
library(kpcalg)

sampsize<-10000
sampsize<-10000
x<-rnorm(sampsize)/sqrt(3)
y<-rnorm(sampsize)/sqrt(3)
a<-runif(sampsize)-0.5
b<-(x+y+rnorm(sampsize)/sqrt(3))
c<-runif(sampsize)-0.5
d<-a+0.08*(b)^3-0.3*b+2*abs(c)+runif(sampsize,-1,1)/20
e<-0.3*(c)^2+0.2*d/2+runif(sampsize,0,0.5)

xybd<-data.frame(cbind(x,y,b,d))

sfStat<-list(data=xybd[,c("x","y","b","d")],ic.method="hsic.gamma")
rfc3<-rfci(suffStat = sfStat,indepTest = kernelCItest,labels=c("x","y","b","d"),
           alpha=0.2,verbose=TRUE)

setwd("/home/robin/causal")
save.image(file=paste("sim_2",Sys.Date(),".RData",sep=""))
