
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

source("~/mcifunction/mcifunction.R")

nsiz<-100000
prodSamp<-function(sampsize){
  x<-rnorm(sampsize)
  a<-rnorm(sampsize)
  b<-a+4*x+rnorm(sampsize)
  c<-2*a-3*b+rnorm(sampsize)
  d<-4*b+5*c+rnorm(sampsize)
  e<-6*d+7*c+rnorm(sampsize)
  f<-8*d+9*e+rnorm(sampsize)
  return(data.frame(cbind(x,a,b,c,d,e,f)))
}

abcdf<-prodSamp(nsiz)[,c("a", "b", "c", "d", "e", "f")]
#pairs(abcde,pch=".")

corm<-cor(abcdf)
suffStat<-list(C = corm, n =nsiz)
rfcf1<-rfci(suffStat=suffStat,indepTest =  gaussCItest, alpha = 0.25, 
           labels =dimnames(corm)[[1]],
           verbose=T)

plot(rfcf1)

xabd<-prodSamp(nsiz)[,c("x","a", "b", "d")]
corm<-cor(xabd)
suffStat<-list(C = corm, n =nsiz)
rfcf2<-rfci(suffStat=suffStat,indepTest =  gaussCItest, alpha = 0.25, 
            labels =dimnames(corm)[[1]],
            verbose=T)

plot(rfcf2)

comb1<-mci(list(rfcf2,rfcf1))
plot(comb1)
# 
# sfStat<-list(data=abcdf[,c("a","b","c","d","f")],ic.method="hsic.gamma")
# rfc2<-rfci(suffStat = sfStat,indepTest = kernelCItest,labels=c("a","b","c","d","f"),
#            alpha=0.2,verbose=TRUE)
# 
# setwd("/home/robin/causal")
