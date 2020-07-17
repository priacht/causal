
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

udsepTest<-function(x,y,S=NULL,suffStat){
  out<-0
  if(graphNEL2adjMAT(suffStat$g)[x,y]==1){
    gprime<-removeEdge(x,y,suffStat$g)
  }else{if(graphNEL2adjMAT(suffStat$g)[x]==1){
    gprime<-removeEdge(x,y,suffStat$g)
  }
  suffStatprime<-list("g"=gprime,"jp"=suffStat$jp)
  if(dsepTest(x,y,S,suffStatprime)) {
      out<-dsepTest(x,y,S,suffStat)
        
        
    }
  }
  return(out)
}


sampsize<-10000
w<-rnorm(sampsize)
x<-w^2+rnorm(sampsize)/sqrt(3)
y<-w^3-4*w^2+rnorm(sampsize)/sqrt(3)
a<-runif(sampsize)-0.5
b<-(x+y+rnorm(sampsize)/sqrt(3))
c<-runif(sampsize)-0.5
d<-a+0.08*(b)^3-0.3*b+2*abs(c)+runif(sampsize,-1,1)/20
e<-1.0*(c)^2+0.2*d/2+runif(sampsize,0,0.5)


simDat<-data.frame(cbind(w,x,y,a,b,c,d,e))
#pairs(simDat,pch=".")

# sfStat<-list(data=abcde[,c("a","b","d","e")],ic.method="hsic.gamma")
# rfc2<-rfci(suffStat = sfStat,indepTest = kernelCItest,labels=c("a","b","d","e"),
#            alpha=0.2,verbose=TRUE)

# setwd("/home/robin/causal")
# save.image(file=paste("sim1_",Sys.Date(),".RData",sep=""))

simDatcut<-0*simDat
for(i in 1:8){
  brks<-quantile(simDat[,i],probs = seq(0,1,0.2))
  brks[1]<-brks[1]-0.1
  brks[6]<-brks[6]+0.1
  simDatcut[,i]<-cut(simDat[,i],brks,labels=1:5)
  
}




B<-4
par(mfcol=c(2,2))
for(b in 1:B){
  
  simC<-simDatcut[sample(1:nrow(simDatcut),nrow(simDatcut)),]
  bn<-tabu(simC)
  bnNEL<-as(amat(bn),"graphNEL")
  plot(bnNEL,attrs=attrs)
  bnlist[[b]]<-bnNEL
  
  suffStat<-list("g"=bnNEL,"jp"=johnson.all.pairs.sp(bnNEL))
  #  con3<-dag2pag(suffStat,indepTest=dsepTest,alpha=0.5,graph=bnNEL,L=NULL)
  con3<-rfci(suffStat,indepTest=dsepTest,alpha=0.5,labels=nodes(bnNEL))
  plot(con3)
}



mrk<-1
mlist<-list()


B<-50

for(b in 1:B){
  vhere<-sample(names(simDatcut),4)
  sdhere<-simDatcut[,vhere]
  bnhere<-tabu(sdhere)
  bnNEL<-as(amat(bnhere),"graphNEL")
  suffStat<-list("g"=bnNEL,"jp"=johnson.all.pairs.sp(bnNEL))
  #  con3<-dag2pag(suffStat,indepTest=dsepTest,alpha=0.5,graph=bnNEL,L=NULL)
  con3<-rfci(suffStat,indepTest=dsepTest,alpha=0.5,labels=vhere)
  mlist[mrk]<-con3
  mrk<-mrk+1
}

bigmod<-mci(mlist)
plot(bigmod)

for(i in 1:ncol(simDatcut)){
simDatcut[,i]<-as.numeric(simDatcut[,i])-1
}




mrk<-1
mlist<-list()


B<-10

for(b in 1:B){
  vhere<-sample(names(simDatcut),5)
  sdhere<-simDatcut[,vhere]
  
  nlv<-apply(sdhere,2,function(x) {length(unique(x))})
  suffStat<-list(dm=sdhere,nlev=nlv,adaptDF=F)
  con3<-rfci(suffStat,indepTest = disCItest,alpha=0.05,labels=names(sdhere))
  
  mlist[mrk]<-con3
  mrk<-mrk+1
}

bigmod<-mci(mlist)
plot(bigmod)


