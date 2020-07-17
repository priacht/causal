


chiZ<-function(lnth){
#  lnth<-length(z)
  w<-rnorm(lnth,0,0.1)
  b<-runif(1,0,2*pi)
  out<-function(z){sqrt(2)*cos(sum(w*z)+b)}
  out
}

xi<-function(d2){
  #output should be a function that takes a d1*d2 matrix
  w<-matrix(rnorm(d2,0,0.1),d2,1)
  b<-runif(1,0,2*pi)
  out<-function(z){sqrt(2)*cos(z%*%w+b)}
  out
}


# x<-(1:100)/10
# plot(NULL,xlim=c(0,10),ylim=c(-2,2))
# for(m in 1:5){
#   rfun<-chiZ(1)
#   for(i in 1:length(x)) y[i]<-rfun(x[i])
#   lines(x,y,col=rgb(0,0,0,0.2))
# }


# rfun<-chiZ(2)
# 
# x1<-(1:100)/10
# x2<-(1:100)/10
# 
# y<-matrix(0,100,100)
# 
# for(i in 1:100){
#   for(j in 1:100){
#     y[i,j]<-rfun(c(x1[i],x2[j]))
#   }
# }
# 
# contour(x1,x2,y)


makeTestDat<-function(nsiz=100000){
  ncola=2
  ncolb=2
  ncolc=2
  nsiz<-100000
  acols<-1:ncola
  bcols<-(ncola+1):(ncola+ncolb)
  ccols<-(ncola+ncolb+1):(ncola+ncolb+ncolc)
  trydat<-matrix(0,nsiz,sum(c(ncola,ncolb,ncolc)))
  trydat[,acols]<-matrix(rnorm(ncola*nsiz,0,10),nsiz)%*%matrix(rnorm(4),2)
  trydat[,ccols]<-sqrt(1+trydat[,acols]^2)+matrix(rnorm(ncolc*nsiz,0,1),nsiz)%*%matrix(rnorm(4),2)
  trydat[,bcols]<-trydat[,ccols]*(1+matrix(rnorm(ncolb*nsiz,0,0.1),nsiz))%*%matrix(rnorm(4),2)
  trydat<-data.frame(trydat)
  names(trydat)<-c("a1","a2","b1","b2","c1","c2")
  return(trydat)
}




makeRBaseDat<-function(trydat,rbasis=xi,ncolv=list(acols,bcols,ccols)){
#  ncolv=list(acols,bcols,ccols)
  acols<-ncolv[[1]]
  bcols<-ncolv[[2]]
  ccols<-ncolv[[3]]
  nsiz=nrow(trydat)
  codtstdat<-matrix(0,nsiz,12)
  codtstdat<-data.frame((codtstdat))
  names(codtstdat)<-list()
  for(i in 1:4){
    rfun<-rbasis(length(acols))
    codtstdat[,i]<-rfun(as.matrix(trydat[,acols]))
    names(codtstdat)[i]<-paste("a",i,sep="")
    
  }
  for(i in 1:4){
    rfun<-rbasis(length(ccols))
    codtstdat[,i+4]<-rfun(as.matrix(trydat[,ccols]))
    names(codtstdat)[i+4]<-paste("c",i,sep="")
  }
  for(i in 1:4){
    rfun<-xi(length(bcols))
    codtstdat[,i+8]<-rfun(as.matrix(trydat[,bcols]))
    names(codtstdat)[i+8]<-paste("b",i,sep="")
  }
  return(codtstdat)
}

trydat<-makeTestDat(100000)

acols<-1:2
bcols<-3:4
ccols<-5:6

codtstdat<-makeRBaseDat(trydat,rbasis=xi,ncolv=list(acols,bcols,ccols))

marg<-apply(codtstdat,2,mean)
codtstdat<-sweep(codtstdat,2,marg,"-")

margsd<-apply(codtstdat,2,sd)
codtstdat<-sweep(codtstdat,2,margsd,"/")

acols<-1:4
ccols<-5:8
bcols<-9:12




#crossp<-data.frame(matrix(0,nsiz,choose(12,2)+12))
makeCrossDat<-function(codtstdat){
  nsiz<-nrow(codtstdat)
  crossp<-data.frame(matrix(0,nsiz,144))
  grpv<-rep(c("a","c","b"),each=4)
  colind<-1
  for(i in 1:ncol(codtstdat)){
    for(j in 1:ncol(codtstdat)){
      crossp[,colind]<-codtstdat[,i]*codtstdat[,j]
      nm<-""
      ip<-i%%4+4*(i%%4==0)
      jp<-j%%4+4*(j%%4==0)
      names(crossp)[colind]<-paste(grpv[i],ip,grpv[j],jp,sep="")
      colind<-colind+1
    }
  }
  return(crossp)
}
crossp<-makeCrossDat(codtstdat)
vsighat<-apply(crossp,2,mean)

pihat<-cov(crossp)
pieig<-eigen(pihat)
pieig$values<-ifelse(pieig$values>1e-20,pieig$values,1e-20)
valroot<-diag(sqrt(pieig$values))

B<-1000
samCondCov<-matrix(0,B,16)
for(b in 1:B){
  vsigstar<-matrix(vsighat+
                     pieig$vectors%*%valroot%*%matrix(rnorm(ncol(crossp))),
                   sqrt(ncol(crossp)))
  #vsigstar<-vsighat+
  #                    pieig$vectors%*%valroot%*%matrix(rnorm(ncol(crossp)))/nsiz
  
#  samCondCov[b,]<-c(solve(solve(vsigstar)[c(acols,bcols),c(acols,bcols)])[1:4,5:8])
  samCondCov[b,]<-c(solve(solve(vsigstar)[c(acols,bcols),c(acols,bcols)])[1:4,5:8])
  
}

mCondCor<-apply(samCondCov,2,mean)
covCondCov<-cov(samCondCov)
# cholCondCov<-chol(covCondCov)
# solve(t(cholCondCov),mCondCor)
eigCovCondCov<-eigen(covCondCov)
normVec<-(eigCovCondCov$vectors)%*%diag(sqrt(1/eigCovCondCov$values))%*%matrix(mCondCor)
sum(normVec^2)
