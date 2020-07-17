


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


x<-(1:100)/10
plot(NULL,xlim=c(0,10),ylim=c(-2,2))
for(m in 1:5){
  rfun<-chiZ(1)
  for(i in 1:length(x)) y[i]<-rfun(x[i])
  lines(x,y,col=rgb(0,0,0,0.2))
}


rfun<-chiZ(2)

x1<-(1:100)/10
x2<-(1:100)/10

y<-matrix(0,100,100)

for(i in 1:100){
  for(j in 1:100){
    y[i,j]<-rfun(c(x1[i],x2[j]))
  }
}

contour(x1,x2,y)


nsiz<-10000
trydat<-matrix(0,nsiz,6)
trydat[,1:2]<-matrix(rnorm(2*nsiz,0,1),nsiz)
trydat[,3:4]<-sqrt(1+trydat[,1:2]^2)+matrix(rnorm(2*nsiz,0,1),nsiz)
trydat[,5:6]<-trydat[,3:4]*(1+matrix(rnorm(2*nsiz,0,0.1),nsiz))
rfun<-xi(2)
codtstdat<-matrix(0,nsiz,12)
for(i in 1:4){
  rfun<-xi(4)
  codtstdat[,i]<-rfun(trydat[,1:4])
}
for(i in 1:4){
  rfun<-xi(2)
  codtstdat[,i+4]<-rfun(trydat[,3:4])
}
for(i in 1:4){
  rfun<-xi(2)
  codtstdat[,i+8]<-rfun(trydat[,5:6])
}

marg<-apply(codtstdat,2,mean)
codtstdat<-sweep(codtstdat,2,marg,"-")

acols<-1:4
ccols<-5:8
bcols<-8:12




#crossp<-data.frame(matrix(0,nsiz,choose(12,2)+12))
crossp<-data.frame(matrix(0,nsiz,144))
names(crossp)
grpv<-rep(c("a","b","c"),each=4)
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

crossp<-apply(crossp,2,mean)
crossp<-matrix(crossp,nrow=sqrt(length(crossp)))
crossp<-crossp-


corm<-cor(codtstdat)
sigmaabdotc<-corm[c(1:4,9:12),c(1:4,9:12)]-
  corm[c(1:4,9:12),5:8]%*%solve(corm[5:8,5:8]+0.001*diag(c(1,1,1,1)))%*%corm[5:8,c(1:4,9:12)]

