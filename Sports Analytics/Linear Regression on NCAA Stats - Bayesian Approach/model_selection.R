#### Bayesian model selection
p<-dim(X)[2]
S<-10000
source("regression_gprior.R")

## Don't run it again if you've already run it
#runmcmc<-!any(system("ls",intern=TRUE)=="diabetesBMA.RData")
#if(!runmcmc){ load("diabetesBMA.RData") }

if(runmcmc){
  
  BETA<-Z<-matrix(NA,S,p)
  z<-rep(1,dim(X)[2] )
  lpy.c<-lpy.X(y,X[,z==1,drop=FALSE])
  for(s in 1:S)
  {
    for(j in sample(1:p))
    {
      zp<-z ; zp[j]<-1-zp[j]
      lpy.p<-lpy.X(y,X[,zp==1,drop=FALSE])
      r<- (lpy.p - lpy.c)*(-1)^(zp[j]==0)
      z[j]<-rbinom(1,1,1/(1+exp(-r)))
      if(z[j]==zp[j]) {lpy.c<-lpy.p}
    }
    
    beta<-z
    if(sum(z)>0){beta[z==1]<-lm.gprior(y,X[,z==1,drop=FALSE],S=1)$beta }
    Z[s,]<-z
    BETA[s,]<-beta
    if(s%%10==0)
    { 
      bpm<-apply(BETA[1:s,],2,mean) ; plot(bpm)
      cat(s,mean(z), mean( (y.te-X.te%*%bpm)^2),"\n")
      Zcp<- apply(Z[1:s,,drop=FALSE],2,cumsum)/(1:s)
      plot(c(1,s),range(Zcp),type="n") ; apply(Zcp,2,lines)
    }
  } 
  save(BETA,Z,file="diabetesBMA.RData")
}
#### Figure 9.7
pdf("fig9_7.pdf",family="Times",height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))

beta.bma<-apply(BETA,2,mean,na.rm=TRUE)
y.te.bma<-X.te%*%beta.bma
mean( (y.te-y.te.bma)^2)

layout( matrix(c(1,1,2),nrow=1,ncol=3) )

plot(apply(Z,2,mean,na.rm=TRUE),xlab="regressor index",ylab=expression(
  paste( "Pr(",italic(z[j] == 1),"|",italic(y),",X)",sep="")),type="h",lwd=2)

plot(y.te,y.te.bma,xlab=expression(italic(y)[test]),
     ylab=expression(hat(italic(y))[test])) ; abline(0,1)

dev.off()

#### Figure 9.7
pdf("fig9_7.pdf",family="Times",height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))

beta.bma<-apply(BETA,2,mean,na.rm=TRUE)
y.te.bma<-X.te%*%beta.bma
mean( (y.te-y.te.bma)^2)

layout( matrix(c(1,1,2),nrow=1,ncol=3) )

plot(apply(Z,2,mean,na.rm=TRUE),xlab="regressor index",ylab=expression(
  paste( "Pr(",italic(z[j] == 1),"|",italic(y),",X)",sep="")),type="h",lwd=2)

plot(y.te,y.te.bma,xlab=expression(italic(y)[test]),
     ylab=expression(hat(italic(y))[test])) ; abline(0,1)

dev.off()
## Removing the burning period to get better estimate of posterior of regressors
beta_prime<-BETA[2000:10000,]
beta_pmean<-apply(beta_prime,2,mean,na.rm=TRUE)

## MCMC Diagnostics Autocorrelation between  models
library(coda)
effsize<-effectiveSize(Z)
acf(Z)

## Calculating the squared residue sum
bayes.srr<-mean((y.te-X.te%*%beta_pmean)^2)