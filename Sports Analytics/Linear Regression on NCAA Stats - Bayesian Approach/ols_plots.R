pdf("fig9_51.pdf",family="Times",height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.5,.5,0))
olsfit<-lm(y~-1+X)
y.te.ols<-X.te%*%olsfit$coef
plot(y.te,y.te.ols,xlab=expression(italic(y)[test]),
     ylab=expression(hat(italic(y))[test])) ; abline(0,1)
mean( (y.te-y.te.ols )^2 )
plot(olsfit$coef,type="h",lwd=2,xlab="regressor index",ylab=expression(hat(beta)[ols]))

dev.off()

## backwards selection
source("backselect.R")

vars<-bselect.tcrit(y,X,tcrit=1.65)
bslfit<-lm(y~-1+X[,vars$remain])
y.te.bsl<-X.te[,vars$remain]%*%bslfit$coef
mean( (y.te-y.te.bsl)^2)
plot(y.te,y.te.bsl,ylim=range( c(y.te.bsl,y.te.ols)),
     xlab=expression(italic(y)[test]),ylab=expression(hat(italic(y))[test]))
abline(0,1)
dev.off()
