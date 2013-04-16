############################################################################################
#
#  Gibbs sampler in the bivariate normal distribution
#
############################################################################################
#
# HEDIBERT FREITAS LOPES
# Associate Professor of Econometrics and Statistics
# The University of Chicago Booth School of Business
# 5807 South Woodlawn Avenue
# Chicago, Illinois, 60637
# Email : hlopes@ChicagoBooth.edu
# URL: http://faculty.chicagobooth.edu/hedibert.lopes
#
############################################################################################
rm(list=ls())
m1  = c(0,0)
S1  = matrix(c(1,-0.95,-0.95,1),2,2)
N   = 100
th1 = seq(-5,5,length=N)
th2 = seq(-5,5,length=N)
den = matrix(0,N,N)
for (i in 1:N)
  for (j in 1:N)
    den[i,j] = dmvnorm(c(th1[i],th2[j]),m1,S1)
    
par(mfrow=c(1,1))
contour(th1,th2,den,xlab=expression(theta[1]),ylab=expression(theta[2]),drawlabels=FALSE)

# Gibbs sampler
set.seed(16168)
lag   = 20
M0    = 1000
M     = 1000
niter = M0+lag*M
ths   = matrix(c(-5,5,-5,-5,5,-5,5,5),4,2,byrow=T)
th.g  = array(0,c(4,niter,2))
th    = rep(0,2)
for (i in 1:4){
  th.g[i,1,] = ths[i,]
  for (iter in 2:niter){
    mean  = m1[1]+S1[1,2]/S1[2,2]*(th[2]-m1[2])
    var   = S1[1,1]-S1[1,2]^2/S1[2,2]
    th[1] = rnorm(1,mean,sqrt(var)) 
    mean  = m1[2]+S1[2,1]/S1[1,1]*(th[1]-m1[1])
    var   = S1[2,2]-S1[2,1]^2/S1[1,1]
    th[2] = rnorm(1,mean,sqrt(var)) 
    th.g[i,iter,] = th
  }
}

par(mfrow=c(1,1))
ind = 1:100
i = 1
plot(th.g[i,ind,],xlim=range(th.g),ylim=range(th.g),
     xlab=expression(theta[1]),ylab=expression(theta[2]),pch=16,type="s")
contour(th1,th2,den,col=2,add=TRUE,drawlabels=FALSE)

ind = 1:100
par(mfrow=c(3,4))
for (i in 1:4){
  plot(th.g[i,ind,],xlim=range(th.g),ylim=range(th.g),
       xlab=expression(theta[1]),ylab=expression(theta[2]),pch=16,type="s")
  contour(th1,th2,den,col=2,add=TRUE,drawlabels=FALSE)
}
ind = 1:niter
for (i in 1:4){
  acf(th.g[i,ind,1],main="")
  acf1=acf(th.g[i,ind,2],col=2,plot=FALSE)
  lines(0:(length(acf1$lag)-1)+0.25,acf1$acf,col=2,type="h")
}
ind = seq(M0+1,niter,by=lag)
for (i in 1:4){
  acf(th.g[i,ind,1],main="")
  acf1=acf(th.g[i,ind,2],col=2,plot=FALSE)
  lines(0:(length(acf1$lag)-1)+0.25,acf1$acf,col=2,type="h")
}

# Posterior inference
ind = seq(M0+1,niter,by=lag)
ths = rbind(th.g[1,ind,],
th.g[2,ind,],
th.g[3,ind,],
th.g[4,ind,])

par(mfrow=c(1,3))
plot(ths,xlab=expression(theta[1]),ylab=expression(theta[2]),pch=16)
contour(th1,th2,den,col=2,add=TRUE,drawlabels=FALSE,lwd=2)
xxx = seq(-5,5,length=1000)
hist(ths[,1],xlab="",main=expression(theta[1]),prob=TRUE,breaks=seq(min(xxx),max(xxx),length=40),
     ylim=c(0,0.45))
lines(xxx,dnorm(xxx),col=2)
hist(ths[,2],xlab="",main=expression(theta[2]),prob=TRUE,breaks=seq(min(xxx),max(xxx),length=40),
     ylim=c(0,0.45))
lines(xxx,dnorm(xxx),col=2)
