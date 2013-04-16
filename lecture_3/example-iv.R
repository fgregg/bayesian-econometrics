############################################################################################
#
#  2-component mixture of bivariate normals:
#
#  Rejection Method
#  Sampling importance resampling
#  Random-walk Metropolis-Hastings algorithm
#  Independent Metropolis-Hastings algorithm
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

post = function(th){pi1*dmvnorm(th,m1,S1)+pi2*dmvnorm(th,m2,S2)+pi3*dmvnorm(th,m3,S3)}

m1  = c(1,4)
m2  = c(4,2)
m3  = c(6.5,2)
S1  = matrix(c(1,-0.9,-0.9,1),2,2)
S2  = matrix(c(1,-0.5,-0.5,1),2,2)
S3  = matrix(c(1,-0.5,-0.5,1),2,2)
pi1 = 1/3
pi2 = 1/3
pi3 = 1-pi1-pi2

N   = 100
th1 = seq(-5,15,length=N)
th2 = seq(-10,15,length=N)
den = matrix(0,N,N)
for (i in 1:N)
  for (j in 1:N)
    den[i,j] = post(c(th1[i],th2[j]))

# Contour plot
par(mfrow=c(1,1))
contour(th1,th2,den,xlab=expression(theta[1]),ylab=expression(theta[2]),drawlabels=FALSE,
xlim=c(-3,10),ylim=c(-1,7.5))

# 3D plot
par(mfrow=c(1,1))
persp(th1,th2,den,theta=30,phi=30,expand=0.5,col="lightblue",
      ltheta=120,shade=0.75,ticktype="detailed",xlab="theta1",
      ylab="theta2",zlab=expression(f(theta)))


# Proposal distribution
mean1 = c(4,2)
Var1  = 9*matrix(c(1,-0.25,-0.25,1),2,2)
den1  = matrix(0,N,N)
for (i in 1:N)
  for (j in 1:N)
    den1[i,j] = dmvnorm(c(th1[i],th2[j]),mean1,Var1)
par(mfrow=c(1,1))
image(th1,th2,den,xlab=expression(theta[1]),ylab=expression(theta[2]),col=topo.colors(12))
contour(th1,th2,den1,add=TRUE,drawlabels=FALSE,lwd=3,col=grey(0.9))

# Accept/reject method
set.seed(1231647)
M     = 10000
draws = rmvnorm(M,mean1,Var1)
prob  = rep(0,M)
for (j in 1:M)
 prob[j]  = post(draws[j,])/(10*dmvnorm(draws[j,],mean1,Var1))
ind = runif(M)<prob
draws1 = draws[ind,]

par(mfrow=c(1,2))
plot(draws,xlab=expression(theta[1]),ylab=expression(theta[2]),pch=16)
contour(th1,th2,den1,col=2,add=TRUE,drawlabels=FALSE)
plot(draws1,xlim=range(draws[,1]),ylim=range(draws[,2]),
     xlab=expression(theta[1]),ylab=expression(theta[2]),pch=16)
contour(th1,th2,den,col=2,add=TRUE,drawlabels=FALSE)

# SIR method
set.seed(1231647)
M      = 10000
M0     = 2000
draws2 = rmvnorm(M,mean1,Var1)
prob   = rep(0,M)
for (j in 1:M)
 prob[j]  = post(draws2[j,])/dmvnorm(draws2[j,],mean1,Var1)
ind = sample(1:M,size=M0,replace=TRUE,prob=prob)
draws3 = draws2[ind,]

par(mfrow=c(1,2))
plot(draws2,xlab=expression(theta[1]),ylab=expression(theta[2]),pch=16)
contour(th1,th2,den1,col=2,add=TRUE,drawlabels=FALSE)
plot(draws3,xlim=range(draws[,1]),ylim=range(draws[,2]),
     xlab=expression(theta[1]),ylab=expression(theta[2]),pch=16)
contour(th1,th2,den,col=2,add=TRUE,drawlabels=FALSE)

# Random walk Metropolis-Hastings
set.seed(16168)
M0    = 1000
M     = 1000
niter = M0+M
Vd    = 0.25*matrix(c(1,-0.5,-0.5,1.0),2,2)
thold = c(8,6)
th.rw = matrix(0,niter,2)
th.rw[1,] = thold
for (iter in 1:niter){
  thnew =  rmvnorm(1,thold,Vd)
  accept = min(1,post(thnew)/post(thold))
  if (runif(1)<accept){
    thold = thnew
  }
  th.rw[iter,] = thold
}
par(mfrow=c(2,2))
plot(th.rw[1:50,],xlim=range(th.rw[,1]),ylim=range(th.rw[,2]),
     xlab=expression(theta[1]),ylab=expression(theta[2]),pch=16,type="s")
contour(th1,th2,den,col=2,add=TRUE,drawlabels=FALSE)
plot(th.rw[51:100,],xlim=range(th.rw[,1]),ylim=range(th.rw[,2]),
     xlab=expression(theta[1]),ylab=expression(theta[2]),pch=16,type="s")
contour(th1,th2,den,col=2,add=TRUE,drawlabels=FALSE)
plot(th.rw[101:150,],xlim=range(th.rw[,1]),ylim=range(th.rw[,2]),
     xlab=expression(theta[1]),ylab=expression(theta[2]),pch=16,type="s")
contour(th1,th2,den,col=2,add=TRUE,drawlabels=FALSE)
plot(th.rw,xlim=range(th.rw[,1]),ylim=range(th.rw[,2]),
     xlab=expression(theta[1]),ylab=expression(theta[2]),pch=16,type="s")
contour(th1,th2,den,col=2,add=TRUE,drawlabels=FALSE)

# Independent Metropolis-Hastings
set.seed(16168)
M0    = 1000
M     = 1000
niter = M0+M
thold = c(8,6)
th.in = matrix(0,niter,2)
th.in[1,] = thold
for (iter in 2:niter){
  thnew =  rmvnorm(1,mean1,Var1)
  accept = min(1,(post(thnew)/post(thold))*(dmvnorm(thold,mean1,Var1)/dmvnorm(thnew,mean1,Var1)))
  if (runif(1)<accept){
    thold = thnew
  }
  th.in[iter,] = thold
}
par(mfrow=c(2,2))
plot(th.in[1:50,],xlim=range(th.in[,1]),ylim=range(th.in[,2]),
     xlab=expression(theta[1]),ylab=expression(theta[2]),pch=16,type="s")
contour(th1,th2,den,col=2,add=TRUE,drawlabels=FALSE)
plot(th.in[51:100,],xlim=range(th.in[,1]),ylim=range(th.in[,2]),
     xlab=expression(theta[1]),ylab=expression(theta[2]),pch=16,type="s")
contour(th1,th2,den,col=2,add=TRUE,drawlabels=FALSE)
plot(th.in[101:150,],xlim=range(th.in[,1]),ylim=range(th.in[,2]),
     xlab=expression(theta[1]),ylab=expression(theta[2]),pch=16,type="s")
contour(th1,th2,den,col=2,add=TRUE,drawlabels=FALSE)
plot(th.in,xlim=range(th.in[,1]),ylim=range(th.in[,2]),
     xlab=expression(theta[1]),ylab=expression(theta[2]),pch=16)
contour(th1,th2,den,col=2,add=TRUE,drawlabels=FALSE)

# Comparing RW Metropolis and independent Metropolis algorithms
par(mfrow=c(2,4))
ts.plot(th.rw[,1],xlab="Iterations",ylab=expression(theta[1]))
acf(th.rw[(M0+1):niter,1],main="")
ts.plot(th.rw[,2],xlab="Iterations",ylab=expression(theta[2]))
acf(th.rw[(M0+1):niter,2],main="")
ts.plot(th.in[,1],xlab="Iterations",ylab=expression(theta[1]))
acf(th.in[(M0+1):niter,1],main="")
ts.plot(th.in[,2],xlab="Iterations",ylab=expression(theta[2]))
acf(th.in[(M0+1):niter,2],main="")
