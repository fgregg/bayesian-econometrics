set.seed(7890)
n = 300
beta = 5
sig2 = 0.20
mu   = 0.00
tau2 = 0.10
true = c(beta,sig2,mu,tau2)
sig = sqrt(sig2)
tau = sqrt(tau2)
theta = rnorm(n,mu,tau)
theta.true = theta
x = runif(n)
y = rnorm(n,x*beta+theta,sig)

par(mfrow=c(1,1))
plot(x,y)
points(x,y-theta,col=2,pch=16)

E=sig2;D=0.1;a=(E/D)^2+2;b=(a-1)*E
E=tau2;D=0.1;c=(E/D)^2+2;d=(c-1)*E
b0=0.0
C0=100
mu0=0.0
V0=100


set.seed(3243)
beta  = true[1]
mu    = true[3]
theta = theta.true
M0  = 50000
L   = 50
M   = 5000
niter = M0+L*M
pars = matrix(0,niter,4)
thetas = matrix(0,niter,n)
for (iter in 1:niter){
  # sampling sig2
  sig2  = 1/rgamma(1,a+n/2,b+sum((y-beta*x-theta)^2)/2)
  # sampling tau2
  tau2  = 1/rgamma(1,c+n/2,d+sum((theta-mu)^2)/2)
  # sampling beta
  var   = 1/(1/C0+sum(x^2)/sig2)
  mean  = var*(b0/C0+sum(x*(y-theta))/sig2)
  beta  = rnorm(1,mean,sqrt(var))
  # sampling mu
  var   = 1/(1/V0+n/tau2)
  mean  = var*(mu0/V0+sum(theta)/tau2)
  mu    = rnorm(1,mean,sqrt(var))
  # sampling theta
  var   = 1/(1/tau2+1/sig2)
  mean  = var*(mu/tau2+(y-beta*x)/sig2)
  theta = rnorm(n,mean,sqrt(var))
  # storing draws
  pars[iter,] = c(beta,sig2,mu,tau2)
  thetas[iter,] = theta
}
ind = seq(M0+1,niter,by=L)
pars1 = pars[ind,]
thetas = thetas[ind,]

pdf(file="posterior.pdf",width=10,height=10)
names = c("beta","sig2","mu","tau2")
par(mfrow=c(4,3))
for (i in 1:4){
  ts.plot(pars1[,i],xlab="iterations",ylab="",main=names[i])
  abline(h=true[i],col=2)
  acf(pars1[,i],main="")
  hist(pars1[,i],prob=TRUE,main="",xlab="",breaks=seq(min(pars1[,i]),max(pars1[,i]),length=30))
  abline(v=true[i],col=2)
}
dev.off()

set.seed(12134)
N = 100
xx = seq(min(x),max(x),length=N)
yy = seq(-3,9,length=N)
den = matrix(0,N,N)
theta1 = rnorm(M,pars[,3],sqrt(pars[,4]))
for (i in 1:N)
  for (j in 1:N)
    den[i,j] = mean(dnorm(yy[i],pars1[,1]*xx[j]+theta1,sqrt(pars1[,2])))

set.seed(12164)
draws = matrix(0,N,M)
means = rep(0,N)
for (j in 1:N){
  draws[j,] = rnorm(M,pars1[,1]*xx[j]+theta1,sqrt(pars1[,2]))
  means[j] = mean(pars1[,1]*xx[j]+theta1)
}
meany = apply(draws,1,mean)  

pdf(file="predictives1.pdf",width=10,height=10)
par(mfrow=c(2,1))
persp(yy,xx,den,theta=30,phi=30,expand=0.5,col="lightblue",ltheta=120,shade=0.75,
      ticktype="detailed",xlab="y",ylab="x",zlab="predictive")

contour(xx,yy,t(den),drawlabels=FALSE,xlab="x",ylab="y")
points(x,y,col=gray(0.75),pch=16)
lines(xx,meany,col=2,lwd=5)
lines(xx,means,col=3,lwd=2)
dev.off()

pdf(file="predictives2.pdf",width=10,height=10)
par(mfrow=c(2,2))
for (i in trunc(seq(1,N,length=4))){
	breaks = seq(min(draws[i,]),max(draws[i,]),length=30)
	hist(draws[i,],xlab="y",ylab="Predictive",prob=TRUE,breaks=breaks,main=paste("x=",round(xx[i],3),sep=""))
	lines(yy,den[,i],col=4)
  abline(v=meany[i],col=2,lwd=5)
  abline(v=means[i],col=3,lwd=2)
}
dev.off()




set.seed(12164)
draws = matrix(0,N,M)
means = matrix(0,N,M)
for (j in 1:N){
  draws[j,] = rnorm(M,pars1[,1]*xx[j]+theta1,sqrt(pars1[,2]))
  means[j,] = pars1[,1]*xx[j]+theta1
}


size=100
k=M/size
meanss = array(0,c(k,N,2))
for (i in 1:k){
  ini = (i-1)*size+1
  fini = i*size
  meanss[i,,1]=apply(draws[,ini:fini],1,mean)
  meanss[i,,2]=apply(means[,ini:fini],1,mean)
}

pdf(file="MCerror.pdf",width=10,height=10)
plot(xx,sqrt(apply(meanss[,,1],2,var)),ylim=c(0.02,0.1),type="b",pch=16,xlab="x",
     ylab="MC standard deviation",col=2)
lines(xx,sqrt(apply(meanss[,,2],2,var)),type="b",pch=16,col=3)
dev.off()
