############################################################
#
#  y(i) = alpha + beta*x(i) + e(i)  
#
#  where e(i) are iid Student t with nu degrees of freedom
#
#
############################################################
set.seed(13454)
n = 500
v = 3
alpha = 0
beta = 1
sig2 = 1
e = rt(n,v)
x = 5*runif(n)
y = alpha + beta*x + e

par(mfrow=c(1,2))
plot(x,y,main="Notice the size of the extreme events")
hist(y,breaks=seq(min(y),max(y),length=20),prob=TRUE)
lines(sort(y),dnorm(sort(y),mean(y),sqrt(var(y))),col=2,lwd=2)

# Log posterior distribution
logpost = function(beta){
-beta^2/200 -5.5*sum(log(1+(y-alpha-beta*x)^2/(10*sig2)))
}

# Enveloping the posterior density with q=N(1,0.05^3)
betas = seq(0.5,1.3,length=500)
lp = rep(0,500)
q  = rep(0,500)
for (i in 1:500){
  lp[i] = logpost(betas[i])
  q[i]  = dnorm(betas[i],1,0.05)
}

par(mfrow=c(1,1))
plot(betas,exp(lp-max(lp)),type="l",xlab=expression(beta),ylab="")
lines(betas,q/max(q),col=2)


# drawing from the proposal - N(0.94,0.12^2)
betas1 = rnorm(10000,1,0.05)

# Computing SIR weights
w = rep(0,10000)
for (i in 1:10000)
  w[i] = logpost(betas1[i])-dnorm(betas1[i],1,0.05,log=TRUE)
w = exp(w-max(w))
w = w/sum(w)

# Resampling
betas2 = sample(betas1,size=1000,replace=TRUE,prob=w)

par(mfrow=c(1,1))
hist(betas2,xlim=range(betas2),breaks=seq(min(betas2),max(betas2),length=20),prob=TRUE)
lines(betas,20*exp(lp-max(lp)),col=2)

