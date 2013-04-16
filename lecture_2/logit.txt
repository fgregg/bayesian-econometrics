prior = function(b){dnorm(b,b0,Sb)}
like  = function(b){prod(exp((a+b*x)*y)/(1+exp(a+b*x)))}
fun   = function(b){prior(b)*like(b)}
post  = function(b){prior(b)*like(b)/py}

# O-ring data
# -----------
y  = c(1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0)
t  = c(53,57,58,63,66,67,67,67,68,69,70,70,70,70,72,73,75,75,76,76,78,79,81)
x  = t - mean(t)
a  = -1.26
b0 = 0.0
Vb = 10
Sb = sqrt(Vb)

# Plotting the data
# -----------------
plot(x,y,xlab="Centered temperature (in Fahrenheit, mean=69.6)",ylab="O-ring failure")

# Plotting prior density and likelihood function
# ----------------------------------------------
par(mfrow=c(1,2))
bs  = seq(-15,15,length=200)
bs1 = seq(-1,0.1,length=200)
p  = NULL
l = NULL
for (i in 1:200){
  p = c(p,prior(bs[i]))
  l = c(l,like(bs1[i]))
}
plot(bs,p,type="l",xlab=expression(beta),ylab="",main=expression(p(beta)))
plot(bs1,l,type="l",xlab=expression(beta),ylab="",main=expression(l(beta,y)))

# Plotting prior*likelihood 
# -------------------------
par(mfrow=c(1,2))
bs = seq(-1,0.1,length=200)
f  = NULL
for (i in 1:200){
  f = c(f,fun(bs1[i]))
}
plot(bs,f,type="l",xlab=expression(beta),ylab="",main="")

# computing p(y)
# --------------
bs = seq(-1,0.1,by=0.0001)
f  = NULL
for (i in 1:length(bs)) f = c(f,fun(bs[i]))
py = (bs[2]-bs[1])*sum(f)

# posterior
# ---------
bs = seq(-1,0.1,length=200)
p = NULL
for (i in 1:200) p = c(p,post(bs[i]))
plot(bs,p,type="l",xlab=expression(beta),ylab="",main="")

# computing E(beta|y) and V(beta|y) - Simpson's rule
# --------------------------------------------------
h  = 0.0001
bs1= seq(-1,0.1,by=h)
m1 = 0.0
m2 = 0.0
for (i in 1:length(bs1)){
  f = fun(bs1[i])
  m1 = m1 + bs1[i]*f
  m2 = m2 + bs1[i]^2*f
}
m1 = m1*h/py
m2 = m2*h/py
V  = m2-m1^2

# computing p(y) by Monte Carlo integration
# -----------------------------------------
M     = 50000
betas = rnorm(M,b0,sqrt(Vb))
l = NULL
for (i in 1:M)
  l = c(l,like(betas[i]))
py.mc = mean(l)
mcerror = sqrt(mean((l-py.mc)^2)/M^2)

# Sampling from p(beta|y) by SIR
# ------------------------------
w = l/sum(l)
ind = sample(1:M,size=M,replace=T,prob=w)
betas1 = betas[ind]
breaks  = seq(min(betas1),max(betas1),length=40)
par(mfrow=c(1,1))
hist(betas1,breaks=breaks,prob=T,xlab=expression(beta),ylab="",
     main="",ylim=c(0,max(p)))
lines(bs,p,col=2)

# Posterior predictive
# --------------------
xs = seq(31,81,by=2)-mean(t)
ps = NULL
for (i in 1:length(xs)){
  A = exp(a+betas1*xs[i])
  ps = c(ps,mean(A/(1+A)))
}
z = xs[(ps<0.51)&(ps>0.49)]+mean(t)
plot(t,y,xlim=range(xs+mean(t)),ylab="O-ring failure",
     xlab="Centered temperature (in Fahrenheit)",axes=F)
axis(2)
axis(1,at=sort(c(seq(30,80,by=10),z)))
lines(xs+mean(t),ps,col=2,lwd=3)
segments(z,0,z,0.5,col=3,lwd=3)
segments(0,0.5,z,0.5,col=3,lwd=3)








