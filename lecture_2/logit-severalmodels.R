# Logit link
like1  = function(b){
  theta = 1/(1+exp(-a-b*x))
  prod((theta^y)*((1-theta)^(1-y)))
}
# Probit link
like2  = function(b){
  theta = pnorm(a+b*x)
  prod((theta^y)*((1-theta)^(1-y)))
}
# Complementary log-log link
like3  = function(b){
  theta = 1-exp(-exp(a+b*x))
  prod((theta^y)*((1-theta)^(1-y)))
}

# Prior distribution for beta (temperature coefficient)
# -----------------------------------------------------
prior = function(b){dnorm(b,b0,Sb)}

# O-ring data - failure when y=1 - x is temperature (F)
# -----------------------------------------------------
y  = c(1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0)
t  = c(53,57,58,63,66,67,67,67,68,69,70,70,70,70,72,73,75,75,76,76,78,79,81)
x  = t - mean(t)
a  = -1.26             # regression intercept
b0 = 0.0               # Prior mean
Vb = 10                # Prior variance
Sb = sqrt(Vb)

# Sampling from p(beta|y) by SIR where q(beta)=prior
# --------------------------------------------------
M     = 50000
betas = rnorm(M,b0,sqrt(Vb))
l     = matrix(0,M,3)
for (i in 1:M){
  l[i,1] = like1(betas[i])
  l[i,2] = like2(betas[i])
  l[i,3] = like3(betas[i])
}

# Predictive densities
py1 = mean(l[,1])
py2 = mean(l[,2])
py3 = mean(l[,3])
py  = c(py1,py2,py3) 

# Posterior model probabilities (assuming prior Pr(M)=1/3)
pm  = py/sum(py)

# sampling from p(beta|y)
w = matrix(0,M,3)
betas1 = matrix(0,M,3)
for (i in 1:3){
  w[,i] = l[,i]/sum(l[,i])
  betas1[,i] = sample(betas,size=M,replace=T,prob=w[,i])
}

# Plotting the posterior distributions of beta
plot(density(betas1[,2],bw=0.02),xlim=c(-1,0.1),xlab="beta",ylab="",main="")
lines(density(betas1[,1],bw=0.02),col=1,lwd=2)
lines(density(betas1[,2],bw=0.02),col=2,lwd=2)
lines(density(betas1[,3],bw=0.02),col=3,lwd=2)
legend(-0.8,5,legend=c("Logit model","Probit Model","Complementary log-log model"),lty=c(1,1,1),col=1:3,bty="n",lwd=c(2,2,2))


# Computing posterior predictive, p(ynew|xnew,yold,xold)
xs = seq(min(t),max(t),by=0.25)-mean(t)
ps1 = rep(0,length(xs))
ps2 = rep(0,length(xs))
ps3 = rep(0,length(xs))
for (i in 1:length(xs)){
  ps1[i] = mean(1/(1+exp(-a-betas1[,1]*xs[i])))
  ps2[i] = mean(pnorm(a+betas1[,2]*xs[i]))
  ps3[i] = mean(1-exp(-exp(a+betas1[,3]*xs[i])))
}

# Average across competing models
ps = cbind(ps1,ps2,ps3)%*%pm

# Graphical summary
par(mfrow=c(1,1))
plot(t,y,xlim=range(xs+mean(t)),ylab="O-ring failure",
     xlab="Temperature (in Fahrenheit)",axes=F)
points(t,y,pch=16)
axis(2)
axis(1,at=seq(min(t),max(t),by=2))
lines(xs+mean(t),ps,lty=1,col=1,lwd=2)
lines(xs+mean(t),ps1,lty=2,col=2,lwd=2)
lines(xs+mean(t),ps2,lty=3,col=3,lwd=2)
lines(xs+mean(t),ps3,lty=4,col=4,lwd=2)
legend(53,0.5,legend=c("Model averaging",
"Logit model","Probit Model","Complementary log-log model"),
lty=1:4,col=1:4,lwd=rep(2,4),bty="n")

