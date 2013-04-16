set.seed(1234)
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
Sb = sqrt(c(1,10,100))

# Sampling from p(beta|y) by SIR where q(beta)=prior
# --------------------------------------------------
M     = 50000
l     = array(0,c(M,3,3))
betas = matrix(0,M,3)
for (j in 1:3){
  betas[,j] = rnorm(M,b0,Sb[j])
  for (i in 1:M){
    l[i,1,j] = like1(betas[i,j])
    l[i,2,j] = like2(betas[i,j])
    l[i,3,j] = like3(betas[i,j])
  }
}

# Predictive densities
py1 = apply(l[,1,],2,mean)
py2 = apply(l[,2,],2,mean)
py3 = apply(l[,3,],2,mean)
py  = c(py1,py2,py3) 

# Posterior model probabilities (assuming prior Pr(M)=1/9)
pm  = py/sum(py)

# sampling from p(beta|y)
w = array(0,c(M,3,3))
betas1 = array(0,c(M,3,3))
for (i in 1:3){
  for (j in 1:3){
    w[,i,j] = l[,i,j]/sum(l[,i,j])
    betas1[,i,j] = sample(betas[,j],size=M,replace=T,prob=w[,i,j])
  }
}

# Computing posterior predictive, p(ynew|xnew,yold,xold)
xs = seq(min(t),max(t),by=0.25)-mean(t)
ps = array(0,c(length(xs),3,3))
for (l in 1:length(xs)){
  for (j in 1:3){
    ps[l,1,j] = mean(1/(1+exp(-a-betas1[,1,j]*xs[l])))
    ps[l,2,j] = mean(pnorm(a+betas1[,2,j]*xs[l]))
    ps[l,3,j] = mean(1-exp(-exp(a+betas1[,3,j]*xs[l])))
  }
}

# Average across competing models
pss = cbind(ps[,1,],ps[,2,],ps[,3,])%*%pm

# Graphical summary
par(mfrow=c(1,1))
plot(t,y,xlim=range(xs+mean(t)),ylab="O-ring failure",
     xlab="Temperature (in Fahrenheit)",axes=F)
points(t,y,pch=16)
axis(2)
axis(1,at=seq(min(t),max(t),by=2))
lines(xs+mean(t),pss,lty=1,col=1,lwd=2)
for (i in 1:3)
  for (j in 1:3)
    lines(xs+mean(t),ps[,i,j],lty=j+1,col=i+1,lwd=2)

