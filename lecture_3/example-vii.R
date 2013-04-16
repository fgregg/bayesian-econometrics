############################################################################################
#
# Simulated annealing 
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
ti    = function(a,b,x){1/(1+exp(-a-b*x))}
like  =  function(k,n,y,a,b,x){sum(y*log(ti(a,b,x)))+sum((n-y)*log(1-ti(a,b,x)))}
like1 = function(theta){-like(k,n,y,theta[1],theta[2],x)}
x     = c(-0.863,-0.296,-0.053,0.727)
n     = c(5,5,5,5)
y     = c(0,1,3,5)
N1    =  50
N2    =  50
as    =  seq(-4,6,length=N1)
bs    =  seq(-10,40,length=N2)
lv    =  matrix(0,N1,N2)
for (i in 1:N1)  
  for (j in 1:N2)
    lv[i,j]  =  like(k,n,y,as[i],bs[j],x)

theta=c(5,30);  mle1  = nlm(like1,theta)
theta=c(-2,40); mle2  = nlm(like1,theta)
theta=c(-4,-10);mle3  = nlm(like1,theta)
theta=c(6,0);   mle4  = nlm(like1,theta)
rbind(c(mle1$estimate,mle1$iterations),c(mle2$estimate,mle2$iterations),
      c(mle3$estimate,mle3$iterations),c(mle4$estimate,mle4$iterations))

# Simulated Annealing
# -------------------
set.seed(923749)
thetas = matrix(c(5,-2,-4,6,30,40,-10,0),4,2)
M      =  5000
ths1   = array(0,c(4,M+1,2))
for (j in 1:4){
  thold  =  thetas[j,]
  ths    =  thold
  for (i in 1:M){
    Ti      =  1/i
    thnew   =  rnorm(2,thold,0.05)
    accept  =  min(1,exp((like(k,n,y,thnew[1],thnew[2],x)-
                         like(k,n,y,thold[1],thold[2],x))/Ti))
    if (runif(1)<accept){
      thold  =  thnew
    }
    ths = rbind(ths,thold)
  }
  ths1[j,,] = ths
}
thetas1 = ths1

set.seed(923749)
thetas = matrix(c(5,-2,-4,6,30,40,-10,0),4,2)
M      =  5000
ths1   = array(0,c(4,M+1,2))
for (j in 1:4){
  thold  =  thetas[j,]
  ths    =  thold
  for (i in 1:M){
    Ti      =  1/(10*log(1+i))
    thnew   =  rnorm(2,thold,0.05)
    accept  =  min(1,exp((like(k,n,y,thnew[1],thnew[2],x)-
                         like(k,n,y,thold[1],thold[2],x))/Ti))
    if (runif(1)<accept){
      thold  =  thnew
    }
    ths = rbind(ths,thold)
  }
  ths1[j,,] = ths
}
thetas2 = ths1


levels = c(2e-4,6e-4,0.0012,0.0018,0.0024)
ind = seq(1,M+1,by=10)
pdf(file="annealing.pdf",width=12,heigh=8)
par(mfrow=c(2,3))
plot(thetas1[1,ind,],xlab=expression(beta[1]),ylab=expression(beta[2]),main="Ti=1/i",xlim=range(as),ylim=range(bs),type="l",axes=F)
axis(1)
axis(2)
contour(as,bs,exp(lv),nlevels=5,levels=levels,drawlabels=F,add=T)
for (i in 1:4)
  lines(thetas1[i,ind,],col=i)
plot(ind,thetas1[1,ind,1],xlab="Iterations",ylab="",main=expression(beta[1]),type="l",ylim=range(as),axes=F)
axis(1)
axis(2)
for (i in 2:4)
  lines(ind,thetas1[i,ind,1],col=i)
abline(h=0.8734387)
plot(ind,thetas1[1,ind,2],xlab="Iterations",ylab="",main=expression(beta[2]),type="l",ylim=range(bs),axes=F)
axis(1)
axis(2)
for (i in 2:4)
  lines(ind,thetas1[i,ind,2],col=i)
abline(h=7.912780)

plot(thetas2[1,ind,],xlab=expression(beta[1]),ylab=expression(beta[2]),
     main="Ti=1/[10log(i+1)]",xlim=range(as),ylim=range(bs),type="l",axes=F)
axis(1)
axis(2)
contour(as,bs,exp(lv),nlevels=5,levels=levels,drawlabels=F,add=T)
for (i in 1:4)
  lines(thetas2[i,ind,],col=i)
plot(ind,thetas2[1,ind,1],xlab="Iterations",ylab="",main=expression(beta[1]),type="l",ylim=range(as),axes=F)
axis(1)
axis(2)
for (i in 2:4)
  lines(ind,thetas2[i,ind,1],col=i)
abline(h=0.8734387)
plot(ind,thetas2[1,ind,2],xlab="Iterations",ylab="",main=expression(beta[2]),type="l",ylim=range(bs),axes=F)
axis(1)
axis(2)
for (i in 2:4)
  lines(ind,thetas2[i,ind,2],col=i)
abline(h=7.912780)
dev.off()

pdf(file="annealing1.pdf",width=8,height=8)
par(mfrow=c(2,1))
plot(3000:5001,thetas1[1,3000:5001,1],xlab="Iterations",ylab="",main=expression(beta[1]),type="l",ylim=range(thetas1[,3000:5001,1]))
for (i in 2:4) lines(3000:5001,thetas1[i,3000:5001,1],col=i)
abline(h=mle1$estimate[1],col=1,lwd=3)
plot(3000:5001,thetas1[1,3000:5001,2],xlab="Iterations",ylab="",main=expression(beta[2]),type="l",ylim=range(thetas1[,3000:5001,2]))
for (i in 2:4) lines(3000:5001,thetas1[i,3000:5001,2],col=i)
abline(h=mle1$estimate[2],col=1,lwd=3)
dev.off()
