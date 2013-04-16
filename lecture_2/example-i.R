################################################################################
#
#  Simple Monte Carlo integration 
#
################################################################################
#
# HEDIBERT FREITAS LOPES
# Associate Professor of Econometrics and Statistics
# The University of Chicago Booth School of Business
# 5807 South Woodlawn Avenue
# Chicago, Illinois, 60637
# Email : hlopes@ChicagoBooth.edu
# URL: http://faculty.chicagobooth.edu/hedibert.lopes
#
################################################################################
set.seed(723497)

par(mfrow=c(2,1))
us = seq(0,1,length=10000)
fu = (cos(50*us)+sin(20*us))^2
plot(us,fu,xlab=expression(theta),ylab=expression(h(theta)),main="",axes=F,type="l",ylim=c(0,4))
axis(1)
axis(2)

M = c(seq(10,100,by=10),
      seq(200,1000,by=100),
      seq(2000,10000,by=1000),
      seq(20000,100000,by=10000),
      seq(200000,1000000,by=100000))
M1 = c(10,100,1000,10000,100000,1000000)

u = runif(1000000)
int = (cos(50*u)+sin(20*u))^2
means  = cumsum(int)[M]/M
means2 = cumsum(int^2)[M]/M
mc.err = sqrt((means2-means^2)/M)
L = means-2*mc.err
U = means+2*mc.err

plot(log10(M),means,xlab="Sample size (log10)",ylab="p",main="",axes=F,type="b",ylim=range(c(L,U)))
axis(2)
axis(1,at=log10(M1),labels=log10(M1))
abline(h=0.965)
j = 0
for (i in M){
  j = j + 1
  points(log10(i),means[j],pch=18)
  segments(log10(i),L[j],log10(i),U[j],col=4)
}




