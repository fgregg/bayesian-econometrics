################################################################################
#
#  Rejection method & Sampling importance resampling (SIR) 
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
xxx = seq(-6,6,length=1000)

# Rejection methods: Uniform proposal
par(mfrow=c(1,2))
set.seed(1231354)
Au     = 7.89
N      = 10000
draw1  = runif(N,-10,10)
accept = dnorm(draw1)/(Au*0.05)
u      = runif(N)
draws1 = draw1[u<accept]
hist(draws1,xlab="",ylab="",main="UNIFORM PROPOSAL",prob=TRUE,ylim=c(0,0.45),xlim=c(-5,5))
lines(xxx,dnorm(xxx),col=2,lwd=2)
length(draws1)/length(draw1)

# Rejection methods: Cauchy proposal
set.seed(1231354)
Ac     = 1.53
N      = 10000
draw1  = rt(N,1)
accept = dnorm(draw1)/(Ac*dt(draw1,1))
u      = runif(N)
draws1 = draw1[u<accept]
hist(draws1,xlab="",ylab="",main="CAUCHY PROPOSAL",prob=TRUE,ylim=c(0,0.45),xlim=c(-5,5))
lines(xxx,dnorm(xxx),col=2,lwd=2)
length(draws1)/length(draw1)


# SIR: Uniform proposal
N      = 10000
m      = 2000
par(mfrow=c(1,2))
set.seed(1231354)
draw1  = runif(N,-10,10)
w      = dnorm(draw1)/0.05
ind    = sample(1:N,size=m,replace=TRUE,prob=w)
draws1 = draw1[ind]
par(mfrow=c(1,2))
hist(draws1,xlab="",ylab="",main="UNIFORM PROPOSAL",prob=TRUE,ylim=c(0,0.45),xlim=c(-5,5))
lines(xxx,dnorm(xxx),col=2,lwd=2)
length(unique(draws1))/m

# SIR: Cauchy proposal
set.seed(1231354)
draw1  = rt(N,1)
w      = dnorm(draw1)/dt(draw1,1)
ind    = sample(1:N,size=m,replace=TRUE,prob=w)
draws1 = draw1[ind]
hist(draws1,xlab="",ylab="",main="CAUCHY PROPOSAL",prob=TRUE,ylim=c(0,0.45),xlim=c(-5,5))
lines(xxx,dnorm(xxx),col=2,lwd=2)
length(unique(draws1))/m




