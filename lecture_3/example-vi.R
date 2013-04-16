############################################################################################
#
#  Random walk Metropolis-Hastings versus Independence Metropolis-Hastings 
# 
#  The target density if a mixture of two 2-dimensional normal densities;
#
#  Part 1 of the example implements the random walk Metropolis-Hasting algorithm:
#    Chain paths for 6 combinations of initial values and tuning parameters;
#    Chain autocorrelations for 6 combinations of initial values and tuning parameters.
#
#  Part 2 of the example implements the independence Metropolis-Hasting algorithm:
#    Chain paths for 6 combinations of initial values and tuning parameters;
#    Chain autocorrelations for 6 combinations of initial values and tuning parameters.
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
dnormm = function(x,m,P){
  exp(-0.5*(t(x-m)%*%P%*%(x-m)))
}
dmixtnormm = function(x,m1,m2,P1,P2,dP1,dP2){
  0.7*dP1*exp(-0.5*(t(x-m1)%*%P1%*%(x-m1)))+
  0.3*dP2*exp(-0.5*(t(x-m2)%*%P2%*%(x-m2)))
}

# True distribution
# -----------------
m1  = c(4,5)
m2  = c(0.7,3.5)
S1  = matrix(c(1,0.7,0.7,1),2,2)
iS1 = solve(S1)
iS1 = (iS1+t(iS1))/2
S2  = matrix(c(1,-0.7,-0.7,1),2,2)
iS2 = solve(S2)
iS2 = (iS2+t(iS2))/2
dP1 = prod(diag(chol(iS1)))
dP2 = prod(diag(chol(iS2)))
x1  = seq(-2,7,by=0.1)
x2  = seq(0,8,by=0.1)
fnorm = matrix(0,length(x1),length(x2))
for (i in 1:length(x1))
  for (j in 1:length(x2))
    fnorm[i,j] = dmixtnormm(c(x1[i],x2[j]),m1,m2,iS1,iS2,dP1,dP2)

par(mfrow=c(1,2))
contour(x1,x2,fnorm,xlab=expression(theta[1]),ylab=expression(theta[2]),main="")
persp(x1,x2,fnorm,theta=30,phi=15,expand=0.5,ltheta=120,shade=0.75,ticktype="detailed",xlab="",ylab="",zlab="")
text(0.35,-0.2,expression(theta[2]),cex=1.5)
text(-0.2,-0.3,expression(theta[1]),cex=1.5)

# Random Walk Metropolis 
# ----------------------
set.seed(1027675)
tunning = c(0.01,1,100)
init  = rbind(c(4,5),c(0,7))
M     = 5000
neff  = NULL
neff1 = NULL
draws = NULL
par(mfrow=c(2,3))
for (j in 1:2){
  for (tun in tunning){
    V = diag(tun,2)
    cV = t(chol(V))
    iV = solve(V)
    x = init[j,]
    xs = x
    accepted = 0
    for (i in 1:M){
      xnew = x+cV%*%rnorm(2)
      acceptance = min(1,(dmixtnormm(xnew,m1,m2,iS1,iS2,dP1,dP2)/
                   dmixtnormm(x,m1,m2,iS1,iS2,dP1,dP2))*(dnormm(x,xnew,iV)/dnormm(xnew,x,iV)))
      u = runif(1)
      if (u<acceptance){
        x = xnew
        accepted = accepted+1
      }
      xs = rbind(xs,t(x))
    }
    # Metropolis paths
    xss = matrix(0,2*M,2)
    xss[1,] = xs[1,1:2]
    xss[2,2] = xs[1,2]
    for (i in 2:M){
      xss[2*(i-1),1] = xs[i,1]
      xss[2*i-1  ,1] = xs[i,1]
      xss[2*i-1,2] = xs[i,2]
      xss[2*i  ,2] = xs[i,2]
    }
    plot(xs,main="",xlim=c(-2,7),ylim=c(0,8),xlab="",ylab="",pch=16,cex=0.4)
    title(paste("tuning=",tun," \n Initial value=(",init[j,1],",",init[j,2],") \n Acceptance rate=",
          round(100*accepted/M,1),"%",sep=""))
    contour(x1,x2,fnorm,main="",xlim=c(-2,7),ylim=c(0,8),drawlabels=FALSE,add=TRUE)
    neff = c(neff,M/(1+0.01*sum((200-2:200)*acf(xs[,1],lag=200,plot=F)$acf[2:200])))
    neff = c(neff,M/(1+0.01*sum((200-2:200)*acf(xs[,2],lag=200,plot=F)$acf[2:200])))
    neff1 = c(neff1,M/(1+2*sum(acf(xs[,1],lag=200,plot=F)$acf[2:200])))
    neff1 = c(neff1,M/(1+2*sum(acf(xs[,2],lag=200,plot=F)$acf[2:200])))
    draws = cbind(draws,xs)
  }
}
matrix(neff,6,2,byrow=T)
matrix(neff1,6,2,byrow=T)

tun = c(tunning,tunning)
ini = rbind(matrix(init[1,],3,2,byrow=T),
            matrix(init[2,],3,2,byrow=T))

par(mfrow=c(2,3))
j = 0
for (i in seq(2,12,2)){
  j = j + 1
  acf(draws[,i],xlab="lag",ylab="",main="",lag=200)
  title(paste("Tuning=",tun[j]," \n Initial value=(",ini[j,1],",",ini[j,2],")",sep=""))
}


# Independent Metropolis
# ----------------------
set.seed(1027675)
init = rbind(c(4,5),c(0,7))
tunning = c(0.5,5,50)
m    = c(3.01,4.55)
M    = 5000
neff = NULL
neff1 = NULL
draws = NULL
par(mfrow=c(2,3))
for (j in 1:2){
  for (tun in tunning){
    V = diag(tun,2)
    cV = t(chol(V))
    iV = solve(V)
    x = init[j,]
    xs = t(x)
    accepted = 0
    for (i in 1:M){
      xnew = m+cV%*%rnorm(2)
      acceptance = min(1,(dmixtnormm(xnew,m1,m2,iS1,iS2,dP1,dP2)/
                   dmixtnormm(x,m1,m2,iS1,iS2,dP1,dP2))*(dnormm(x,m,iV)/dnormm(xnew,m,iV)))
      u = runif(1)
      if (u<acceptance){
        x = xnew
        accepted = accepted+1
      }
      xs = rbind(xs,t(x))
    }
    # Metropolis paths
    xss = matrix(0,2*M,2)
    xss[1,] = xs[1,1:2]
    xss[2,2] = xs[1,2]
    for (i in 2:M){
      xss[2*(i-1),1] = xs[i,1]
      xss[2*i-1  ,1] = xs[i,1]
      xss[2*i-1,2] = xs[i,2]
      xss[2*i  ,2] = xs[i,2]
    }
    plot(xs,main="",xlim=c(-2,7),ylim=c(0,8),xlab="",ylab="",pch=16,cex=0.4)
    title(paste("tuning=",tun," \n Initial value=(",init[j,1],",",init[j,2],") \n Acceptance rate=",
          round(100*accepted/M,1),"%",sep=""))
    text(2.1,7.85,paste("Acceptance rate=",round(100*accepted/M,1),"%",sep=""))
    contour(x1,x2,fnorm,main="",xlim=c(-2,7),ylim=c(0,8),drawlabels=FALSE,add=TRUE)
    neff = c(neff,M/(1+0.01*sum((200-2:200)*acf(xs[,1],lag=200,plot=F)$acf[2:200])))
    neff = c(neff,M/(1+0.01*sum((200-2:200)*acf(xs[,2],lag=200,plot=F)$acf[2:200])))
    neff1 = c(neff1,M/(1+2*sum(acf(xs[,1],lag=200,plot=F)$acf[2:200])))
    neff1 = c(neff1,M/(1+2*sum(acf(xs[,2],lag=200,plot=F)$acf[2:200])))
    draws = cbind(draws,xs)
  }
}
matrix(neff,6,2,byrow=T)
matrix(neff1,6,2,byrow=T)
tun = c(tunning,tunning)
ini = rbind(matrix(init[1,],3,2,byrow=T),matrix(init[2,],3,2,byrow=T))

par(mfrow=c(2,3))
j = 0
for (i in seq(2,12,2)){
  j = j + 1
  acf(draws[,i],xlab="lag",ylab="",main="",lag=200)
  title(paste("Tuning=",tun[j]," \n Initial value=(",ini[j,1],",",ini[j,2],")",sep=""))
}


