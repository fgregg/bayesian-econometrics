#################################################################
#
#  Using the bivariate normal to illustrate the effect
#  of unnecessary MC steps, or simply the Raoblackwellization
#  of MC methods.
#
#################################################################
set.seed(1547)
mux   = 1
muy   = 2
sig2x = 1
sig2y = 1
sigxy = 0.75
sigx  = sqrt(sig2x)
sigy  = sqrt(sig2y)
sdy.x = sqrt(sig2y-sigxy^2/sig2x)


# Showing the difference in MC errors for both estimators 
# of the marginal mean
par(mfrow=c(1,3))
for (M in c(100,1000,10000)){
  E1 = rep(0,1000)
  E2 = rep(0,1000)
  for (i in 1:1000){
    x = rnorm(M,mux,sigx)
    y = muy + sigxy/sig2x*(x-mux) + sdy.x*rnorm(M)
    E1[i] = muy + sigxy/sig2x*(mean(x)-mux)
    E2[i] = mean(y)
  }
  boxplot(E2,E1,names=c("No Raoblackwellization","Raoblackwellization"),outline=FALSE,
          main=paste("MC size=",M,sep=""))
}



# Showing the difference in MC errors for both estimators 
# of the marginal density
M     = 1000
ys    = seq(-3,7,length=100)
deny  = matrix(0,1000,100)
y     = matrix(0,1000,M)
for (j in 1:1000){
   x = rnorm(M,mux,sigx)
   y[j,] = muy + sigxy/sig2x*(x-mux) + sdy.x*rnorm(M)
   for (i in 1:100)
     deny[j,i] = mean(dnorm(ys[i],muy + sigxy/sig2x*(x-mux),sdy.x))
}

par(mfrow=c(1,2))
plot(density(y[1,]),col=grey(0.8),ylim=c(0,0.5),main="No Raoblackwellization",xlab="")
for (j in 2:1000)
  lines(density(y[j,]),col=grey(0.8))
lines(ys,dnorm(ys,2,1),lwd=2)

plot(ys,deny[1,],type="l",col=grey(0.8),ylim=c(0,0.5),
     main="Raoblackwellization",xlab="",ylab="Density")
for (j in 2:1000)
  lines(ys,deny[j,],col=grey(0.8))
lines(ys,dnorm(ys,2,1),lwd=2)



