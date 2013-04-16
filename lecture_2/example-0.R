# Approximating pi via Monte Carlo
# --------------------------------
pdf(file="pi0.pdf",width=6,height=6)
xxx = seq(-1,1,length=1000)
fx = sqrt(1-xxx^2)
plot(xxx,fx,xlab="",ylab="",main="",xlim=c(-1,1),ylim=c(-1,1),
     type="l",axes=FALSE)
lines(xxx,-fx)
segments(1,-1,1,1)
segments(1,1,-1,1)
segments(-1,1,-1,-1)
segments(-1,-1,1,-1)
points(0,0,pch=16)
segments(-1,0,0,0)
text(-0.5,0.05,"R=1")
dev.off()


set.seed(1245)
pdf(file="pi4.pdf",width=6,height=6)
Ns = c(100,1000,10000,100000)
par(mfrow=c(2,2))
for (N in Ns){
xxx = seq(-1,1,length=1000)
fx = sqrt(1-xxx^2)
plot(xxx,fx,xlab="",ylab="",main="",xlim=c(-1,1),ylim=c(-1,1),
     type="l",axes=FALSE)
lines(xxx,-fx)
segments(1,-1,1,1)
segments(1,1,-1,1)
segments(-1,1,-1,-1)
segments(-1,-1,1,-1)
u = -1+2*matrix(runif(2*N),N,2)
ind = apply(u^2,1,sum)<1
for (i in 1:N)
  points(u[i,1],u[i,2],pch=16,col=1+ind[i],cex=0.5)
title(paste("MC approximation to pi=",round(4*mean(ind),4),
      "\n Based on ",N," draws",sep=""))
}
dev.off()

set.seed(1245)
pis = NULL
sds = NULL
Ns  = trunc(10^(seq(1,5,by=0.5)))
for (N in Ns){
  u   = -1+2*matrix(runif(2*N),N,2)
  ind = apply(u^2,1,sum)<1
  h   = 4*ind
  pis = c(pis,mean(h))
  sds = c(sds,sqrt(var(h)/N))
}

pdf(file="pi5.pdf",width=6,height=6)
par(mfrow=c(1,1))
L = min(pis-2*sds)
U = max(pis+2*sds)
plot(pis,ylab="MC approximation",xlab="Draws (log10)",type="b",
     axes=FALSE,ylim=c(L,U),pch=16)
lines(pis+2*sds,col=2,type="b",pch=16)
lines(pis-2*sds,col=2,type="b",pch=16)
axis(2);box();axis(1,label=round(log10(Ns),1),at=1:length(Ns))
abline(h=pi,col=4,lty=2)
dev.off()

# Approximating pi via simple Riemann Sum
Ns  = seq(5,20,by=1)
pid = NULL
for (N in Ns){
  xxx = seq(0,10,length=N+1)
  h = xxx[2]-xxx[1]
  xxx = seq(h/2,10,by=h)
  pid = c(pid,2*(h*sum(exp(-0.5*xxx^2)))^2)
}

pdf(file="pi-det.pdf",width=6,height=6)
plot(Ns,pid,xlab="N",ylab="RS Approximation",pch=16,ylim=range(pid,pi))
abline(h=pi,col=4,lty=2)
dev.off()






