set.seed(1255)
n    = 1000
mu   = 0
sig2 = 1
nu   = 5
ns  = trunc(10^(seq(1,5,by=0.5)))
I1s  = NULL
I2s  = NULL
sd1s = NULL
sd2s = NULL
for (n in ns){
  y    = rgamma(n,nu/2,nu/2)
  x    = rnorm(n,mu,sqrt(sig2/y))
  I1   = mean(exp(-x^2))
  I2   = mean(1/sqrt(2*sig2/y+1))
  sd1  = sqrt(var(exp(-x^2))/n)
  sd2  = sqrt(var(1/sqrt(2*sig2/y+1))/n)
  I1s  = c(I1s,I1)
  I2s  = c(I2s,I2)
  sd1s = c(sd1s,sd1)
  sd2s = c(sd2s,sd2)
}
L1=I1s-2*sd1s
U1=I1s+2*sd1s
L2=I2s-2*sd2s
U2=I2s+2*sd2s

pdf(file="raoblackwell.pdf",width=6,height=6)
plot(I1s,xlab="Draws (log10)",ylab="Approximation",
     ylim=range(L1,L2,U1,U2),type="b",pch=16,axes=FALSE)
axis(2);box();axis(1,at=1:length(ns),label=log10(ns))
lines(I2s,col=2,pch=16,type="b")
lines(I1s-2*sd1s,col=1)
lines(I1s+2*sd1s,col=1)
lines(I2s-2*sd2s,col=2)
lines(I2s+2*sd2s,col=2)
dev.off()



