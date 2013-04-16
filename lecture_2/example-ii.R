################################################################################
#
#  Monte Carlo integration 
#  Monte Carlo via Importance function integration
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
set.seed(123456)

x = seq(-10,10,length=1000)
u = seq(0,0.5,length=1000)
Ms = c(100,1000,10000,100000,1000000)
stats = NULL
for (M in Ms){
  # Simple MC
  theta = rt(M,1)
  h1    = rep(0,M)
  h1[theta>2]=1
  p1    = mean(h1)
  var1  = mean((h1-p1)^2)/M
  sd1   = sqrt(var1)
  # Not so simple MC
  u    = 0.5*runif(M)
  h2   = u^(-2)/(2*pi*(1+u^(-2)))
  p2   = mean(h2)
  var2 = mean((h2-p2)^2)/M
  sd2  = sqrt(var2)
  stats = rbind(stats,c(c(p1,p2),c(round(sd1,6),round(sd2,6))))
}
stats = cbind(Ms,stats)

stats
