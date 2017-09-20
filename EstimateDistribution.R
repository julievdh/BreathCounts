n <- 500
mean <- 1.54
sd <- 0.22
lowerBound <- 1.08
upperBound <- 2.18 

RoosEstdist <- distEst(n,mean,sd,lowerBound,upperBound)

hist(RoosEstdist)
