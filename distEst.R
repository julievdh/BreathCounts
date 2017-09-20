distEst<-function(n,mean,sd,lowerBound,upperBound){
  range <- upperBound - lowerBound
  m <- (mean-lowerBound) / range #mapping mean to 0-1 range
  s <- sd / range #mapping sd to 0-1 range
  a <- (m^2 - m^3 - m*s^2)/s^2 #calculating alpha for rbeta 
  b <- (m-2*m^2+m^3-s^2+m*s^2)/s^2 #calculating beta for rbeta
  data <- rbeta(n,a,b)  #generating data
  data <- lowerBound + data * range #remaping to given bounds
  return(data)
}