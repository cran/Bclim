.NIGB <-
function(delta, IG.sum, tb.points, c.start, c.end){
  total.points <- length(tb.points)
  NIG.bridge <- rep(NA, total.points)
  NIG.bridge[1] <- c.start
  NIG.bridge[total.points] <- c.end
  IG.increment <- rep(NA, total.points-2)
  z <- IG.sum
  
  j <- 2
  while(j < total.points){
    # Generate chi-square random variate
    q <- rnorm(1,0,1)^2  
    
    # Reparameterise
    mu <- (tb.points[total.points]-tb.points[j]) / (tb.points[j]-tb.points[j-1])   
    lambda <- (delta^2 * (tb.points[total.points]-tb.points[j])^2) / z
    
    if(z==0) {
      print(c(delta,IG.sum,tb.points,c.start,c.end))
      stop("Problem in Inverse Gaussian Bridge")
    }
    
    # Compute the roots of the chi-square random variate
    s1 <- mu + (mu^2*q)/(2*lambda) - (mu/(2*lambda)) * sqrt(4*mu*lambda*q + mu^2*q^2)
    if(lambda<1e-5) s1 <- mu
    s2 <- mu^2 / s1
    
    # Acceptance/rejection of root
    #s <- ifelse(runif(1,0,1) < mu*(1+s1)/((1+mu)*(mu+s1)), s2, s1)
    s <- ifelse(runif(1,0,1) < mu*(1+s1)/((1+mu)*(mu+s1)), s1, s2)
    
    # Compute the IG incrrement
    IG.increment[j-1] <- z / (1 + s)
    
    # Rescale the sum of left-over distance of the IG bridge
    z <- z - IG.increment[j-1]

    #Compute the current value of the NIG bridge
    NIG.bridge[j] <- (c.start*z + c.end*(IG.sum-z)) / IG.sum + 
                        rnorm(1, 0, 1) * (IG.sum-z)*z / IG.sum
    j <- j +1
  }
  return(list(IGB = c(IG.increment, (IG.sum - sum(IG.increment))), NIGB = NIG.bridge))
}
