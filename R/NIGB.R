NIGB <-
function(delta, IG.sum, tb.points, c.start, c.end){
  total.points <- length(tb.points)
  NIG.bridge <- rep(NaN, nrow=total.points)
  NIG.bridge[1] <- c.start
  NIG.bridge[total.points] <- c.end
  IG.increment <- rep(NaN, total.points-2)
  z <- IG.sum
  
  l <- 2
  while(l < total.points){
    # Generate chi-square random variate
    q <- rnorm(1,0,1)^2  
    
    # Reparameterise
    mu <- (tb.points[total.points]-tb.points[l]) / (tb.points[l]-tb.points[l-1])   
    lambda <- (1/delta^2 * (tb.points[total.points]-tb.points[l])^2) / z
    
    if(z==0) {
      print(c(delta,IG.sum,tb.points,c.start,c.end))
      stop("Problem in Inverse Gaussian Bridge")
    }

    
    # Compute the roots of the chi-square random variate
    s1 <- mu + (mu^2*q)/(2*lambda) - (mu/(2*lambda)) * sqrt(4*mu*lambda*q + mu^2*q^2)
    if(lambda<1e-5) s1 <- mu
    s2 <- mu^2 / s1
    
    # Acceptance/rejection of root
    s <- ifelse(runif(1,0,1) < mu*(1+s1)/((1+mu)*(mu+s1)), s1, s2)
    
    
    # Compute the IG incrrement
    IG.increment[l-1] <- z / (1 + s)
    
    # Rescale the sum of left-over distance of the IG bridge
    z <- z - IG.increment[l-1]
    
    #Compute the current value of the NIG bridge
    NIG.bridge[l] <- (c.start*z + c.end*(IG.sum-z)) / IG.sum +
                        rnorm(1, 0, 1) * (IG.sum-z)*z / IG.sum
    l <- l +1
  }
  return(list(IGB = c(IG.increment, (IG.sum - sum(IG.increment))), NIGB = t(NIG.bridge)))
}
