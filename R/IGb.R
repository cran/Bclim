.IGb <-
function(delta, startValue, endValue, totalCount, distance)
{
  X <- rep(NA,totalCount)
  X[1] <- startValue
  X[totalCount] <- endValue
  z <- (endValue - startValue)
  dt <- distance/(totalCount-1)
  for (i in 2:(totalCount - 1))
  {  
    #generate chi-square random variate
    q <- rnorm(1,0,1)^2    

    #reparameterise
    mu <- (totalCount - i) / (1)   #multiply both top and bottom by dt, which cancled out
    lambda <- (delta^2 * ((totalCount-i)*dt)^2) / z

    #compute the roots of the chi-square random variate
    s1 <- mu + (mu^2*q)/(2*lambda) - (mu/(2*lambda))*sqrt(4*mu*lambda*q + mu^2*q^2)
    s2 <- mu^2/s1
    
    if (runif(1,0,1) > mu*(1+s1)/((1+mu)*(mu+s1))) 
    {                  
      s <- s2
    }else
    {
      s <- s1
    }
    X[i] <- X[i-1] + (X[totalCount] - X[i-1]) / (1 + s)
    z <- z - (X[totalCount] - X[i-1]) / (1 + s)
  }
  return(X)
}
