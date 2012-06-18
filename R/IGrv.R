.IGrv <-
function(mu, lambda, totalCount)
{
	IG.rv <- rep(NA,totalCount)
	s1 <- rep(NA,totalCount)
	s2 <- rep(NA,totalCount)
	for (i in 1:totalCount)
	{
	  v <- rnorm(1,0,1)^2
    s1[i] <- mu + mu^2 * v / (2 * lambda) - mu/ (2 * lambda) * sqrt(4 * mu * lambda * v + mu^2 * v^2)
  	s2[i] <- mu^2 / s1[i]
  	u <- runif(1,0,1)
    if(u > mu*(1+s1[i]) / ((1+mu)*(mu+s1[i])))
    {
      IG.rv[i] <- s2[i]
  	}else
  	{
      IG.rv[i] <- s1[i]
    }
  }
  return(IG.rv)
}
