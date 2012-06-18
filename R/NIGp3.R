.NIGp3 <-
function(begin, alpha, delta, totalCount, distance)
{
  IGprocess <- .IGp(mu = 1/(delta*alpha), lambda=1, totalCount, distance)
  NIGprocess <- B_dt <- IG_dt <- rep(0,totalCount)
  NIGprocess[1] <- begin
  for (i in 2:totalCount)
  {
    IG_dt[i] <- IGprocess[i] - IGprocess[i-1]
    z <- rnorm(1,0,1)
    B_dt[i] <- B_dt[i-1] + sqrt(IG_dt[i])*z  #Time change of a standard Brownian motion
    NIGprocess[i] <- sum(begin,delta*B_dt[i])
    
    #See page 9 of http://iriaf.univ-poitiers.fr/colloque2011/article/v1s1a3.pdf
    #A two factor Levy model for stochastic mortality - Viou Ainou    
  }
  return(list(IGB = IG_dt[-1], NIGB = NIGprocess))
}
