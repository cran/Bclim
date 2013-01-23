NIGp3 <-
function(begin, alpha, delta, totalCount, distance)
{
  dt <- distance/(totalCount-1)
  ig_dt <- rinvgauss(totalCount-1, mu = 1*dt/alpha, lambda=1*dt^2)
  NIGprocess <- B_dt <- rep(0,totalCount)
  NIGprocess[1] <- begin
  
  for (i in 2:totalCount) {
    z <- rnorm(1,0,1)
    B_dt[i] <- B_dt[i-1] + sqrt(ig_dt[i-1])*z  #Time change of a standard Brownian motion
    NIGprocess[i] <- sum(NIGprocess[i-1], 1*B_dt[i]) #delta*B_dt[i]
    
    #See page 9 of http://iriaf.univ-poitiers.fr/colloque2011/article/v1s1a3.pdf
    #A two factor Levy model for stochastic mortality - Viou Ainou    
  }
  return(list(IGB = cumsum(ig_dt), NIGB = NIGprocess))
}
