.IGp <-
function(mu, lambda, totalCount, distance)
{
  IG.process <- IG.rv <- rep(0,totalCount)
  dt <- distance/(totalCount-1)
  #IGvariate <- IGrv(mu*dt, lambda*dt^2, totalCount)
  #IGprocess <- cumsum(IGvariate)
  for (i in 2:totalCount)
  {
    IG.rv[i] <- .IGrv(mu*dt, (lambda*dt)^2, 1)
    IG.process[i] <- IG.process[i-1] + IG.rv[i]  
  }
  return(IG.process)
}
