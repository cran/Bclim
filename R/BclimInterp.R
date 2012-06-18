BclimInterp <-
function(Bclim.data, Bclim.res, time.grid=seq(0,14,by=0.1)) {
  # Note that this is transforming from the IG2 parameterisation
  alpha <- sqrt(Bclim.res$phi/Bclim.res$eta)
	delta <- sqrt(Bclim.res$eta*Bclim.res$phi)
	chron <- Bclim.res$chron.store
	clim <- Bclim.res$c.store
	vol <- Bclim.res$v.store
  
  n.interp <- dim(Bclim.res$c.store)[1] # Number of samples?
  n.sections <- dim(Bclim.res$c.store)[2] # Number of layers?

  L <- n.sections + length(time.grid)
  t.all.mat <- matrix(NA, n.interp, L)
  c.all.mat <- matrix(NA, n.interp, L)
  v.all.mat <- matrix(NA, n.interp, L-1) 
  t.all.mat2 <- matrix(NA, n.interp, L)
  clim.interp <- array(NA, dim = c(n.interp, length(time.grid), Bclim.data$m))
  vol.interp <- array(NA, dim = c(n.interp, length(time.grid)-1, Bclim.data$m))
  
  for (k in 1:Bclim.data$m){
    cat("Performing interpolation for climate dimension",k, "\n")
    for (i in 1:n.interp) {
      t.all <- c()
      v.all <- c() #IGB
      c.all <- c() #NIGB
      
      if(i%%10==0) {
        cat("\r")
        #cat("Completed:",round(100*iter/iterations,2),"%")
        cat("Completed:",format(round(100*i/n.interp,2), nsmall = 2),"%")
        if(i<n.interp) cat("\r")
      }
  
      t.start <- chron[i, 1]
  		t.end <- chron[i, 2]
  		c.start <- clim[i, 1, k]
  		c.end <- clim[i, 2, k]
  
      #Move to 0 if first t value =! 0
      if(t.start != 0) {t.start <- 0}  
      t.b.points <- c(t.start, time.grid[time.grid > t.start & time.grid < t.end], t.end)
      temp <- .NIGB(delta[k], vol[i, 1, k], t.b.points, c.start, c.end)
      t.all <- c(t.all, t.b.points[-length(t.b.points)])
      v.all <- c(v.all, temp$IGB)
      c.all <- c(c.all, temp$NIGB[-length(temp$NIGB)])
  
      for (j in 2:(n.sections-2)) {
        t.start <- chron[i, j]
      	t.end <- chron[i, j+1]
    		c.start <- clim[i, j, k]
    		c.end <- clim[i, j+1, k]
        
        if(t.start > max(time.grid)) {t.start <- max(time.grid)-0.000001} 
        if(t.end > max(time.grid)) {t.end <- max(time.grid)-0.000001}  
        #if use 14, problem t.all.mat[i,] %in% t.grids (too many t points at 14)
        
        t.b.points <- c(t.start, time.grid[time.grid > t.start & time.grid < t.end], t.end)
        temp <- .NIGB(delta[k], vol[i, j, k], t.b.points, c.start, c.end)
        t.all <- c(t.all, t.b.points[-length(t.b.points)])
        v.all <- c(v.all, temp$IGB)
        c.all <- c(c.all, temp$NIGB[-length(temp$NIGB)])
        if(any(is.na(c.all))) stop("NAs in interp!")
      }
        
      t.start <- chron[i, n.sections-1]
      t.end <- chron[i, n.sections]
  		c.start <- clim[i, n.sections-1, k]
  		c.end <- clim[i, n.sections, k]
      #if the last of t < 14, do extrapolation
      if(t.end < max(time.grid)){
        t.b.points <- c(t.start, time.grid[time.grid > t.start & time.grid < max(time.grid)], max(time.grid))
        temp <- .NIGp3(c.start, alpha[ k], delta[k], length(t.b.points), (max(time.grid)-t.start))
      } else {
          t.end <- max(time.grid)
          if (t.start > max(time.grid)) {t.start <- max(time.grid)-0.000001}
          t.b.points <- c(t.start, time.grid[time.grid > t.start & time.grid < t.end], t.end)
          temp <- .NIGB(delta[k], vol[ i, j, k], t.b.points, c.start, c.end)
      }
      t.all <- c(t.all, t.b.points)
      v.all <- c(v.all, temp$IGB)
      c.all <- c(c.all, temp$NIGB)
      if(any(is.na(c.all))) stop("NAs in extrap!")
      
      t.all.mat[i, 1: length(t.all)] <- t.all  
      v.all.mat[i, 1: length(v.all)] <- v.all
      c.all.mat[i, 1: length(c.all)] <- c.all
      
    }
  
    for (i in 1:n.interp) {
      t.all.mat2[i, ] <- as.numeric(t.all.mat[i, ] %in% time.grid) # indicator if T same as t.grids
      if(sum(t.all.mat2[i, ]) > length(time.grid)) {t.all.mat2[i, 1] <- 0} #check for duplication of o's
      temp1 <- 0
      temp2 <- c()  #first volatility
      temp3 <- c()  #first climate history
      
      for(j in 1:(L-1)) {
        if (t.all.mat2[i, j] == 1) {
          temp2 <- c(temp2, temp1 + v.all.mat[i,j]) 
          temp1 <- 0
          temp3 <- c(temp3, c.all.mat[i,j])
        } else {
          temp1 <- temp1 + v.all.mat[i,j]
        }
      }
      vol.interp[i, , k] <- temp2[-length(temp2)]
      clim.interp[i, , k] <- temp3
    }
    cat("\n")
  }
  for(k in 1:Bclim.data$m) {
    clim.interp[,, k] <- clim.interp[,, k] * sqrt(Bclim.data$ScVar)[k] + Bclim.data$ScMean[k]
  }
  #cat("Completed! \n")
  return(list(clim.interp = clim.interp, vol.interp = vol.interp,time.grid=time.grid))
}
