BclimInterp <-
function(Bclim.data, Bclim.res, tgrid=seq(0,14,by=0.1)) {
	chron <- Bclim.res$chron.store  # 1000-by-150
  alpha <- sqrt(Bclim.res$phi/Bclim.res$eta)
  delta <- sqrt(Bclim.res$eta*Bclim.res$phi)
	clim <- Bclim.res$c.store   # 1000-by-150-by-3
	vol <- Bclim.res$v.store  # 1000-by-149-by-3
  
  n.interp <- dim(clim)[1] #1000
  n.sections <- dim(clim)[2] #150

	L <- n.sections + length(tgrid)  #150+141
	t.all.mat <- matrix(NA, nrow = n.interp, ncol = L) #1000-by-(150+141)
	c.all.mat <- matrix(NA, n.interp, L) #1000-by-(150+141)
	v.all.mat <- matrix(NA, n.interp, L-1) #1000-by-(150+140)
	clim.interp <- array(NA, dim = c(n.interp, length(tgrid), Bclim.data$m)) # r=3 clim dim
	vol.interp <- array(NA, dim = c(n.interp, length(tgrid)-1, Bclim.data$m))

  
  for (k in 1:Bclim.data$m){  #loop through clim dim
    cat("Perform interpolation for climate dimension", k, "\n")
    for (i in 1:n.interp) { #loop through each sample
      t.all <- c()
      v.all <- c() #IGB
      c.all <- c() #NIGB
      
      if(i%%10==0) {
        cat("\r")
        #cat("Completed:",round(100*iter/iterations,2),"%")
        cat("Completed:",format(round(100*i/n.interp,2), nsmall = 2),"%")
        if(i<n.interp) cat("\r")
      }
  
  		chron.diff.i <-  diff(chron[i,])  # rep(1, length(chron[i,]))
      t.start <- chron[i, 1]
  		t.end <- chron[i, 2]
  		c.start <- clim[i, 1, k]
  		c.end <- clim[i, 2, k]
        
      i.tmp <- length( (tgrid < chron[i,1])[(tgrid < chron[i,1])==TRUE] )
      tmp.vstart <- NULL
      if (tgrid[1]==t.start) {
        t.b.points <- c(t.start, tgrid[tgrid > t.start & tgrid < t.end], t.end)
        temp <- NIGB(delta[k]*chron.diff.i[1], vol[i, 1, k], t.b.points, c.start, c.end)
        t.all <- c(t.all, t.b.points[-length(t.b.points)]) #the last t point is the starting point for the next bridge
        v.all <- c(v.all, temp$IGB)
        c.all <- cbind(c.all, temp$NIGB[-length(temp$NIGB)]) #the last clim point is the starting point for the next bridge
        j.tmp <- 1
      } else if(i.tmp>0) { #MIRROR RESULT OF EXTRAPOLATION (NIG PROCESS), I.E. EXTRAPOLATE BACKWARD
        t.b.points <- c(tgrid[1], tgrid[tgrid > tgrid[1] & tgrid < tgrid[i.tmp+1]], tgrid[i.tmp+1])
        temp <- NIGp3(c.start, alpha[k], delta[k]*chron.diff.i[i.tmp], length(t.b.points), (t.start-tgrid[1]) )
        temp$IGB <- temp$IGB
        t.all <- c(t.all, t.b.points[-length(t.b.points)] )
        tmp.vstart <- v.all <- c(v.all, rev(temp$IGB) )
        c.all <- c(c.all, rev(temp$NIGB)[-length(temp$NIGB)] )
        j.tmp <- 0
      } else
      {
        j.tmp <- length( (chron[i,] < tgrid[1])[(chron[i,] < tgrid[1])==TRUE] ) 
        t.b.points <- c(chron[i, j.tmp], tgrid[tgrid < chron[i, j.tmp+1]], chron[i, j.tmp+1])
        temp <- NIGB(delta[k]*chron.diff.i[j.tmp], vol[i, j.tmp, k], t.b.points, c.start, c.end)
        if(chron[i, j.tmp]==tgrid[1]) {
          t.all <- c(t.all, t.b.points[-length(t.b.points)])
        } else
        {
          t.all <- c(t.all, t.b.points[-c(1, length(t.b.points))])
        }
        v.all <- c(v.all, temp$IGB)
        c.all <- c(c.all, temp$NIGB[-length(temp$NIGB)])
      }
       
      j.tmp2 <- length( (chron[i, ] <= max(tgrid))[(chron[i, ] <= max(tgrid))==T] )
      
      for (j in (j.tmp+1):(j.tmp2-1) ) {
        t.start <- chron[i, j]
      	t.end <- chron[i, j+1]
    	  c.start <- clim[i, j, k]
    	  c.end <- clim[i, j+1, k]
        
        if(t.start > max(tgrid)) {t.start <- max(tgrid)-1E-6} 
        if(t.end > max(tgrid)) {t.end <- max(tgrid)-1E-6}  
        #if use 14, problem t.all.mat[i,] %in% t.grids (too many t points at 14)
        
        t.b.points <- c(t.start, tgrid[tgrid > t.start & tgrid < t.end], t.end)
        temp <- NIGB(delta[k]*chron.diff.i[j], vol[i, j, k], t.b.points, c.start, c.end)
        t.all <- c(t.all, t.b.points[-length(t.b.points)])
        v.all <- c(v.all, temp$IGB)
        if(!is.null(tmp.vstart)) {
          v.all[1:length(tmp.vstart)] <- v.all[1:length(tmp.vstart)] + v.all[length(tmp.vstart)+1]
          tmp.vstart <- NULL
        }
        c.all <- c(c.all, temp$NIGB[-length(temp$NIGB)])
        if(any(is.na(c.all))) stop("NAs in interp!")
      }
        

      t.start <- chron[i, j.tmp2]
      c.start <- clim[i, j.tmp2, k]
      
      if(t.start < max(tgrid)){
        t.b.points <- c(t.start, tgrid[tgrid > t.start & tgrid < max(tgrid)], max(tgrid))
        tmp1 <- NIGp3(c.start, alpha[k], delta[k]*chron.diff.i[j.tmp2], length(t.b.points), (max(tgrid)-t.start))
        tmp1$IGB <- tmp1$IGB + temp$IGB[length(temp$IGB)]
        t.all <- c(t.all, t.b.points)
        v.all <- c(v.all, tmp1$IGB)
        c.all <- c(c.all, tmp1$NIGB)
      } else
      {
        t.all <- c(t.all, t.b.points[length(t.b.points)])
      	c.all <- c(c.all, temp$NIGB[length(temp$NIGB)])
      }
      if(any(is.na(c.all))) stop("NAs in extrap!")

      
      t.all.mat[i, 1: length(t.all)] <- c(t.all[-length(t.all)], max(tgrid))  
      v.all.mat[i, 1: length(v.all)] <- v.all
      c.all.mat[i, 1: length(c.all)] <- c.all
    }
  
    
    for (i in 1:n.interp) {
      tmp <- as.numeric(t.all.mat[i, ] %in% tgrid) # indicate where T () are same as time grids (141 == 140 centuries)
      temp1 <- 0
      temp2 <- c()  #first volatility
      temp3 <- c.all.mat[i, 1] # NIGB  #first climate history
      
      for(j in 2:length(tmp)) { #loop though all the grid points and observed time points
        if (tmp[j] == 1) {
          if(!is.null(temp2)) {
            temp2 <- c(temp2, temp1 + v.all.mat[i,j-1])
          } else
          {
            temp2 <- temp1 + v.all.mat[i, j-1]
          }
          temp1 <- 0
          temp3 <- c(temp3, c.all.mat[i, j])
        } else 
        {
          temp1 <- temp1 + v.all.mat[i, j-1]
        }
      }
      
      
       if(length(vol.interp[i, , k])==length(temp2)) {vol.interp[i, , k] <- temp2
          } else {
            vol.interp[i, , k] <- temp2[-length(temp2)]
          }
      if(length(clim.interp[i, , k])==length(temp3)) {clim.interp[i, , k] <- temp3
        } else {
        clim.interp[i, , k] <- c(temp3, c.all.mat[i,][!is.na(c.all.mat[i,])][length(c.all.mat[i,][!is.na(c.all.mat[i,])])]) 
        }
    }
  }
  
  for(k in 1:Bclim.data$m) {
    clim.interp[,, k] <- clim.interp[,, k] * sqrt(Bclim.data$ScVar)[k] + Bclim.data$ScMean[k]
  }
 
  cat("Completed! \n")
  return(list(clim.interp = clim.interp, vol.interp = vol.interp, time.grid=tgrid))
}
