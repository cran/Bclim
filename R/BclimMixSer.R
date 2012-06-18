BclimMixSer <-
function(MDP,G=10,mixwarnings=FALSE) {
  
#   registerDoMC(cores=num.cores)  
  
  # Calculate n.samp = number of samples, n = number of layers, m = number of climate dimensions
  n.samp <- dim(MDP)[1]
  n <- dim(MDP)[2]
  m <- 3
  
  ScMean <- rep(0,m)
  ScVar <- rep(1,m)
  MDP2 <- MDP
  for(i in 1:m) {
    ScMean[i] <- mean(MDP[,,i])
    ScVar[i] <- max(var(MDP[,,i]))
    MDP2[,,i] <- (MDP[,,i]-ScMean[i])/sqrt(ScVar[i])
  }
  
  ################# MIXTURE ESTIMATION #################
  
  # Set up mixture components
  mu.mat <- array(NA,dim=c(n,m,G))
  tau.mat <- array(NA,dim=c(n,m,G))
  p.mat <- matrix(NA,nrow=n,ncol=G)
  
#   ans.all <- foreach(i = icount(n)) %dopar% {
#     cat(i,"\r")
#     Mclust(MDP2[,i,],G=G,modelName="EII",warn=mixwarnings)  
#   }
  ans.all <- list()
  for(i in 1:n) {
    cat("\r")
    cat("Completed:",format(round(100*i/n,2), nsmall = 2),"%")
    ans.all[[i]] <- Mclust(MDP2[,i,],G=G,modelNames="EII",warn=mixwarnings)
  }
  
  for(i in 1:n) {
    
    #len <- G-ncol(ans$parameters$mean)
    #if(length(len)==0) len <- 1
    #  if(len==0) {     
    mu.mat[i,,] <- ans.all[[i]]$parameters$mean
    for(g in 1:G) {
      tau.mat[i,,g] <- 1/diag(ans.all[[i]]$parameters$variance$sigma[,,g])
    }
    p.mat[i,] <- ans.all[[i]]$parameters$pro
    #  } else {
    #    ans <- Mclust(MDP2[,i,],G=nummix,modelName=ifelse(r>1,"EII","E"),warn=TRUE)
    #    MixMeans[,i,] <- ans$parameters$mean
    #    MixVar[i,] <- rep(ans$parameters$variance$sigmasq,nummix)
    #    MixPro[i,] <- ans$parameters$pro
    #}
    
    #rm(ans)
    # End of loop through n layers
  }
  
  # Now output everything to a nice neat list
  Bclimdata <- list(MDP=MDP2,n=n,m=m,n.samp=n.samp,ScMean=ScMean,ScVar=ScVar,G=G,mu.mat=mu.mat,tau.mat=tau.mat,p.mat=p.mat)
  cat("\n")
  return(Bclimdata)
  
}
