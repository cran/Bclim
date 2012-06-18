plotBclimVol <-
function(x,dim=1,title=NULL,presentleft=TRUE,denscol="red",denstransp=0.5,leg=TRUE,dolog10=FALSE,med=FALSE,legloc="topleft",...) {
  #dim=1;title=NULL;presentleft=TRUE;denscol="red";denstransp=0.5;leg=TRUE;legloc="topleft"
  
  if(class(x)!="Bclim") stop("Needs a Bclim output object")
  
  # Create HDRs for volatility
  errorbar <- matrix(NA,nrow=length(x$time.grid)-1,ncol=3)
  if(dolog10) for(i in 1:(length(x$time.grid)-1)) errorbar[i,] <- quantile(log10(abs(x$vol.interp[,i,dim]*sqrt(x$ScVar[dim]))),probs=c(0.025,0.5,0.975))    
  if(!dolog10) for(i in 1:(length(x$time.grid)-1)) errorbar[i,] <- quantile(abs(x$vol.interp[,i,dim]*sqrt(x$ScVar[dim])),probs=c(0.025,0.5,0.975))
  
  # Sort out colours
  tmp <- col2rgb(denscol)
  mycol <- rgb(tmp[1,1]/255,tmp[2,1]/255,tmp[3,1]/255)
  mycol2 <- paste(mycol,as.character(as.hexmode(round(denstransp*255,0))),sep="")
  
  # Set up plot 
  par(mar=c(4,4,3,1))
  xrange <- range(c(0,x$time.grid))
  if(!presentleft) xrange <- rev(xrange)
  yrange <- range(c(0,as.vector(errorbar)))
  mytitle <- title
  if(is.null(title)) mytitle <- paste(x$core.name,": ",x$clim.dims[dim],sep="")
  plot(1,1,type="n",xlim=xrange,ylim=yrange,xlab="Age (k cal years BP)",ylab=paste(x$clim.dims[dim],ifelse(dolog10,"log10 volatility","volatility")),las=1,bty="n",main=mytitle)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="lightgray",border="NA")
  grid(col="white")  
  
  # Draw lines    
  #lines(x$time.grid[-1],errorbar[,1],col="red")
  #lines(x$time.grid[-1],errorbar[,2],col="blue")
  #lines(x$time.grid[-1],errorbar[,3],col="green")
  polygon(c(x$time.grid[-1],rev(x$time.grid[-1])),c(errorbar[,1],rev(errorbar[,3])),col=mycol2,border=mycol2)
  if(med) lines(x$time.grid[-1],apply(errorbar,1,"median"),col=denscol)
  
  # Finally draw a legend
  if(leg==TRUE) {
    legend(legloc,legend=c("95% credibility interval"),fill=c(denscol),bty="n")
  }
    
# End of function   
}
