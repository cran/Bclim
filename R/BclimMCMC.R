BclimMCMC <-
function(Bclimdata,chron.loc,nchron=10000,control.mcmc=list(iterations=100000,burnin=20000,thinby=40,report=100),control.chains=list(v.mh.sd=0.8,vstart=rinvgauss(Bclimdata$n,2,1),Kstart=rep(4,Bclimdata$n)),control.priors=list(eta=rep(2.66,Bclimdata$m),phi=rep(15.33,Bclimdata$m))) {

# Create output matrices
remaining <- (control.mcmc$iterations-control.mcmc$burnin)/control.mcmc$thinby
if(remaining!=as.integer(remaining)) stop("Iterations minus burnin divided by thinby must be an integer")

vout <- rep(0,length=Bclimdata$m*(Bclimdata$n-1)*remaining)
chronout <- rep(0,length=Bclimdata$n*remaining)
cout <- rep(0,length=Bclimdata$m*(Bclimdata$n)*remaining)

# Re-dim the precisions matrix
Bclimprec <- Bclimdata$tau.mat[,1,]

#cat("\n")
#cat("Running MCMC...\n")

# Run C code
out <- .C("BclimMCMC3D", 
        as.integer(Bclimdata$G),
        as.integer(Bclimdata$n),
        as.integer(Bclimdata$m),
        as.integer(nchron),
        as.double(Bclimdata$p.mat),
        as.double(Bclimdata$mu.mat),
        as.double(Bclimprec),
        as.character(chron.loc),
        as.integer(control.mcmc$iterations),
        as.integer(control.mcmc$burnin),
        as.integer(control.mcmc$thinby),
        as.integer(control.mcmc$report),
        as.double(control.chains$v.mh.sd),
        as.double(control.chains$vstart),
        as.integer(control.chains$Kstart),
        as.double(control.priors$eta),
        as.double(control.priors$phi),
        as.double(vout),
        as.double(chronout),
        as.double(cout)
        ,PACKAGE="Bclim")

vout  <- array(NA,dim=c(remaining,Bclimdata$n-1,Bclimdata$m))
cout  <- array(NA,dim=c(remaining,Bclimdata$n,Bclimdata$m))
for(i in 1:remaining) {
    for(j in 1:Bclimdata$m) {
        vout[i,,j] <- out[[18]][seq(1,Bclimdata$n-1)+(j-1)*(Bclimdata$n-1)+(i-1)*(Bclimdata$n-1)*Bclimdata$m]
        cout[i,,j] <- out[[20]][seq(1,Bclimdata$n)+(j-1)*(Bclimdata$n)+(i-1)*(Bclimdata$n)*Bclimdata$m]
    }   
}
chronout <- matrix(out[[19]],ncol=Bclimdata$n,nrow=remaining,byrow=TRUE)

return(list(v.store=vout,chron.store=chronout,c.store=cout,eta=control.priors$eta,phi=control.priors$phi,chron.loc=chron.loc,nchron=nchron))

}
