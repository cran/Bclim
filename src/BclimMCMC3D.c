// This function runs the Bclim MCMC functions

#include<R.h>
#include<Rmath.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include"use.h"
#include <R_ext/Utils.h>

#define LARGE 1.0e20
   // Define your own tolerance

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
////////////////////// Main MCMC function //////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void BclimMCMC3D(int *G,int *n, int *m,int *nchrons, double *MixPro, double *MixMeans, double *MixPrec, char **ChronsFILE, int *iterations, int *burnin, int *thinby, int *reportevery, double *vmhsd, double *vstart, int *Kstart, double *muval, double *phival, double *vstore, double *chronstore, double *cstore)
{

// G is number of mixture groups, n is number of layers, m is number of climate dimensions (always 3 here), nchrons is number of chronologies
// MixPro is mixture proportions, MixMeans is mixture means, MixPrec is mixture precisions
// ChronsFILE is the path to the chronologies file
// iterations is number of iterations, burnin is size of burnin, thinby is thinning amount, reportevery is how often to report
// vmhsd is the standard deviation of the truncated random walk for proposing new values of v
// vstart and Kstart are starting values for v and K respectively
// mu and phi are 3-vectors of the values used for mu and phi
// vstore is variance output, chronstore are used chronologies, cstore is used for climates
    
// Declare indicators
int i,j,k,l,j2,j3,iter;
	
// Get a new seed
GetRNGstate();

// Set up arrays to enable simpler mixture read-ins
double ***MyMixMean,**MyMixPrec,**MyMixPro,*CurrMixPrec;
MyMixMean = (double ***)calloc(*n, sizeof(double **));
for(i=0;i<*n;i++) MyMixMean[i] = (double **)calloc(*m,sizeof(double *));
for(i=0;i<*n;i++) for(j=0;j<*m;j++) MyMixMean[i][j] = (double *)calloc(*G, sizeof(double));

MyMixPrec = (double **)calloc(*n, sizeof(double *));
for(i=0;i<*n;i++) MyMixPrec[i] = (double *)calloc(*G, sizeof(double));
MyMixPro = (double **)calloc(*n, sizeof(double *));
for(i=0;i<*n;i++) MyMixPro[i] = (double *)calloc(*G, sizeof(double));
CurrMixPrec = (double *)calloc(*n, sizeof(double));
if(MyMixMean==NULL) error("Can't allocate memory");
if(MyMixPrec==NULL) error("Can't allocate memory");
if(MyMixPro==NULL) error("Can't allocate memory");
if(CurrMixPrec==NULL) error("Can't allocate memory");

// Turn MixPrec, MixMeans and MixPro into arrays for easy lookup
for(i=0;i<*n;i++) for(k=0;k<*G;k++) MyMixPrec[i][k] = MixPrec[k*(*n)+i];
for(i=0;i<*n;i++) {
	for(j=0;j<*m;j++) {
		for(k=0;k<*G;k++) {
			MyMixMean[i][j][k] = MixMeans[k*(*m)*(*n)+j*(*n)+i];
		}
	}
}
for(i=0;i<*n;i++) for(j=0;j<*G;j++) MyMixPro[i][j] = MixPro[j*(*n)+i];

// Set up matrices
double **v,*mu,*phi,*mixpro,*choldiag,*chollower,*rannormal,*cholzeros,*mvnvar;
int *K,*mixseq;
v = (double **)calloc(*n-1, sizeof(double *));
for(i=0;i<*n-1;i++) v[i] = (double *)calloc(*m, sizeof(double));
mu = (double *)calloc(*m, sizeof(double));
phi = (double *)calloc(*m, sizeof(double));
K = (int *)calloc(*n, sizeof(int));
mixpro = (double *)calloc(*G, sizeof(double));
for(i=0;i<*G;i++) mixpro[i]=1/((double)*G);
mixseq = (int *)calloc(*G, sizeof(int));
for(i=0;i<*G;i++) mixseq[i]=i;
choldiag = (double *)calloc(*n, sizeof(double));
chollower = (double *)calloc(*n, sizeof(double));
rannormal = (double *)calloc(*n, sizeof(double));
cholzeros = (double *)calloc(*n, sizeof(double));
mvnvar = (double *)calloc(*n, sizeof(double));

if(v==NULL) error("Can't allocate memory");
if(mu==NULL) error("Can't allocate memory");
if(phi==NULL) error("Can't allocate memory");
if(K==NULL) error("Can't allocate memory");
if(mixseq==NULL) error("Can't allocate memory");
if(choldiag==NULL) error("Can't allocate memory");
if(chollower==NULL) error("Can't allocate memory");
if(rannormal==NULL) error("Can't allocate memory");
if(cholzeros==NULL) error("Can't allocate memory");
if(mvnvar==NULL) error("Can't allocate memory");

// Give starting values
for(i=0;i<*n-1;i++) for(j=0;j<*m;j++) v[i][j] =  vstart[i];
for(i=0;i<*m;i++) mu[i] = muval[i];
for(i=0;i<*m;i++) phi[i] = phival[i];
for(i=0;i<*n;i++) K[i] = Kstart[i];    

// Set up the matrices I'll need;
double *mujk,*Djk,*VDmu,*Vz,*mujkstar,*Djkstar,*VDmustar,*Dmu,*Dmustar;
mujk = (double *)calloc(*n, sizeof(double));
Djk = (double *)calloc(*n, sizeof(double));
VDmu = (double *)calloc(*n, sizeof(double));
Vz = (double *)calloc(*n, sizeof(double));
mujkstar = (double *)calloc(*n, sizeof(double));
Djkstar = (double *)calloc(*n, sizeof(double));
VDmustar = (double *)calloc(*n, sizeof(double));
Dmu = (double *)calloc(*n, sizeof(double));
Dmustar = (double *)calloc(*n, sizeof(double));

if(mujk==NULL) error("Can't allocate memory");
if(Djk==NULL) error("Can't allocate memory");
if(VDmu==NULL) error("Can't allocate memory");
if(Vz==NULL) error("Can't allocate memory");
if(mujkstar==NULL) error("Can't allocate memory");
if(Djkstar==NULL) error("Can't allocate memory");
if(VDmustar==NULL) error("Can't allocate memory");

// Set up somewhere to store output
double **allchrons;
allchrons = (double **)calloc(*nchrons, sizeof(double *));
for(i=0;i<*nchrons;i++) allchrons[i] = (double *)calloc(*n, sizeof(double));
if(allchrons==NULL) error("Can't allocate memory");

// Read chronologies into matrix to make it easier to randomly select one
FILE *chrons;
chrons = fopen(*ChronsFILE,"r");
if(chrons==NULL) error("Error: can't open chronologies file %s.\n",*ChronsFILE);		
double tmp;	
for(i=0;i<*nchrons;i++) for(j=0;j<*n;j++) tmp=fscanf(chrons,"%lf",&allchrons[i][j]);                       
fclose(chrons);

// Print one line of chronologies if required
//for(i=1;i<*n;i++) Rprintf("%lf \n",allchrons[0][i]);

// Create chron differences
double *diffchron,*currchron;
diffchron = (double *)calloc(*n-1, sizeof(double));
currchron = (double *)calloc(*n, sizeof(double));
if(diffchron==NULL) error("Can't allocate memory");
if(currchron==NULL) error("Can't allocate memory");

// Set up MCMC variables I'll need
double *triupper,*tridiag,*trilower,vsqrtstar,vsqrtstarrat,temp,*z,tzVz,logdetratio,expnum,expden,igratio,logRv,*ones,*Vone,oneVone,exp1bit,exp2bit,*tempv,u1,u2,u3,sumv,sumxv,*tempc;
int accept,countchron=0,countv=0,countmu=0,countphi=0,currentpos,countK=0,countc=0;
triupper = (double *)calloc(*n, sizeof(double));
tridiag = (double *)calloc(*n, sizeof(double));
trilower = (double *)calloc(*n, sizeof(double));
tempv = (double *)calloc(*n-1, sizeof(double));
z = (double *)calloc(*n, sizeof(double));
ones = (double *)calloc(*n, sizeof(double));
Vone = (double *)calloc(*n, sizeof(double));
tempc = (double *)calloc(*n, sizeof(double));

// Start off the chronologies
for(i=0;i<*n;i++) currchron[i] = allchrons[0][i];
diff(currchron,n,diffchron);

// Start iterations loop
double progress=0;
for(iter=0;iter<*iterations;iter++) {
//for(i=0;i<1;i++) {

	if(iter%*reportevery==0) {
	    progress = (double) 100*iter/ *iterations;
        Rprintf("\r");
        Rprintf("Completed: %4.2f %%",progress);
        //Rprintf("Completed: %i ",iter);
	    Rprintf("\r");
	    R_FlushConsole();
	}

	// Check if someone's pressed escape
	R_CheckUserInterrupt();

    // Loop through climate dimensions;
    for(j=0;j<*m;j++) {
    
        // Create mujk,Djk and VDmu
        for(i=0;i<*n;i++) {
            mujk[i] =  MyMixMean[i][j][K[i]];
            Djk[i] = MyMixPrec[i][K[i]];
            Dmu[i] = Djk[i]*mujk[i];

        }

        for(l=0;l<*n-1;l++) tempv[l] = v[l][j];
        maketri(tempv,*n,Djk,triupper,tridiag,trilower);       
        trisolve (*n, trilower,tridiag, triupper, Dmu, VDmu);

        // Sample a new v       
        for(i=0;i<*n-1;i++) {
            
            vsqrtstar = truncatedwalk(sqrt(v[i][j]), *vmhsd, 0, LARGE);
            vsqrtstarrat = truncatedrat(sqrt(v[i][j]), *vmhsd, 0, LARGE, vsqrtstar);
            
            for(l=0;l<*n;l++) z[l] = 0.0;
            z[i] = 1.0;
            z[i+1] = -1.0;
            trisolve (*n, trilower,tridiag, triupper, z, Vz); 
            tzVz = Vz[i]-Vz[i+1];
            
            logdetratio = -0.5*log(1+(1/pow(vsqrtstar,2)-1/v[i][j])*tzVz);
            expnum = -pow(VDmu[i]-VDmu[i+1],2);
            expden = 2*((1/(1/pow(vsqrtstar,2)-1/v[i][j]))+tzVz);
            igratio = dlinvgauss2(pow(vsqrtstar,2),mu[j]*diffchron[j],phi[j]*diffchron[j])-dlinvgauss2(v[i][j],mu[j]*diffchron[j],phi[j]*diffchron[j]);
            logRv = logdetratio + expnum/expden + igratio + 0.5*(log(v[i][j])-log(pow(vsqrtstar,2)));            
                                                                                                                                                            
            accept = (int)UpdateMCMC(logRv,1,0,vsqrtstarrat);
			if(accept==1) {
                v[i][j] = pow(vsqrtstar,2);
                for(l=0;l<*n-1;l++) tempv[l] = v[l][j];
                maketri(tempv,*n,Djk,triupper,tridiag,trilower);       
                trisolve (*n, trilower,tridiag, triupper, Dmu, VDmu);
            }

 
        } // End of i loop for v update
/*

        // Create some useful things for updating mu and phi
        sumv=0; sumxv=0;
        for(i=0;i<(*n-1);i++) {
            sumv += v[i][j];
            sumxv += pow(diffchron[i],2)/v[i][j];
        }
        
        u1 = sumv;
        u2 = diffchron[*n-1]-diffchron[0]-*phipar2;
        u3 = sumxv;
         
        // Update mu
        mustar = truncatedwalk(mu[j], *mumhsd, 0, LARGE);
        mustarrat = truncatedrat(mu[j], *mumhsd, 0, LARGE, mustar);
        logRmu = ((((double)*n-1)/2)+*mupar1-1)*(log(mustar)-log(mu[j]))-0.5*(phi[j]*u1*(1/mustar-1/mu[j]))-0.5*(phi[j]*u3+2*(*mupar2))*(mustar-mu[j]); 
        mu[j] = UpdateMCMC(logRmu,mustar,mu[j],mustarrat);

        // Update phi
        //phi[j] = rgamma(((double)*n-1)/2+*phipar1,1/(u1/(2*mu[j])-u2+u3*mu[j]/2));
*/        
    } // End of j loop updating v mu and phi

	// Store everything
	if((iter%*thinby==0) & (iter>=*burnin)) {
        

        currentpos = (int)((iter-*burnin)/(*thinby));
        // Sample a chronology
        for(i=0;i<*n;i++) currchron[i] = allchrons[(int)currentpos][i];
        diff(currchron,n,diffchron);
		// Sample K
	    for(i=0;i<*n;i++) {
    	    for(l=0;l<*G;l++) mixpro[l] = MyMixPro[i][l];
    	    K[i] = sample(mixseq, *G, mixpro);
	    }        

        for(i=0;i<*n;i++) chronstore[countchron+i] = currchron[i];
        countchron += *n;

		for(j=0;j<*m;j++) {
			for(i=0;i<*n-1;i++) vstore[countv+i] = v[i][j];
			countv += *n-1;

            // Calculate some climates
			for(i=0;i<*n;i++) {
                mujk[i] =  MyMixMean[i][j][K[i]];
                Djk[i] = MyMixPrec[i][K[i]];
                Dmu[i] = Djk[i]*mujk[i];
            }
            
			for(l=0;l<*n-1;l++) tempv[l] = v[l][j];
            maketri(tempv,*n,Djk,triupper,tridiag,trilower);       

            // Create cholesky decomposition of the tridiag
            CholTriDiag(tridiag,triupper,*n,choldiag,chollower);

            // Generate n independent random normals
            for (l=0;l<*n;l++) rannormal[l] = rnorm(0,1);

            // Solve the cholesky decomposition
            trisolve (*n, chollower,choldiag, cholzeros, rannormal, mvnvar);

            // Get the mean
            trisolve (*n, trilower,tridiag, triupper, Dmu, VDmu);

            // Output the climates 
            for(i=0;i<*n;i++) tempc[i] = VDmu[i] + mvnvar[i];       
                                                                                                                        
            // Store the climates in the right place
            for(i=0;i<*n;i++) cstore[countc+i] = tempc[i];
            countc += *n;

		}
    
    // End of storage if statement
    }
        
// End of iterations loop
}

// Change the RNG state
PutRNGstate();

Rprintf("\r");
R_FlushConsole();
Rprintf("Completed: 100.00 %%");
Rprintf("\n");
R_FlushConsole();
// End of function
}


/*
double *temppro,*tempphi,*temptheta,*diffthetas,sigmasqnew,alphanew,deltanew,alpharat;
int accept,countth=0,countsig=0,counta=0,countd=0,countchron=0;
diffthetas = (double *)calloc(*m-1, sizeof(double));
temppro = (double *)calloc(*nummix, sizeof(double));
tempphi = (double *)calloc(*m-1, sizeof(double));
temptheta = (double *)calloc(*m, sizeof(double));
if(diffthetas==NULL) error("Can't allocate memory");
if(temppro==NULL) error("Can't allocate memory");
*/

        
      //  if(v[148][0]==100) {
        //    for(l=0;l<*n;l++) Rprintf("trilower=%lf tridiag=%lf triupper=%lf Dmu=%lf VDmu=%lf \n",trilower[l],tridiag[l],triupper[l],Dmu[l],VDmu[l]);
          //  error("BAD");
        //}

    /*if(iter=20) {
        Rprintf("\n v[i][j] = \n");
        for(i=0;i<*n-1;i++) {
            for(j=0;j<*m;j++) Rprintf("%lf ",v[i][j]);
            Rprintf("\n");
        }
        error("Got to here");
    }*/
        /*if(i==148 & j==0) {
            Rprintf("\n v=%lf, vstar=%lf VDmu[i]%lf VDmu[i+1]=%lf \n",v[i][j],vstar,VDmu[i],VDmu[i+1]);
            //error("Got to here");
        }*/
/*
Rprintf("mu[j]=%lf, mustar=%lf, mumhsd=%lf, logRmu=%lf \n",mu[j],mustar,*mumhsd,logRmu);
        error("Good?\n");
        for(i=0;i<*n-1;i++) {
    for(j=0;j<*m;j++) Rprintf("%lf ",v[i][j]);
    Rprintf("\n");
    }
    if(iter=500) error("V's ok?\n");
     for(i=0;i<(*n-1);i++) {
            Rprintf("%i v=%lf \n",i,v[i][j]);
        } 
       Rprintf("mu[j]=%lf, mustar=%lf, mumhsd=%lf, logRmu=%lf \n",mu[j],mustar,*mumhsd,logRmu);
       Rprintf("j=%i, phi[j]=%lf \n",j,phi[j]);
               Rprintf("j=%i, mu[j]=%lf, mustar=%lf, mumhsd=%lf, logRmu=%lf, bit1=%lf, bit2=%lf, bit3=%lf \n",j,mu[j],mustar,*mumhsd,logRmu,((((double)*n-1)/2)+*mupar1-1)*(log(mustar)-log(mu[j])),0.5*(phi[j]*u1*(1/mustar-1/mu[j])),0.5*(phi[j]*u3+2*(*mupar2))*(mustar-mu[j]));
        Rprintf("phi=%lf, par1=%lf, par2=%lf \n",phi[j],((double)*n-1)/2+*phipar1,1/(u1/(2*mu[j])-u2+u3*mu[j]/2));
Rprintf("Got to here \n");
     for(l=0;l<*n;l++) Rprintf("trilower=%lf tridiag=%lf triupper=%lf Dmu=%lf VDmu=%lf \n",trilower[l],tridiag[l],triupper[l],Dmu[l],VDmu[l]);
	    // for(l=0;l<*n;l++) Rprintf("tempv=%lf *n=%i Djk=%lf diffchron=%lf \n",tempv[l],*n,Djk[l],diffchron[l]);
		for(i=0;i<*n-1;i++) {
			for(l=0;l<*m;l++)Rprintf("v=%lf ",v[i][l]);
			Rprintf("\n");
		}

*/
/*
    // Update K

    for(i=0;i<*n;i++) {
        temp = sample(mixseq, *G, mixpro);
        kstar = mixseq[(int)temp];
        for(l=0;l<*n;l++) Kstar[l] = K[l];
        Kstar[i]=kstar;


        logRk = 0;
        for(j=0;j<*m;j++) {
            for(l=0;l<*n;l++) {
                mujk[l] =  MyMixMean[l][j][(int)K[l]];

                Djk[l] = MyMixPrec[l][(int)K[l]];
            }
            for(l=0;l<*n;l++) Dmu[l] = mujk[l]*Djk[l];
            for(l=0;l<*n-1;l++) tempv[l] = v[l][j];


            maketri(tempv,*n,Djk,triupper,tridiag,trilower);       
			
            trisolve (*n, trilower,tridiag, triupper, Dmu, VDmu);
			
            for(l=0;l<*n;l++) {
                mujkstar[l] =  MyMixMean[l][j][(int)Kstar[l]];

                Djkstar[l] = MyMixPrec[l][(int)Kstar[l]];
            }
			
            for(l=0;l<*n;l++) Dmustar[l] = mujkstar[l]*Djkstar[l];
            for(l=0;l<*n-1;l++) tempv[l] = v[l][j];

            maketri(tempv,*n,Djkstar,triupper,tridiag,trilower);       
			
                                  
            trisolve (*n, trilower,tridiag, triupper, Dmustar, VDmustar);
            for(l=0;l<*n;l++) ones[l] = 0.0;

            ones[i] = 1.0;
            trisolve (*n, trilower,tridiag, triupper, ones, Vone);

            oneVone = Vone[i];
            exp1bit = 0.0;
            for(l=0;l<*n;l++) exp1bit = mujkstar[l]*Djkstar[l]*mujkstar[l]-mujk[l]*Djk[l]*mujk[l];

            exp2bit = 0.0;
            for(l=0;l<*n;l++) exp2bit = mujk[l]*Djk[l]*VDmu[l]-mujkstar[l]*Djkstar[l]*VDmustar[l];
            logdetratio = -log(1+(Djkstar[i]-Djk[i])*oneVone);
            
            logRk += - 0.5*exp1bit-0.5*exp2bit+0.5*logdetratio + log(MyMixPro[i][(int)kstar])-log(MyMixPro[i][(int)K[i]])+0.5*log(Djkstar[i])-0.5*log(Djk[i]);

            
        } // End of j loop updating k
    
        K[i] = UpdateMCMC(logRk,kstar,K[i],1.0); 
            

    } // End of i loop updating k
    
    for(i=*n-10;i<*n;i++) {
                for(l=*n-10;l<*n;l++) Rprintf("%4.2lf ",bigQ[i][l]);
                //Rprintf("%4.2lf ",trilower[i]);
                Rprintf("\n");
            }

//            for(l=0;l<*n;l++) Rprintf("diag=%lf, upp=%lf, choldiag=%lf, chollower=%lf \n",tridiag[l],triupper[l],choldiag[l],chollower[l]);
            for(l=0;l<*n;l++) Rprintf("tempc = %lf \n",tempc[l]);
            
*/    
