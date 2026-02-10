// Written by Ryan Hayes and Matt Speranza 2026-02-02

// Usage:
// impcons [N] [Mode Options] [ImpCons Options]
// [N] Number of samples
// [Mode] gimp or sample (gimp populates the G_imp directory, sample creates mc_sample.dat)
// [gimp Options] square root of number of 1D bins, or number of 2D bins in each direction
// [ImpCons] style of implicit constraints: fnex2011, fnexdozen2024, or fnpwise2026
// [fnex2011 Options] fnexp_c
// [fnexdozen2024 Options] fnexp_c
// [fnpwise2026 Options] pwise_k pwise_w
// nsubs is required in execution directory

/*
Old G_imp documentation:
These G_imp values are subtracted from the free energy profiles before
flattening. This is done because the implicit constraints focus sampling
on the alchemical endpoints without adding barriers in alchemical space,
and we desire to preserve this feature rather than flatten it out too.

This directory contains the free energy for the implicit constraints (in
units of kT) under standard conditions so that the casual user does not
need to calculate it, and these values are used by default, unless the
user passes the path for another G_imp directory to runflat and
postprocess.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "impcons.h"

#define kB 0.00198614L

#define MAXSTRING 1024
#define BLOCK 512

#define USEGPU

#ifdef USEGPU
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
// From http://stackoverflow.com/questions/16077464/atomicadd-for-real-on-gpu
// And https://stackoverflow.com/questions/37566987/cuda-atomicadd-for-doubles-definition-error
__device__ static inline
double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                                          (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                        __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
#endif
#endif

double randDouble()
{
  return (rand()+0.5)/(RAND_MAX+1.0);
}

// Start fnexp functions----
void initializeFnexp(struct_plmd* plmd, real* thetas, int nSub)
{
  int i;
  thetas[0]=M_PI/2;
  for (i=1; i<nSub; i++) {
    thetas[i]=3*M_PI/2;
  }
}

void metropolisFnexp(struct_plmd* plmd, real* thetas, int nSub)
{
  int i;
  real thetaNew;
  for (i=0; i<nSub; i++) {
    thetaNew=2*M_PI*randDouble();
    thetas[i]=thetaNew;
  }
}

void initializeDozen(struct_plmd* plmd, real* thetas, int nSub)
{
  int i;
  real b=1;
  for (i=0; i<50; i++) {
    b=0.5*log(0.25*b*nSub*nSub*M_PI/2);
    if (!(b>0)) b=0;
  }
  plmd->dozen_b=b;

  thetas[0]=M_PI/2;
  for (i=1; i<nSub; i++) {
    thetas[i]=3*M_PI/2;
  }
}

void metropolisDozen(struct_plmd* plmd, real* thetas, int nSub)
{
  int i;
  real b=plmd->dozen_b;
  real st, eOld, eNew, thetaNew;
  for (i=0; i<nSub; i++) {
    st=(-0.5*sin(thetas[i])+0.5);
    eOld=-b*st*st*st*st;

    thetaNew=2*M_PI*randDouble();
    st=(-0.5*sin(thetaNew)+0.5);
    eNew=-b*st*st*st*st;

    if (exp(eOld-eNew)>randDouble()) {
      thetas[i]=thetaNew;
    }
  }
}

void thetasToLambdasFnexp(struct_plmd* plmd, real* thetas, real* lambdas, int nSub)
{
  int i;
  real c=plmd->fnexp_c;
  real norm=0;
  for (i=0; i<nSub; i++) {
    norm+=exp(c*sin(thetas[i]));
  }
  for (i=0; i<nSub; i++) {
    lambdas[i]=exp(c*sin(thetas[i]))/norm;
  }
}
// End fnexp functions ----

// Start pwise functions----
void initializePwise(struct_plmd* plmd, real* thetas, int nSub)
{
  int i;

  thetas[0]=1;
  for (i=1; i<nSub; i++) {
    thetas[i]=0;
  }
}

real findMaxTheta(real* thetas, int nSub)
{
  real maxTheta = thetas[0];
  for (int i = 1; i < nSub; i++) {
    if (thetas[i] > maxTheta) {
      maxTheta = thetas[i];
    }
  }
  return maxTheta;
}

real fPwise(real* x, int nSub)
{
  real sum = -1;

  for (int i = 0; i < nSub; ++i) {
    if (x[i] >= 0) {
      x[i] = x[i]*x[i]*x[i];
      sum += x[i];
    } else {
      x[i] = 0;
    }
  }
  return sum;
} 

real fPrimePwise(real* xP, int nSub)
{
  real fPrimeSum = 0;

  for (int i = 0; i < nSub; ++i) {
    if (xP[i] >= 0) {
      xP[i] = 3.0 * xP[i] * xP[i];
      fPrimeSum += xP[i];
    } else {
      xP[i] = 0;
    }
  }
  return -1 * fPrimeSum;
}

real newtonsMethod(real x0, real* thetas, int nSub, real tolerance = 1e-12, real epsilon = 1e-7, int maxIter = 100)
{
  int newtonsIter = 0; 
  real y = 0; 
  real yPrime = 0; 

  real *deltaThetas = (real*)calloc(nSub, sizeof(real));
  real *deltaThetasP = (real*)calloc(nSub, sizeof(real));

  while (newtonsIter < maxIter) {
    newtonsIter++; 

    for (int i = 0; i < nSub; i++) {
      deltaThetas[i] = thetas[i] - x0;
      deltaThetasP[i] = thetas[i] - x0; 
    }

    y = fPwise(deltaThetas, nSub);
    yPrime = fPrimePwise(deltaThetasP, nSub); 

    real x1 = x0 - (y / yPrime); 
    real difference = x1 - x0; 

    if (abs(difference) < tolerance) {
      free(deltaThetas); 
      free(deltaThetasP); 
      //printf("Iter: %d\n", newtonsIter);
      return x1; 
    }

    x0 = x1; 
  }

  free(deltaThetas); 
  free(deltaThetasP); 
  return x0;
}

void thetasToLambdasPwise(struct_plmd* plmd, real* thetas, real* lambdas, int nSub)
{
  real x0 = findMaxTheta(thetas, nSub) - 1;
  real theta0 = newtonsMethod(x0, thetas, nSub); 

  for (int i = 0; i < nSub; i++) {
    if (thetas[i] >= theta0) {
      real diff = (thetas[i] - theta0);
      lambdas[i] = diff * diff * diff;
    } else {
      lambdas[i] = 0;
    }
  }
}

real getBias(real thetas, real n, real w, int nSub)
{
  real bias = 0;
  if (thetas <= (-nSub * w)) {
    real dist = (thetas + (nSub * w));
    bias = 0.5 * n * dist * dist;
  } else if (thetas >= (nSub * w)) {
    real dist = (thetas - (nSub * w));
    bias = 0.5 * n * dist * dist;
  } else {
    bias = 0;
  }
  return bias;
}

// updates thetas 
void metropolisPwise(struct_plmd* plmd, real* currentThetas, int nSub)
{
  real w = plmd->pwise_w;
  real n = plmd->pwise_k;
  real stepSize = nSub*w;
  for (int i = 0; i < nSub; i++) {
    real delta = (randDouble()-.5)*2*stepSize;
    real trial = currentThetas[i] + delta;
    real eOld = getBias(currentThetas[i], n, w, nSub);
    real eNew = getBias(trial, n, w, nSub);
    if (exp(eOld - eNew) > randDouble()) {
      currentThetas[i] = trial;
    } 
  } 
}
// End pwise functions ----

void monte_carlo_Z(struct_plmd* plmd)
{
  int ibeg,iend,Ns;
  int Neq=plmd->B/10;
  int Nmc=plmd->B;
  real *theta;
  int s,i;

  theta=(real*) calloc(plmd->nblocks,sizeof(real));

  for (s=0; s<plmd->nsites; s++) {
    ibeg=plmd->block0[s];
    iend=plmd->block0[s+1];
    Ns=iend-ibeg;

    if (plmd->impcons==fnex2011) {
      initializeFnexp(plmd,theta,Ns);
    } else if (plmd->impcons==fnexdozen2024) {
      initializeDozen(plmd,theta,Ns);
    } else if (plmd->impcons==fnpwise2026) {
      initializePwise(plmd,theta,Ns);
    }

    for (i=-Neq; i<Nmc; i++) {
      if (i%Neq==0) {
        fprintf(stdout,"Partition Function Sample Step %d\n", i); 
      }

      if (plmd->impcons==fnex2011) {
        metropolisFnexp(plmd, theta, Ns);
      } else if (plmd->impcons==fnexdozen2024) {
        metropolisDozen(plmd, theta, Ns);
      } else if (plmd->impcons==fnpwise2026) {
        metropolisPwise(plmd, theta, Ns);
      }
      if (i>=0) {
        if (plmd->impcons==fnex2011) {
          thetasToLambdasFnexp(plmd, theta, &plmd->mc_lambda[plmd->nblocks*i + ibeg], Ns);
        } else if (plmd->impcons==fnexdozen2024) {
          thetasToLambdasFnexp(plmd, theta, &plmd->mc_lambda[plmd->nblocks*i + ibeg], Ns); // Dozen and Fnexp have same function
        } else if (plmd->impcons==fnpwise2026) {
          thetasToLambdasPwise(plmd, theta, &plmd->mc_lambda[plmd->nblocks*i + ibeg], Ns);
        }
      }
    }
  }

  free(theta);
}

struct_plmd* setup(int argc, char *argv[])
{
  struct_plmd *plmd;
  int i,j,k,c;
  FILE *fp;
  char token[MAXSTRING];

  plmd=(struct_plmd*) malloc(sizeof(struct_plmd));

  fp=fopen("nsubs","r");
  if (fp==NULL) {
    fprintf(stderr,"Error, file nsubs is required\n");
    exit(1);
  }
  plmd->nsites=0;
  while (fscanf(fp,"%d",&i)==1) {
    plmd->nsites++;
  }
  fclose(fp);

  fp=fopen("nsubs","r");
  i=0;
  plmd->nblocks=0;
  plmd->nsubs=(int*) calloc(plmd->nsites,sizeof(int));
  for(i=0; i<plmd->nsites; i++) {
    fscanf(fp,"%d",&(plmd->nsubs[i]));
    plmd->nblocks+=plmd->nsubs[i];
  }
  fclose(fp);

  plmd->block0=(int*) calloc(plmd->nsites+1,sizeof(int));
  k=0;
  for(i=0; i<plmd->nsites; i++) {
    plmd->block0[i]=k;
    for(j=0; j<plmd->nsubs[i]; j++) {
      k++;
    }
  }
  plmd->block0[i]=k;

  cudaMalloc(&(plmd->block0_d),(plmd->nsites+1)*sizeof(int));
  cudaMemcpy(plmd->block0_d,plmd->block0,(plmd->nsites+1)*sizeof(int),cudaMemcpyHostToDevice);

  c=1;

  if (argc>c && sscanf(argv[c],"%d",&plmd->B)==1) {
    c++;
  } else {
    fprintf(stderr,"Error, argument %d must be a number of Monte Carlo samples\n",c);
    exit(1);
  }

  if (argc>c && sscanf(argv[c],"%s",token)==1) {
    if (strcmp("gimp",token)==0) {
      plmd->mode=gimp;
      c++;
      if (argc>c && sscanf(argv[c],"%d",&plmd->NBINS2)==1) {
        plmd->NBINS=plmd->NBINS2*plmd->NBINS2;
        c++;
      } else {
        fprintf(stderr,"Error, argument %d after gimp must be number of 2D bins and square root of number of 1D bins\n",c);
        exit(1);
      }
    } else if (strcmp("sample",token)==0) {
      plmd->mode=sample;
      c++;
    } else {
      fprintf(stderr,"Error, mode in argument %d is not gimp or sample\n",c);
      exit(1);
    }
  } else {
    fprintf(stderr,"Error, argument %d must specify gimp or sample as the mode\n",c);
    exit(1);
  }

  if (argc>c && sscanf(argv[c],"%s",token)==1) {
    if (strcmp("fnex2011",token)==0) {
      plmd->impcons=fnex2011;
      c++;
    } else if (strcmp("fnexdozen2024",token)==0) {
      plmd->impcons=fnexdozen2024;
      c++;
    } else if (strcmp("fnpwise2026",token)==0) {
      plmd->impcons=fnpwise2026;
      c++;
    } else {
      fprintf(stderr,"Error, impcons in argument %d is not fnex2011, fnexdozen2024, or fnpwise2026\n",c);
      exit(1);
    }
  } else {
    fprintf(stderr,"Error, argument %d must specify fnex2011, fnexdozen2024, or fnpwise2026 as the implicit constraint type\n",c);
    exit(1);
  }

  if (plmd->impcons==fnex2011 || plmd->impcons==fnexdozen2024) {
    if (argc>c && sscanf(argv[c],"%lg",&plmd->fnexp_c)==1) {
      c++;
    } else {
      fprintf(stderr,"Error, argument %d must contain fnex c value\n",c);
      exit(1);
    }
  } else if (plmd->impcons==fnpwise2026) {
    if (argc>c+1 && sscanf(argv[c],"%lg",&plmd->pwise_k)==1 && sscanf(argv[c+1],"%lg",&plmd->pwise_w)==1) {
      c+=2;
    } else {
      fprintf(stderr,"Error, arguments %d and %d must contain pwise k and w values\n",c,c+1);
      exit(1);
    }
  }

  plmd->mc_lambda=(real*) calloc(plmd->nblocks*plmd->B,sizeof(real));

#ifdef USEGPU
  plmd->mc_ensweight=(real*)calloc(plmd->B,sizeof(real));
  for (i=0; i<plmd->B; i++) {
    plmd->mc_ensweight[i]=1;
  }

  cudaMalloc(&(plmd->mc_lambda_d),plmd->B*plmd->nblocks*sizeof(real));
  cudaMalloc(&(plmd->mc_ensweight_d),plmd->B*sizeof(real));

  cudaMemcpy(plmd->mc_ensweight_d,plmd->mc_ensweight,plmd->B*sizeof(real),cudaMemcpyHostToDevice);

  // count nprof
  plmd->nprof=0;
  for (i=0; i<plmd->nsites; i++) {
    for (j=i; j<plmd->nsites; j++) {
      if (i==j) {
        if (plmd->nsubs[i]==2) {
          plmd->nprof+=plmd->nsubs[i]+plmd->nsubs[i]*(plmd->nsubs[i]-1)/2;
        } else {
          plmd->nprof+=plmd->nsubs[i]+2*plmd->nsubs[i]*(plmd->nsubs[i]-1)/2;
        }
      } else { // if (plmd->msprof)
        plmd->nprof+=plmd->nsubs[i]*plmd->nsubs[j];
      }
    }
  }
#endif
  cudaMalloc(&(plmd->mc_Zprofile_d),plmd->nprof*plmd->NBINS*sizeof(real));
  return plmd;
}


#ifdef USEGPU
/*
__device__
void reduceBitonicSort(int itmp,real Ztmp,int* iloc,real* Zloc,real* Zloc2,real* Zglobal)
{
  int i1,i2;
  int direction,otherThreadIdx,iother,bswitch;

  if (threadIdx.x<NBINS) {
    Zloc2[threadIdx.x]=0;
  }
  // Bitonic sort
  for (i1=1; i1<BLOCK; i1*=2) {
    direction=(((2*i1)&threadIdx.x)!=0); // 0 ascending, 1 descending
    for (i2=i1; i2>0; i2/=2) {
      otherThreadIdx=(threadIdx.x^i2);
      if (i2<32) {
        iother=__shfl_xor_sync(-1,itmp,i2);
        bswitch=(((otherThreadIdx>threadIdx.x)==(iother>itmp))==direction);
        bswitch=(iother==itmp?0:bswitch);
        itmp=__shfl_sync(-1,itmp,threadIdx.x^(i2*bswitch));
        Ztmp=__shfl_sync(-1,Ztmp,threadIdx.x^(i2*bswitch));
      } else {
        iloc[threadIdx.x]=itmp;
        Zloc[threadIdx.x]=Ztmp;
        __syncthreads();
        iother=iloc[otherThreadIdx];
        bswitch=(((otherThreadIdx>threadIdx.x)==(iother>itmp))==direction);
        bswitch=(iother==itmp?0:bswitch);
        itmp=iloc[threadIdx.x^(i2*bswitch)];
        Ztmp=Zloc[threadIdx.x^(i2*bswitch)];
        __syncthreads();
      }
    }
  }
  iloc[threadIdx.x]=itmp;
  Zloc[threadIdx.x]=Ztmp;
  __syncthreads();
  // Reduction
  for (i1=1; i1<BLOCK; i1*=2) {
    if ((threadIdx.x&i1) && (threadIdx.x&(i1-1))==0) {
      if (itmp==iloc[threadIdx.x-i1]) {
        Zloc[threadIdx.x-i1]+=Zloc[threadIdx.x];
      } else {
        Zloc2[itmp]+=Zloc[threadIdx.x];
      }
    }
    __syncthreads();
  }
  if (threadIdx.x==0) {
    Zloc2[itmp]+=Zloc[threadIdx.x];
  }
  __syncthreads();
  if (threadIdx.x<NBINS) {
    atomicAdd(&Zglobal[threadIdx.x],Zloc2[threadIdx.x]);
  }
}
*/
__global__
void profilekernel(struct_plmd plmd,real* lambda,real* inweight,real* weightprofile,real* outweight,real* Zprofile)
{
  int b=blockIdx.x*blockDim.x+threadIdx.x;
  int s1,s2;
  int i1,i2;
  int k;
  real q1,q2;
  real w, wout;
  int itmp,jtmp;
  real Ztmp;
  // __shared__ int iloc[BLOCK];
  // __shared__ real Zloc[BLOCK];
  // __shared__ real Zloc2[NBINS+1];
  int NBINS=plmd.NBINS;
  int NBINS2=plmd.NBINS2;

  lambda+=plmd.nblocks*b;

  wout=0;

  if (b<plmd.B) {
    w=inweight[b];
  } else {
    w=0;
  }

  k=0;
  for (s1=0; s1<plmd.nsites; s1++) {
    for (s2=s1; s2<plmd.nsites; s2++) {
      if (s1==s2) { // Same site
        for (i1=plmd.block0_d[s1]; i1<plmd.block0_d[s1+1]; i1++) {
          __syncthreads();
          itmp=NBINS;
          Ztmp=w;
          if (b<plmd.B) {
            q1=lambda[i1];
            itmp=(int)floor(q1*NBINS);
            itmp=(itmp<0)?0:itmp; // sanity check for bad cuda floor
            itmp=(itmp>=NBINS)?(NBINS-1):itmp; // sanity check for q1=1
            // if (weightprofile) assert(w<=plmd.Zprofile_d[k*NBINS+itmp]);
            if (weightprofile) Ztmp*=weightprofile[k*NBINS+itmp];
          }
          if (outweight) wout+=Ztmp;
          // if (Zprofile) reduceBitonicSort(itmp,Ztmp,iloc,Zloc,Zloc2,&Zprofile[k*NBINS]);
          if (Zprofile) atomicAdd(&Zprofile[k*NBINS+itmp],Ztmp);
          k++;
        }

        for (i1=plmd.block0_d[s1]; i1<plmd.block0_d[s1+1]; i1++) {
          for (i2=i1+1; i2<plmd.block0_d[s2+1]; i2++) {
            __syncthreads();
            itmp=NBINS;
            Ztmp=w;
            if (b<plmd.B) {
              q1=lambda[i1];
              q2=lambda[i2];
              if (q1+q2>0.8) {
                itmp=(int)floor(q1/(q1+q2)*NBINS);
                itmp=(itmp<0)?0:itmp; // sanity check for bad cuda floor
                itmp=(itmp>=NBINS)?(NBINS-1):itmp; // sanity check for q1=1
                // if (weightprofile) assert(w<=plmd.Zprofile_d[k*NBINS+itmp]);
                if (weightprofile) Ztmp*=weightprofile[k*NBINS+itmp];
              } else { // WORKING - testing the next line
                Ztmp*=0;
              }
            }
            if (outweight) wout+=Ztmp;
            // if (Zprofile) reduceBitonicSort(itmp,Ztmp,iloc,Zloc,Zloc2,&Zprofile[k*NBINS]);
            if (Zprofile) atomicAdd(&Zprofile[k*NBINS+itmp],Ztmp);
            k++;
          }
        }

        if (plmd.block0_d[s1+1]-plmd.block0_d[s1]>2) {
          for (i1=plmd.block0_d[s1]; i1<plmd.block0_d[s1+1]; i1++) {
            for (i2=i1+1; i2<plmd.block0_d[s2+1]; i2++) {
              __syncthreads();
              itmp=NBINS;
              Ztmp=w;
              if (b<plmd.B) {
                q1=lambda[i1];
                q2=lambda[i2];
                // itmp=NBINS2*((int)floor(q1*NBINS2))+(int)floor(q2*NBINS2);
                itmp=(int)floor(q1*NBINS2);
                itmp=(itmp<0)?0:itmp; // sanity check for bad cuda floor
                itmp=(itmp>=NBINS2)?(NBINS2-1):itmp; // sanity check for q1=1
                jtmp=(int)floor(q2*NBINS2);
                jtmp=(jtmp<0)?0:jtmp; // sanity check for bad cuda floor
                jtmp=(itmp+jtmp>=NBINS2)?(NBINS2-1-itmp):jtmp; // sanity check for q1+q2=1
                itmp=NBINS2*itmp+jtmp;
                // if (weightprofile) assert(w<=plmd.Zprofile_d[k*NBINS+itmp]);
                if (weightprofile) Ztmp*=weightprofile[k*NBINS+itmp];
              }
              if (outweight) wout+=Ztmp;
              // if (Zprofile) reduceBitonicSort(itmp,Ztmp,iloc,Zloc,Zloc2,&Zprofile[k*NBINS]);
              if (Zprofile) atomicAdd(&Zprofile[k*NBINS+itmp],Ztmp);
              k++;
            }
          }
        }
      } else { // only needed if plmd.msprof
        for (i1=plmd.block0_d[s1]; i1<plmd.block0_d[s1+1]; i1++) {
          for (i2=plmd.block0_d[s2]; i2<plmd.block0_d[s2+1]; i2++) {
            __syncthreads();
            itmp=NBINS;
            Ztmp=w;
            if (b<plmd.B) {
              q1=lambda[i1];
              q2=lambda[i2];
              // itmp=NBINS2*((int)floor(q1*NBINS2))+(int)floor(q2*NBINS2);
              itmp=(int)floor(q1*NBINS2);
              itmp=(itmp<0)?0:itmp; // sanity check for bad cuda floor
              itmp=(itmp>=NBINS2)?(NBINS2-1):itmp; // sanity check for q1=1
              jtmp=(int)floor(q2*NBINS2);
              jtmp=(jtmp<0)?0:jtmp; // sanity check for bad cuda floor
              jtmp=(jtmp>=NBINS2)?(NBINS2-1):jtmp; // sanity check for q2=1
              itmp=NBINS2*itmp+jtmp;
              // if (weightprofile) assert(w<=plmd.Zprofile_d[k*NBINS+itmp]);
              if (weightprofile) Ztmp*=weightprofile[k*NBINS+itmp];
            }
            if (outweight) wout+=Ztmp;
            // if (Zprofile) reduceBitonicSort(itmp,Ztmp,iloc,Zloc,Zloc2,&Zprofile[k*NBINS]);
            if (Zprofile) atomicAdd(&Zprofile[k*NBINS+itmp],Ztmp);
            k++;
          }
        }
      }
    }
  }

  if (outweight) {
    if (b<plmd.B) {
      outweight[b]=wout;
    }
  }
}

void evaluateGimp(struct_plmd *plmd)
{
    int s1,s2,i,j,k,kn;
    int NBINS, NBINS2;
    real *mc_Zprofile, *Gimp;
    char fnm[MAXSTRING];
    FILE *fp;

    NBINS=plmd->NBINS;
    NBINS2=plmd->NBINS2;

    cudaMemset(plmd->mc_Zprofile_d,0,plmd->nprof*NBINS*sizeof(real));
    profilekernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0],plmd->mc_lambda_d,plmd->mc_ensweight_d,NULL,NULL,plmd->mc_Zprofile_d);

    mc_Zprofile=(real*)calloc(NBINS*plmd->nprof,sizeof(real));
    Gimp=(real*)calloc(NBINS,sizeof(real));
    cudaMemcpy(mc_Zprofile,plmd->mc_Zprofile_d,NBINS*plmd->nprof*sizeof(real),cudaMemcpyDeviceToHost);
    k=0;
    for (s1=0; s1<plmd->nsites; s1++) {
      for (s2=s1; s2<plmd->nsites; s2++) {
        if (s1==s2) { // Same site
          sprintf(fnm,"../G_imp/G1_%d.dat",plmd->nsubs[s1]);
          fp=fopen(fnm,"w");
          kn=k+plmd->nsubs[s1];
          for (i=0; i<NBINS; i++) {
            Gimp[i]=0;
            for (j=k; j<kn; j++) {
              Gimp[i]+=mc_Zprofile[NBINS*j+i];
            }
            if (Gimp[i]==0) {
              fprintf(stdout,"Warning, empty 1D Gimp[%d]\n",i);
            }
            Gimp[i]=-log(Gimp[i]);
            fprintf(fp,"%lg\n",Gimp[i]);
          }
          // for (j=k; j<kn; j++) {
          //   cudaMemcpy(&plmd->Gimp_d[NBINS*j],Gimp,NBINS*sizeof(real),cudaMemcpyHostToDevice);
          // }
          k=kn;
          fclose(fp);

          sprintf(fnm,"../G_imp/G12_%d.dat",plmd->nsubs[s1]);
          fp=fopen(fnm,"w");
          kn=k+(plmd->nsubs[s1]*(plmd->nsubs[s1]-1))/2;
          for (i=0; i<NBINS; i++) {
            Gimp[i]=0;
            for (j=k; j<kn; j++) {
              Gimp[i]+=mc_Zprofile[NBINS*j+i];
            }
            if (Gimp[i]==0) {
              fprintf(stdout,"Warning, empty transition Gimp[%d]\n",i);
            }
            Gimp[i]=-log(Gimp[i]);
            fprintf(fp,"%lg\n",Gimp[i]);
          }
          // for (j=k; j<kn; j++) {
          //   cudaMemcpy(&plmd->Gimp_d[NBINS*j],Gimp,NBINS*sizeof(real),cudaMemcpyHostToDevice);
          // }
          k=kn;
          fclose(fp);

          if (plmd->nsubs[s1]>2) {
            sprintf(fnm,"../G_imp/G2_%d.dat",plmd->nsubs[s1]);
            fp=fopen(fnm,"w");
            kn=k+(plmd->nsubs[s1]*(plmd->nsubs[s1]-1))/2;
            for (i=0; i<NBINS; i++) {
              Gimp[i]=0;
              for (j=k; j<kn; j++) {
                Gimp[i]+=mc_Zprofile[NBINS*j+i];
              }
              fprintf(stdout,"i=%d i/%d=%d i%%%d=%d Gimp[%d]=%lf\n",i,NBINS2,i/NBINS2,NBINS2,i%NBINS2,i,Gimp[i]);
              if (Gimp[i]==0 && (i/NBINS2)+(i%NBINS2)<NBINS2) {
                fprintf(stdout,"Warning, empty 2D intrasite Gimp[%d]\n",i);
              } else if (Gimp[i]!=0 && (i/NBINS2)+(i%NBINS2)>=NBINS2) {
                fprintf(stdout,"Warning, nonempty 2D intrasite Gimp[%d] should be empty\n",i);
              }
              Gimp[i]=-log(Gimp[i]);
              fprintf(fp,"%lg\n",Gimp[i]);
            }
            // for (j=k; j<kn; j++) {
            //   cudaMemcpy(&plmd->Gimp_d[NBINS*j],Gimp,NBINS*sizeof(real),cudaMemcpyHostToDevice);
            // }
            k=kn;
            fclose(fp);
          }
        } else { // only needed if plmd->msprof
          sprintf(fnm,"../G_imp/G1_%d_%d.dat",plmd->nsubs[s1],plmd->nsubs[s2]);
          fp=fopen(fnm,"w");
          kn=k+plmd->nsubs[s1]*plmd->nsubs[s2];
          for (i=0; i<NBINS; i++) {
            Gimp[i]=0;
            for (j=k; j<kn; j++) {
              Gimp[i]+=mc_Zprofile[NBINS*j+i];
            }
            if (Gimp[i]==0) {
              fprintf(stdout,"Warning, empty 2D intersite Gimp[%d]\n",i);
            }
            Gimp[i]=-log(Gimp[i]);
            fprintf(fp,"%lg\n",Gimp[i]);
          }
          // for (j=k; j<kn; j++) {
          //   cudaMemcpy(&plmd->Gimp_d[NBINS*j],Gimp,NBINS*sizeof(real),cudaMemcpyHostToDevice);
          // }
          k=kn;
          fclose(fp);
        }
      }
    }
    free(mc_Zprofile);
    free(Gimp);
}
#else

void profile1(real *Gimp,int NBINS,int nblocks,int B,int i1,real *lambda)
{
  int b, itmp;
  real q1;
  for (b=0; b<B; b++) {
    q1=lambda[b*nblocks+i1];
    itmp=(int)floor(q1*NBINS);
    itmp=(itmp<0)?0:itmp; // sanity check for bad cuda floor
    itmp=(itmp>=NBINS)?(NBINS-1):itmp; // sanity check for q1=1
    Gimp[itmp]+=1;
  }
}

void profile12(real *Gimp,int NBINS,int nblocks,int B,int i1,int i2,real *lambda)
{
  int b, itmp;
  real q1, q2;
  for (b=0; b<B; b++) {
    q1=lambda[b*nblocks+i1];
    q2=lambda[b*nblocks+i2];
    if (q1+q2>0.8) {
      itmp=(int)floor(q1/(q1+q2)*NBINS);
      itmp=(itmp<0)?0:itmp; // sanity check for bad cuda floor
      itmp=(itmp>=NBINS)?(NBINS-1):itmp; // sanity check for q1=1
      Gimp[itmp]+=1;
    }
  }
}

void profile2(real *Gimp,int NBINS2,int nblocks,int B,int i1,int i2,real *lambda)
{
  int b, itmp, jtmp;
  real q1, q2;
  for (b=0; b<B; b++) {
    q1=lambda[b*nblocks+i1];
    q2=lambda[b*nblocks+i2];
    itmp=(int)floor(q1*NBINS2);
    itmp=(itmp<0)?0:itmp; // sanity check for bad cuda floor
    itmp=(itmp>=NBINS2)?(NBINS2-1):itmp; // sanity check for q1=1
    jtmp=(int)floor(q2*NBINS2);
    jtmp=(jtmp<0)?0:jtmp; // sanity check for bad cuda floor
    jtmp=(itmp+jtmp>=NBINS2)?(NBINS2-1-itmp):jtmp; // sanity check for q1+q2=1
    itmp=NBINS2*itmp+jtmp;
    Gimp[itmp]+=1;
  }
}

void profile11(real *Gimp,int NBINS2,int nblocks,int B,int i1,int i2,real *lambda)
{
  int b, itmp, jtmp;
  real q1, q2;
  for (b=0; b<B; b++) {
    q1=lambda[b*nblocks+i1];
    q2=lambda[b*nblocks+i2];
    itmp=(int)floor(q1*NBINS2);
    itmp=(itmp<0)?0:itmp; // sanity check for bad cuda floor
    itmp=(itmp>=NBINS2)?(NBINS2-1):itmp; // sanity check for q1=1
    jtmp=(int)floor(q2*NBINS2);
    jtmp=(jtmp<0)?0:jtmp; // sanity check for bad cuda floor
    jtmp=(jtmp>=NBINS2)?(NBINS2-1):jtmp; // sanity check for q2=1
    itmp=NBINS2*itmp+jtmp;
  }
}

void evaluateGimp(struct_plmd *plmd)
{
  int s1,s2,i1,i2,i;
  real *Gimp;
  char fnm[MAXSTRING];
  FILE *fp;

  Gimp=(real*)calloc(plmd->NBINS,sizeof(real));
  for (s1=0; s1<plmd->nsites; s1++) {
    for (s2=s1; s2<plmd->nsites; s2++) {
      if (s1==s2) { // Same site
        for (i=0; i<plmd->NBINS; i++) {
          Gimp[i]=0;
        }
        for (i1=plmd->block0[s1]; i1<plmd->block0[s1+1]; i1++) {
          profile1(Gimp,plmd->NBINS,plmd->nblocks,plmd->B,i1,plmd->mc_lambda);
        }
        sprintf(fnm,"../G_imp/G1_%d.dat",plmd->nsubs[s1]);
        fp=fopen(fnm,"w");
        if (fp==NULL) {
          fprintf(stderr,"Error opening file %s\n",fnm);
        }
        for (i=0; i<plmd->NBINS; i++) {
          if (Gimp[i]==0) {
            fprintf(stdout,"Warning, empty 1D Gimp[%d]\n",i);
          }
          fprintf(fp,"%lg\n",-log(Gimp[i]));
        }
        fclose(fp);

        for (i=0; i<plmd->NBINS; i++) {
          Gimp[i]=0;
        }
        for (i1=plmd->block0[s1]; i1<plmd->block0[s1+1]; i1++) {
          for (i2=i1+1; i2<plmd->block0[s1+1]; i2++) {
            profile12(Gimp,plmd->NBINS,plmd->nblocks,plmd->B,i1,i2,plmd->mc_lambda);
            profile12(Gimp,plmd->NBINS,plmd->nblocks,plmd->B,i2,i1,plmd->mc_lambda);
          }
        }
        sprintf(fnm,"../G_imp/G12_%d.dat",plmd->nsubs[s1]);
        fp=fopen(fnm,"w");
        for (i=0; i<plmd->NBINS; i++) {
          if (Gimp[i]==0) {
            fprintf(stdout,"Warning, empty transition Gimp[%d]\n",i);
          }
          fprintf(fp,"%lg\n",-log(Gimp[i]));
        }
        fclose(fp);

        if (plmd->nsubs[s1]>2) {
          for (i=0; i<plmd->NBINS; i++) {
            Gimp[i]=0;
          }
          for (i1=plmd->block0[s1]; i1<plmd->block0[s1+1]; i1++) {
            for (i2=i1+1; i2<plmd->block0[s1+1]; i2++) {
              profile2(Gimp,plmd->NBINS2,plmd->nblocks,plmd->B,i1,i2,plmd->mc_lambda);
              profile2(Gimp,plmd->NBINS2,plmd->nblocks,plmd->B,i2,i1,plmd->mc_lambda);
            }
          }
          sprintf(fnm,"../G_imp/G2_%d.dat",plmd->nsubs[s1]);
          fp=fopen(fnm,"w");
          for (i=0; i<plmd->NBINS; i++) {
            if (Gimp[i]==0 && (i/plmd->NBINS2)+(i%plmd->NBINS2)<plmd->NBINS2) {
              fprintf(stdout,"Warning, empty 2D intrasite Gimp[%d]\n",i);
            } else if (Gimp[i]!=0 && (i/plmd->NBINS2)+(i%plmd->NBINS2)>=plmd->NBINS2) {
              fprintf(stdout,"Warning, nonempty 2D intrasite Gimp[%d] should be empty\n",i);
            }
            fprintf(fp,"%lg\n",-log(Gimp[i]));
          }
          fclose(fp);
        }
      } else { // only needed if (plmd->msprof)
        for (i=0; i<plmd->NBINS; i++) {
          Gimp[i]=0;
        }
        for (i1=plmd->block0[s1]; i1<plmd->block0[s1+1]; i1++) {
          for (i2=plmd->block0[s2]; i2<plmd->block0[s2+1]; i2++) {
            profile11(Gimp,plmd->NBINS2,plmd->nblocks,plmd->B,i1,i2,plmd->mc_lambda);
          }
        }
        sprintf(fnm,"../G_imp/G1_%d_%d.dat",plmd->nsubs[s1],plmd->nsubs[s2]);
        fp=fopen(fnm,"w");
        for (i=0; i<plmd->NBINS; i++) {
          if (Gimp[i]==0) {
            fprintf(stdout,"Warning, empty 2D intersite Gimp[%d]\n",i);
          }
          fprintf(fp,"%lg\n",-log(Gimp[i]));
        }
        fclose(fp);
      }
    }
  }
  free(Gimp);
}
#endif

void printSamples(struct_plmd *plmd)
{
  int i,j;
  FILE *fp;

  fp=fopen("mc_Lambda.dat","w");
  for (i=0; i<plmd->B; i++) {
    for (j=0; j<plmd->nblocks; j++) {
      fprintf(fp," %lg",(double) plmd->mc_lambda[i*plmd->nblocks+j]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}

void run(struct_plmd *plmd)
{
  monte_carlo_Z(plmd);
#ifdef USEGPU
  cudaMemcpy(plmd->mc_lambda_d,plmd->mc_lambda,plmd->B*plmd->nblocks*sizeof(real),cudaMemcpyHostToDevice);
#endif

  if (plmd->mode==gimp) {
    evaluateGimp(plmd);
  } else if (plmd->mode==sample) {
    printSamples(plmd);
  }
}

void finish(struct_plmd *plmd,int argc, char *argv[])
{
  free(plmd->mc_lambda);

  free(plmd->nsubs);
  free(plmd->block0);

  free(plmd);
}

int main(int argc, char *argv[])
{
  struct_plmd *plmd;

  plmd = setup(argc,argv);
 
  run(plmd);

  finish(plmd,argc,argv);

  return 0;
}
