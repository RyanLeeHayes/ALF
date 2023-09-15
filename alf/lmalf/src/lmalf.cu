// Written by Ryan Hayes 2017-06-20
// plmDCA algorithm from R470 - DOI: 10.1016/j.jcp.2014.07.024
// Quasi newton equations from https://www.rose-hulman.edu/~bryan/lottamath/quasinewton.pdf

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <omp.h>

// #include <assert.h>

//   ID=omp_get_thread_num();
//   NID=omp_get_max_threads();

#include "lmalf.h"

#define MAXLENGTH 1024
#define BLOCK 512
#define kB 0.00198614L

#define NBINS 256
#define NBINS2 16

#define PROFILE true
#define MOMENT true

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

void printarray(real* x,int N,char* fnm)
{
  FILE *fp;
  real *xcpu;
  int i;

  fp=fopen(fnm,"w");
  xcpu=(real*)malloc(N*sizeof(real));
  cudaMemcpy(xcpu,x,N*sizeof(real),cudaMemcpyDeviceToHost);
  for (i=0; i<N; i++) {
    // fprintf(fp,"%g\n",xcpu[i]);
    fprintf(fp,"%22.15e\n",xcpu[i]);
  }
  free(xcpu);
  fclose(fp);
}

double randDouble()
{
  return (rand()+0.5)/(RAND_MAX+1.0);
}

void monte_carlo_Z(struct_plmd plmd)
{
  int ibeg,iend,Ns;
  int Neq=plmd.B/10;
  int Nmc=plmd.B;
  real *theta;
  int s,i,j;
  real b, st, norm;
  real thetaNew,eOld,eNew;

  theta=(real*) calloc(plmd.nblocks,sizeof(real));

  for (s=0; s<plmd.nsites; s++) {
    ibeg=plmd.block0[s];
    iend=plmd.block0[s+1];
    Ns=iend-ibeg;

    b=1;
    for (i=0; i<50; i++) {
      b=0.5*log(0.25*b*Ns*Ns*M_PI/2);
      if (!(b>0)) b=0;
    }

    theta[ibeg]=M_PI/2;
    for (i=ibeg+1; i<iend; i++) {
      theta[i]=3*M_PI/2;
    }

    for (i=-Neq; i<Nmc; i++) {
      if (i%Neq==0) {
        fprintf(stdout,"Partition Function Sample Step %d\n",i);
      }

      for (j=ibeg; j<iend; j++) {
        st=(-0.5*sin(theta[j])+0.5);
        eOld=-b*st*st*st*st;

        thetaNew=2*M_PI*randDouble();
        st=(-0.5*sin(thetaNew)+0.5);
        eNew=-b*st*st*st*st;

        if (exp(eOld-eNew)>randDouble()) {
          theta[j]=thetaNew;
        }
      }

      if (i>=0) {
        norm=0;
        for (j=ibeg; j<iend; j++) {
          norm+=exp(5.5*sin(theta[j]));
        }
        for (j=ibeg; j<iend; j++) {
          plmd.mc_lambda[plmd.nblocks*i+j]=exp(5.5*sin(theta[j]))/norm;
        }
      }
    }
  }

  free(theta);
}

struct_plmd* setup(int argc, char *argv[])
{
  struct_plmd *plmd;
  int si,sj,i,j,k,l;
  real k0;
  FILE *fp;
  char line[MAXLENGTH];

  if (argc<7) {
    fprintf(stderr,"Error: not enough input arguments\n");
    exit(1);
  }

  plmd=(struct_plmd*) malloc(sizeof(struct_plmd));

  fp=fopen("nsubs","r");
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
  plmd->block2site=(int*) calloc(plmd->nblocks,sizeof(int));
  k=0;
  for(i=0; i<plmd->nsites; i++) {
    plmd->block0[i]=k;
    for(j=0; j<plmd->nsubs[i]; j++) {
      plmd->block2site[k]=i;
      k++;
    }
  }
  plmd->block0[i]=k;

  cudaMalloc(&(plmd->block0_d),(plmd->nsites+1)*sizeof(int));
  cudaMemcpy(plmd->block0_d,plmd->block0,(plmd->nsites+1)*sizeof(int),cudaMemcpyHostToDevice);

  i=sscanf(argv[2],"%d",&plmd->ms);
  if (i!=1) {
    fprintf(stderr,"Error, first argument must be a boolean flag for whether to use multisite coupling\n");
    exit(1);
  }

  i=sscanf(argv[3],"%d",&plmd->msprof);
  if (i!=1) {
    fprintf(stderr,"Error, second argument should indicate whether to use multisite profiles.\n");
    exit(1);
  }

  fp=fopen(argv[4],"r");
  for (plmd->B=0; fgets(line,MAXLENGTH,fp) != NULL; plmd->B++) {
    ;
  }
  fprintf(stdout,"%d frames\n",plmd->B); // DEBUG
  fclose(fp);

  plmd->lambda=(real*) calloc(plmd->B*plmd->nblocks,sizeof(real));
  plmd->ensweight=(real*) calloc(plmd->B,sizeof(real));
  plmd->mc_lambda=(real*) calloc(plmd->B*plmd->nblocks,sizeof(real));
  plmd->mc_ensweight=(real*) calloc(plmd->B,sizeof(real));

  fp=fopen(argv[4],"r");
  for (i=0;i<plmd->B;i++) {
    for (j=0;j<plmd->nblocks;j++) {
      double buffer;
      fscanf(fp,"%lf",&buffer);
      plmd->lambda[i*plmd->nblocks+j]=buffer;
    }
  }
  fclose(fp);

  monte_carlo_Z(plmd[0]);

  fp=fopen(argv[5],"r");
  for (i=0; i<plmd->B; i++) {
    double buffer;
    fscanf(fp,"%lf",&buffer);
    plmd->ensweight[i]=buffer;
    plmd->mc_ensweight[i]=1;
  }
  fclose(fp);

  if (argc>=8) {
    double criteria;
    i=sscanf(argv[7],"%lg",&criteria);
    fprintf(stdout,"Note, found seventh argument indicating halting criteria. Overwriting default value of 1.25e-3\n");
    if (i!=1) {
      fprintf(stderr,"Error, seventh argument should indicate halting criteria.\n");
      exit(1);
    }
    plmd->criteria=criteria;
  } else {
    plmd->critetia=1.25e-3;
  }

  cudaMalloc(&(plmd->lambda_d),plmd->B*plmd->nblocks*sizeof(real));
  cudaMalloc(&(plmd->mc_lambda_d),plmd->B*plmd->nblocks*sizeof(real));
  cudaMalloc(&(plmd->ensweight_d),plmd->B*sizeof(real));
  cudaMalloc(&(plmd->mc_ensweight_d),plmd->B*sizeof(real));

  cudaMemcpy(plmd->lambda_d,plmd->lambda,plmd->B*plmd->nblocks*sizeof(real),cudaMemcpyHostToDevice);
  cudaMemcpy(plmd->mc_lambda_d,plmd->mc_lambda,plmd->B*plmd->nblocks*sizeof(real),cudaMemcpyHostToDevice);
  cudaMemcpy(plmd->ensweight_d,plmd->ensweight,plmd->B*sizeof(real),cudaMemcpyHostToDevice);
  cudaMemcpy(plmd->mc_ensweight_d,plmd->mc_ensweight,plmd->B*sizeof(real),cudaMemcpyHostToDevice);



  // count nbias
  plmd->nbias=0;
  for (i=0; i<plmd->nsites; i++) {
    for (j=i; j<plmd->nsites; j++) {
      if (i==j) {
        plmd->nbias+=plmd->nsubs[i]+(5*plmd->nsubs[i]*(plmd->nsubs[i]-1))/2;
      } else if (plmd->ms==1) {
        plmd->nbias+=5*plmd->nsubs[i]*plmd->nsubs[j];
      } else if (plmd->ms==2) {
        plmd->nbias+=plmd->nsubs[i]*plmd->nsubs[j];
      }
    }
  }

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
      } else if (plmd->msprof) {
        plmd->nprof+=plmd->nsubs[i]*plmd->nsubs[j];
      }
    }
  }

  // plmd->nx=plmd->nbias+plmd->nprof;
  plmd->nx=plmd->nbias;

  double temperature;
  sscanf(argv[1],"%lg",&temperature);
  plmd->kT=kB*temperature;

  // Set regularization constants
  real kp=1.0/(plmd->kT*plmd->kT);
  plmd->kx=(real*) calloc(plmd->nx,sizeof(real));
  plmd->xr=(real*) calloc(plmd->nx,sizeof(real));
  real *xr_x, *xr_s;
  // load starting values if needed for ms==1
  if (plmd->ms==1) {
    xr_x=(real*) calloc(plmd->nblocks*plmd->nblocks,sizeof(real));
    fp=fopen("x_prev.dat","r");
    for(i=0; i<plmd->nblocks; i++) {
      for(j=0; j<plmd->nblocks; j++) {
        double buffer;
        fscanf(fp,"%lf",&buffer);
        xr_x[i*plmd->nblocks+j]=buffer;
      }
    }
    fclose(fp);
    xr_s=(real*) calloc(plmd->nblocks*plmd->nblocks,sizeof(real));
    fp=fopen("s_prev.dat","r");
    for(i=0; i<plmd->nblocks; i++) {
      for(j=0; j<plmd->nblocks; j++) {
        double buffer;
        fscanf(fp,"%lf",&buffer);
        xr_s[i*plmd->nblocks+j]=buffer;
      }
    }
    fclose(fp);
  }
  // k0=1e-2; // 1.0/400;
  k0=kp/400;
  // k0=1;
  k=0;
  for (si=0; si<plmd->nsites; si++) {
    for (sj=si; sj<plmd->nsites; sj++) {
      if (si==sj) {
        for (i=0; i<plmd->nsubs[si]; i++) {
          plmd->kx[k++]=k0/4; // b
          for (j=i+1; j<plmd->nsubs[sj]; j++) {
            plmd->kx[k++]=k0/64; // c
            plmd->kx[k++]=k0/4; // x
            plmd->kx[k++]=k0/4; // x
            plmd->kx[k++]=k0/1; // s
            plmd->kx[k++]=k0/1; // s
          }
        }
      } else if (plmd->ms) {
        for (i=0; i<plmd->nsubs[si]; i++) {
          for (j=0; j<plmd->nsubs[sj]; j++) {
            plmd->kx[k++]=k0/4; // c
            if (plmd->ms==1) {
              plmd->xr[k]=xr_x[(plmd->block0[si]+i)*plmd->nblocks+plmd->block0[sj]+j]; // x
              plmd->kx[k++]=k0/0.25; // x
              plmd->xr[k]=xr_x[(plmd->block0[sj]+j)*plmd->nblocks+plmd->block0[si]+i]; // x
              plmd->kx[k++]=k0/0.25; // x
              plmd->xr[k]=xr_s[(plmd->block0[si]+i)*plmd->nblocks+plmd->block0[sj]+j]; // s
              plmd->kx[k++]=k0/0.25; // s
              plmd->xr[k]=xr_s[(plmd->block0[sj]+j)*plmd->nblocks+plmd->block0[si]+i]; // s
              plmd->kx[k++]=k0/0.25; // s
            }
          }
        }
      }
    }
  }
  if (plmd->ms==1) {
    free(xr_x);
    free(xr_s);
  }
  // No restraints on average profile values - treated implicitly now
  /*for (i=0; i<plmd->nprof; i++) {
    plmd->kx[k++]=0;
  }*/
  cudaMalloc(&plmd->kx_d,plmd->nx*sizeof(real));
  cudaMemcpy(plmd->kx_d,plmd->kx,plmd->nx*sizeof(real),cudaMemcpyHostToDevice);
  cudaMalloc(&plmd->xr_d,plmd->nx*sizeof(real));
  cudaMemcpy(plmd->xr_d,plmd->xr,plmd->nx*sizeof(real),cudaMemcpyHostToDevice);

  // plmd->kprofile=1.0/NBINS;
  real kp0=kp/NBINS;
  plmd->kprofile=(real*)calloc(NBINS*plmd->nprof,sizeof(real));
  k=0;
  for (si=0; si<plmd->nsites; si++) {
    for (sj=si; sj<plmd->nsites; sj++) {
      if (si==sj) { // Same site
        for (i=plmd->block0[si]; i<plmd->block0[si+1]; i++) {
          for (l=0; l<NBINS; l++) {
            plmd->kprofile[NBINS*k+l]=kp0;
            if (l==NBINS-1) plmd->kprofile[NBINS*k+l]*=(NBINS/4.0);
          }
          k++;
        }
        for (i=plmd->block0[si]; i<plmd->block0[si+1]; i++) {
          for (j=i+1; j<plmd->block0[sj+1]; j++) {
            for (l=0; l<NBINS; l++) {
              plmd->kprofile[NBINS*k+l]=kp0/((plmd->nsubs[si]-1)/2.0);
            }
            k++;
          }
        }
        if (plmd->nsubs[si]>2) {
          for (i=plmd->block0[si]; i<plmd->block0[si+1]; i++) {
            for (j=i+1; j<plmd->block0[sj+1]; j++) {
              for (l=0; l<NBINS; l++) {
                plmd->kprofile[NBINS*k+l]=kp0/((plmd->nsubs[si]-1)/2.0);
              }
              k++;
            }
          }
        }
      } else if (plmd->msprof) {
        for (i=plmd->block0[si]; i<plmd->block0[si+1]; i++) {
          for (j=plmd->block0[sj]; j<plmd->block0[sj+1]; j++) {
            for (l=0; l<NBINS; l++) {
              plmd->kprofile[NBINS*k+l]=kp0/(plmd->nsubs[si]*plmd->nsubs[sj]);
              if (l==NBINS-1) plmd->kprofile[NBINS*k+l]*=(NBINS/4.0);
            }
            k++;
          }
        }
      }
    }
  }
  cudaMalloc(&plmd->kprofile_d,NBINS*plmd->nprof*sizeof(real));
  cudaMemcpy(plmd->kprofile_d,plmd->kprofile,NBINS*plmd->nprof*sizeof(real),cudaMemcpyHostToDevice);

  plmd->L=(real*) calloc(1,sizeof(real));
  plmd->dLds=(real*) calloc(1,sizeof(real));

  cudaMalloc(&plmd->L_d,sizeof(real));
  cudaMalloc(&plmd->dLds_d,sizeof(real));

  plmd->x=(real*) calloc(plmd->nx,sizeof(real));
  plmd->dLdx=(real*) calloc(plmd->nx,sizeof(real));
  cudaMalloc(&plmd->dLdx_d,plmd->nx*sizeof(real));
  // plmd->Hinv=(real*) calloc(plmd->nx*plmd->nx,sizeof(real));
  plmd->Nmemax=50;
  plmd->Nmem=0;
  plmd->d_x=(real*) calloc(plmd->nx*plmd->Nmemax,sizeof(real));
  plmd->d_dLdx=(real*) calloc(plmd->nx*plmd->Nmemax,sizeof(real));
  plmd->rho=(real*) calloc(plmd->Nmemax,sizeof(real));
  plmd->alpha=(real*) calloc(plmd->nx*plmd->Nmemax,sizeof(real));
  plmd->beta=(real*) calloc(plmd->nx*plmd->Nmemax,sizeof(real));

  plmd->hi=(real*) calloc(plmd->nx,sizeof(real));
  plmd->x0=(real*) calloc(plmd->nx,sizeof(real));
  plmd->dLdx0=(real*) calloc(plmd->nx,sizeof(real));



  cudaMalloc(&(plmd->E_d),plmd->B*sizeof(real));
  cudaMalloc(&(plmd->dEds_d),plmd->B*sizeof(real));
  cudaMalloc(&(plmd->mc_E_d),plmd->B*sizeof(real));
  cudaMalloc(&(plmd->mc_dEds_d),plmd->B*sizeof(real));

  cudaMalloc(&(plmd->weight_d),plmd->B*sizeof(real));
  cudaMalloc(&(plmd->mc_weight_d),plmd->B*sizeof(real));

  cudaMalloc(&(plmd->x_d),plmd->nx*sizeof(real));
  cudaMalloc(&(plmd->dxds_d),plmd->nx*sizeof(real));

  cudaMalloc(&(plmd->Z_d),plmd->nprof*NBINS*sizeof(real));
  cudaMalloc(&(plmd->mc_Z_d),plmd->nprof*NBINS*sizeof(real));
  cudaMalloc(&(plmd->Zprofile_d),plmd->nprof*NBINS*sizeof(real));
  cudaMalloc(&(plmd->mc_Zprofile_d),plmd->nprof*NBINS*sizeof(real));

  cudaMalloc(&(plmd->dLdZprofile_d),plmd->nprof*NBINS*sizeof(real));
  cudaMalloc(&(plmd->mc_dLdZprofile_d),plmd->nprof*NBINS*sizeof(real));
  cudaMalloc(&(plmd->dLdE_d),plmd->B*sizeof(real));
  cudaMalloc(&(plmd->Gimp_d),plmd->nprof*NBINS*sizeof(real));
  cudaMalloc(&(plmd->G_d),plmd->nprof*NBINS*sizeof(real));

  cudaMalloc(&(plmd->Esum_d),sizeof(real));
  cudaMalloc(&(plmd->dEdssum_d),sizeof(real));
  cudaMalloc(&(plmd->mc_dEdssum_d),sizeof(real));
  cudaMalloc(&(plmd->moments_d),plmd->nbias*sizeof(real));
  cudaMalloc(&(plmd->mc_moments_d),plmd->nbias*sizeof(real));
  cudaMalloc(&(plmd->sumensweight_d),sizeof(real));

  plmd->fplog=fopen("log.log","w");

  return plmd;
}

__device__
void reduce(real local,real* shared,real* global)
{
  int k;

  shared[threadIdx.x]=local;

  __syncthreads();

  for (k=1; k<BLOCK; k*=2) {
    if ((threadIdx.x % (2*k)) == 0) {
      shared[threadIdx.x]+=shared[threadIdx.x+k];
    }
    __syncthreads();
  }

  if (threadIdx.x==0) {
    atomicAdd(global,shared[0]);
  }
}

__device__
void reduceNBINS(real local,real* shared,real* global)
{
  int k;

  shared[threadIdx.x]=local;

  __syncthreads();

  for (k=1; k<NBINS; k*=2) {
    if ((threadIdx.x % (2*k)) == 0) {
      shared[threadIdx.x]+=shared[threadIdx.x+k];
    }
    __syncthreads();
  }

  if (threadIdx.x==0) {
    atomicAdd(global,shared[0]);
  }
}

__device__
void reduceBroadcast(real local,real* shared)
{
  int k;
  real buf;

  shared[threadIdx.x]=local;

  __syncthreads();

  for (k=1; k<NBINS; k*=2) {
    buf=0;
    if ((threadIdx.x^k)<NBINS) buf=shared[threadIdx.x^k];
    __syncthreads();
    shared[threadIdx.x]+=buf;
    __syncthreads();
  }
}

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
/*
{
  if (itmp<NBINS) {
    atomicAdd(&Zglobal[itmp],Ztmp);
  }
}*/

__global__
void energykernel(struct_plmd plmd,real* x,real* lambda,real* energy)
{
  int b=blockIdx.x*blockDim.x+threadIdx.x;
  int s1,s2;
  int i1,i2;
  int k;
  real q1,q2;
  real E;

  lambda+=plmd.nblocks*b;

  if (b<plmd.B) {
    k=0;
    E=0;
    for (s1=0; s1<plmd.nsites; s1++) {
      for (s2=s1; s2<plmd.nsites; s2++) {
        if (s1==s2) { // Same site
          for (i1=plmd.block0_d[s1]; i1<plmd.block0_d[s1+1]; i1++) {
            q1=lambda[i1];
            E+=x[k]*q1;
            k++;
            for (i2=i1+1; i2<plmd.block0_d[s1+1]; i2++) {
              q2=lambda[i2];
              E+=x[k]*q1*q2;
              k++;
              E+=x[k]*q2*(1-exp(-q1/0.18));
              k++;
              E+=x[k]*q1*(1-exp(-q2/0.18));
              k++;
              E+=x[k]*q2*(1-1/(q1/0.017+1));
              k++;
              E+=x[k]*q1*(1-1/(q2/0.017+1));
              k++;
            }
          }
        } else if (plmd.ms) { // Different sites
          for (i1=plmd.block0_d[s1]; i1<plmd.block0_d[s1+1]; i1++) {
            q1=lambda[i1];
            for (i2=plmd.block0_d[s2]; i2<plmd.block0_d[s2+1]; i2++) {
              q2=lambda[i2];
              E+=x[k]*q1*q2;
              k++;
              if (plmd.ms==1) { // include extra terms
                E+=x[k]*q2*(1-exp(-q1/0.18));
                k++;
                E+=x[k]*q1*(1-exp(-q2/0.18));
                k++;
                E+=x[k]*q2*(1-1/(q1/0.017+1));
                k++;
                E+=x[k]*q1*(1-1/(q2/0.017+1));
                k++;
              }
            }
          }
        }
      }
    }
    energy[b]=E;
  }
}

__global__
void dotenergykernel(struct_plmd plmd,real sign,real* x,real* y,real* z)
{
  int b=blockIdx.x*blockDim.x+threadIdx.x;
  real xtmp;
  __shared__ real xloc[BLOCK];

  if (b<plmd.B) {
    xtmp=sign*x[b]*y[b];
  } else {
    xtmp=0;
  }
  reduce(xtmp,xloc,z);
}

__global__
void weightedenergykernel(struct_plmd plmd,real sign,real* lambda,real* weight,real* dEdx)
{
  int b=blockIdx.x*blockDim.x+threadIdx.x;
  int s1,s2;
  int i1,i2;
  int k;
  real q1,q2;
  real w,E;
  __shared__ real Eloc[BLOCK];

  lambda+=plmd.nblocks*b;

  w=0;
  q1=0;
  q2=0;

  if (b<plmd.B) {
    w=sign*weight[b];
  }

  k=0;
  for (s1=0; s1<plmd.nsites; s1++) {
    for (s2=s1; s2<plmd.nsites; s2++) {
      if (s1==s2) { // Same site
        for (i1=plmd.block0_d[s1]; i1<plmd.block0_d[s1+1]; i1++) {
          if (b<plmd.B) q1=lambda[i1];
          E=w*q1;
          reduce(E,Eloc,&dEdx[k]);
          k++;
          for (i2=i1+1; i2<plmd.block0_d[s1+1]; i2++) {
            if (b<plmd.B) q2=lambda[i2];
            E=w*q1*q2;
            reduce(E,Eloc,&dEdx[k]);
            k++;
            E=w*q2*(1-exp(-q1/0.18));
            reduce(E,Eloc,&dEdx[k]);
            k++;
            E=w*q1*(1-exp(-q2/0.18));
            reduce(E,Eloc,&dEdx[k]);
            k++;
            E=w*q2*(1-1/(q1/0.017+1));
            reduce(E,Eloc,&dEdx[k]);
            k++;
            E=w*q1*(1-1/(q2/0.017+1));
            reduce(E,Eloc,&dEdx[k]);
            k++;
          }
        }
      } else if (plmd.ms) { // Different sites
        for (i1=plmd.block0_d[s1]; i1<plmd.block0_d[s1+1]; i1++) {
          q1=lambda[i1];
          for (i2=plmd.block0_d[s2]; i2<plmd.block0_d[s2+1]; i2++) {
            q2=lambda[i2];
            E=w*q1*q2;
            reduce(E,Eloc,&dEdx[k]);
            k++;
            if (plmd.ms==1) { // include extra terms
              E=w*q2*(1-exp(-q1/0.18));
              reduce(E,Eloc,&dEdx[k]);
              k++;
              E=w*q1*(1-exp(-q2/0.18));
              reduce(E,Eloc,&dEdx[k]);
              k++;
              E=w*q2*(1-1/(q1/0.017+1));
              reduce(E,Eloc,&dEdx[k]);
              k++;
              E=w*q1*(1-1/(q2/0.017+1));
              reduce(E,Eloc,&dEdx[k]);
              k++;
            }
          }
        }
      }
    }
  }
}

__global__
void boltzmannkernel(struct_plmd plmd,real sign,real* energy,real s,real* denergyds,real* inweight,real* outweight,real* Z)
{
  int b=blockIdx.x*blockDim.x+threadIdx.x;
  real E;
  real w;
  __shared__ real Zloc[BLOCK];

  if (b<plmd.B) {
    w=inweight[b];
    E=energy[b];
    if (s) { // add the displacement if it is non-zero
      E+=s*denergyds[b];
    }
    w*=exp(-sign*E/plmd.kT);
    outweight[b]=w;
  } else {
    w=0;
  }

  if (Z) { // calculate the partition function if requested
    __syncthreads();
    reduce(w,Zloc,Z);
  }
}

__global__
void profilekernel(struct_plmd plmd,real* lambda,real* inweight,real* weightprofile,real* outweight,real* Zprofile)
{
  int b=blockIdx.x*blockDim.x+threadIdx.x;
  int s1,s2;
  int i1,i2;
  int k;
  real q1,q2;
  real w, wout;
  int itmp;
  real Ztmp;
  __shared__ int iloc[BLOCK];
  __shared__ real Zloc[BLOCK];
  __shared__ real Zloc2[NBINS+1];

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
            // if (weightprofile) assert(w<=plmd.Zprofile_d[k*NBINS+itmp]);
            if (weightprofile) Ztmp*=weightprofile[k*NBINS+itmp];
          }
          if (outweight) wout+=Ztmp;
          if (Zprofile) reduceBitonicSort(itmp,Ztmp,iloc,Zloc,Zloc2,&Zprofile[k*NBINS]);
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
                // if (weightprofile) assert(w<=plmd.Zprofile_d[k*NBINS+itmp]);
                if (weightprofile) Ztmp*=weightprofile[k*NBINS+itmp];
              } else { // WORKING - testing the next line
                Ztmp*=0;
              }
            }
            if (outweight) wout+=Ztmp;
            if (Zprofile) reduceBitonicSort(itmp,Ztmp,iloc,Zloc,Zloc2,&Zprofile[k*NBINS]);
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
                itmp=NBINS2*((int)floor(q1*NBINS2))+(int)floor(q2*NBINS2);
                // if (weightprofile) assert(w<=plmd.Zprofile_d[k*NBINS+itmp]);
                if (weightprofile) Ztmp*=weightprofile[k*NBINS+itmp];
              }
              if (outweight) wout+=Ztmp;
              if (Zprofile) reduceBitonicSort(itmp,Ztmp,iloc,Zloc,Zloc2,&Zprofile[k*NBINS]);
              k++;
            }
          }
        }
      } else if (plmd.msprof) {
        for (i1=plmd.block0_d[s1]; i1<plmd.block0_d[s1+1]; i1++) {
          for (i2=plmd.block0_d[s2]; i2<plmd.block0_d[s2+1]; i2++) {
            __syncthreads();
            itmp=NBINS;
            Ztmp=w;
            if (b<plmd.B) {
              q1=lambda[i1];
              q2=lambda[i2];
              itmp=NBINS2*((int)floor(q1*NBINS2))+(int)floor(q2*NBINS2);
              // if (weightprofile) assert(w<=plmd.Zprofile_d[k*NBINS+itmp]);
              if (weightprofile) Ztmp*=weightprofile[k*NBINS+itmp];
            }
            if (outweight) wout+=Ztmp;
            if (Zprofile) reduceBitonicSort(itmp,Ztmp,iloc,Zloc,Zloc2,&Zprofile[k*NBINS]);
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

__global__
void freeenergykernel(struct_plmd plmd,real* Zprofile,real* G,real* L,real* dLdZprof)
{
  int i=blockIdx.x*blockDim.x+threadIdx.x;
  // int iprof=blockIdx.x;
  real Ztmp, Gtmp, dGtmp, kprof, Gavg, Ltmp;
  __shared__ real Lloc[NBINS];

  if (i<plmd.nprof*NBINS) {
    Ztmp=Zprofile[i];
    Gtmp=-plmd.kT*log(Ztmp);
  } else {
    Ztmp=0;
    Gtmp=0;
  }

  if (G) {
    if (i<plmd.nprof*NBINS) {
      G[i]=Gtmp;
    }
  }
  if (L || dLdZprof) {
    if (Ztmp) {
      kprof=plmd.kprofile_d[i];
      dGtmp=Gtmp-plmd.Gimp_d[i];
    } else {
      kprof=0;
      dGtmp=0;
    }
    reduceBroadcast(kprof*dGtmp,Lloc);
    Gavg=Lloc[threadIdx.x];
    reduceBroadcast(kprof,Lloc);
    if (Lloc[threadIdx.x]) Gavg/=Lloc[threadIdx.x];
    if (Ztmp) dGtmp-=Gavg; // Contribution of Gavg to derivatives of L with respect to Gtmp magically cancels out
  }
  if (L) {
    Ltmp=0.5*kprof*dGtmp*dGtmp;
    reduceNBINS(Ltmp,Lloc,L);
  }
  if (dLdZprof) {
    if (Ztmp) {
      dLdZprof[i]=kprof*dGtmp/Ztmp;
    }
  }
}

__global__
void regularizeLkernel(struct_plmd plmd)
{
  int i=blockIdx.x*blockDim.x+threadIdx.x;
  real deltax, L;
  __shared__ real Lloc[BLOCK];

  if (i<plmd.nx) {
    deltax=plmd.x_d[i]-plmd.xr_d[i];
    L=0.5*plmd.kx_d[i]*deltax*deltax;
  } else {
    L=0;
  }

  reduce(L,Lloc,plmd.L_d);
}

__global__
void regularizelinekernel(struct_plmd plmd,real s)
{
  int i=blockIdx.x*blockDim.x+threadIdx.x;
  real deltax, dxds, L, dLds;
  __shared__ real Lloc[BLOCK];

  if (i<plmd.nx) {
    dxds=plmd.dxds_d[i];
    deltax=plmd.x_d[i]+s*dxds-plmd.xr_d[i];
    L=0.5*plmd.kx_d[i]*deltax*deltax;
    dLds=plmd.kx_d[i]*deltax*dxds;
  } else {
    L=0;
    dLds=0;
  }

  reduce(L,Lloc,plmd.L_d);
  reduce(dLds,Lloc,plmd.dLds_d);
}

__global__
void regularizedLdxkernel(struct_plmd plmd)
{
  int i=blockIdx.x*blockDim.x+threadIdx.x;

  if (i<plmd.nx) {
    plmd.dLdx_d[i]=plmd.kx_d[i]*(plmd.x_d[i]-plmd.xr_d[i]);
  }
}

void evaluateGimp(struct_plmd *plmd)
{
  if (PROFILE) {
    cudaMemset(plmd->mc_Zprofile_d,0,plmd->nprof*NBINS*sizeof(real));
    // void profilekernel(struct_plmd plmd,real* lambda,real* inweight,real* weightprofile,real* outweight,real* Zprofile)
    profilekernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0],plmd->mc_lambda_d,plmd->mc_ensweight_d,NULL,NULL,plmd->mc_Zprofile_d);
    // freeenergykernel<<<(plmd->nprof*NBINS+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0],plmd->mc_Zprofile_d,plmd->Gimp_d,0,NULL,NULL,NULL);
    // Too noisy:
    // freeenergykernel<<<plmd->nprof,NBINS>>>(plmd[0],plmd->mc_Zprofile_d,plmd->Gimp_d,NULL,NULL);
    int s1,s2,i,j,k,kn;
    real *mc_Zprofile, *Gimp;
    mc_Zprofile=(real*)calloc(NBINS*plmd->nprof,sizeof(real));
    Gimp=(real*)calloc(NBINS,sizeof(real));
    cudaMemcpy(mc_Zprofile,plmd->mc_Zprofile_d,NBINS*plmd->nprof*sizeof(real),cudaMemcpyDeviceToHost);
    k=0;
    for (s1=0; s1<plmd->nsites; s1++) {
      for (s2=s1; s2<plmd->nsites; s2++) {
        if (s1==s2) { // Same site
          kn=k+plmd->nsubs[s1];
          for (i=0; i<NBINS; i++) {
            Gimp[i]=0;
            for (j=k; j<kn; j++) {
              Gimp[i]+=mc_Zprofile[NBINS*j+i];
            }
            Gimp[i]=-plmd->kT*log(Gimp[i]);
          }
          for (j=k; j<kn; j++) {
            cudaMemcpy(&plmd->Gimp_d[NBINS*j],Gimp,NBINS*sizeof(real),cudaMemcpyHostToDevice);
          }
          k=kn;

          kn=k+(plmd->nsubs[s1]*(plmd->nsubs[s1]-1))/2;
          for (i=0; i<NBINS; i++) {
            Gimp[i]=0;
            for (j=k; j<kn; j++) {
              Gimp[i]+=mc_Zprofile[NBINS*j+i];
            }
            Gimp[i]=-plmd->kT*log(Gimp[i]);
          }
          for (j=k; j<kn; j++) {
            cudaMemcpy(&plmd->Gimp_d[NBINS*j],Gimp,NBINS*sizeof(real),cudaMemcpyHostToDevice);
          }
          k=kn;

          if (plmd->nsubs[s1]>2) {
            kn=k+(plmd->nsubs[s1]*(plmd->nsubs[s1]-1))/2;
            for (i=0; i<NBINS; i++) {
              Gimp[i]=0;
              for (j=k; j<kn; j++) {
                Gimp[i]+=mc_Zprofile[NBINS*j+i];
              }
              Gimp[i]=-plmd->kT*log(Gimp[i]);
            }
            for (j=k; j<kn; j++) {
              cudaMemcpy(&plmd->Gimp_d[NBINS*j],Gimp,NBINS*sizeof(real),cudaMemcpyHostToDevice);
            }
            k=kn;
          }
        } else if (plmd->msprof) {
          kn=k+plmd->nsubs[s1]*plmd->nsubs[s2];
          for (i=0; i<NBINS; i++) {
            Gimp[i]=0;
            for (j=k; j<kn; j++) {
              Gimp[i]+=mc_Zprofile[NBINS*j+i];
            }
            Gimp[i]=-plmd->kT*log(Gimp[i]);
          }
          for (j=k; j<kn; j++) {
            cudaMemcpy(&plmd->Gimp_d[NBINS*j],Gimp,NBINS*sizeof(real),cudaMemcpyHostToDevice);
          }
          k=kn;
        }
      }
    }
    free(mc_Zprofile);
    free(Gimp);
  }

  if (MOMENT) {
    int i;
    real sum;
    sum=0;
    for (i=0; i<plmd->B; i++) {
      sum+=plmd->ensweight[i];
    }
    cudaMemcpy(plmd->sumensweight_d,&sum,sizeof(real),cudaMemcpyHostToDevice);
    cudaMemset(plmd->moments_d,0,plmd->nbias*sizeof(real));
// void weightedenergykernel(struct_plmd plmd,real sign,real* lambda,real* weight,real* dEdx)
    weightedenergykernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0],-1,plmd->lambda_d,plmd->ensweight_d,plmd->moments_d);
  }
}

__global__
void likelihoodkernel(struct_plmd plmd,real s,real* L,real* dLds)
{
  if (L) {
    atomicAdd(L,plmd.Esum_d[0]/(plmd.sumensweight_d[0]*plmd.kT));
    if (s) atomicAdd(L,s*plmd.dEdssum_d[0]/(plmd.sumensweight_d[0]*plmd.kT));
    atomicAdd(L,log(plmd.mc_Z_d[0]));
  }
  if (dLds) {
    atomicAdd(dLds,plmd.dEdssum_d[0]/(plmd.sumensweight_d[0]*plmd.kT));
    atomicAdd(dLds,-plmd.mc_dEdssum_d[0]/(plmd.mc_Z_d[0]*plmd.kT));
  }
}

__global__
void gradientlikelihoodkernel(struct_plmd plmd,real* norm,real* dLdxin)
{
  int i=blockIdx.x*blockDim.x+threadIdx.x;

  if (i<plmd.nbias) {
    atomicAdd(&plmd.dLdx_d[i],-dLdxin[i]/(norm[0]*plmd.kT));
  }
}

void evaluateL(struct_plmd *plmd)
{
  /*{ // DEBUG
  FILE *fp;
  int i;
  fp=fopen("xbad.dat","r");
  for (i=0;i<plmd->nx;i++) {
    double buffer;
    fscanf(fp,"%lf",&buffer);
    plmd->x[i]=buffer;
  }
  fclose(fp);
  }*/

  cudaMemcpy(plmd->x_d,plmd->x,plmd->nx*sizeof(real),cudaMemcpyHostToDevice);
  cudaMemset(plmd->L_d,0,sizeof(real));

  regularizeLkernel<<<(plmd->nx+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0]);

  energykernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0],plmd->x_d,plmd->lambda_d,plmd->E_d);

  // cudaMemcpy(plmd->L,plmd->E_d,sizeof(real),cudaMemcpyDeviceToHost); // DEBUG
  // fprintf(stderr,"Debug    energy[0]=%lg\n",plmd->L[0]); // DEBUG

  if (PROFILE) {
    cudaMemset(plmd->Z_d,0,sizeof(real));
// void boltzmannkernel(struct_plmd plmd,real sign,real* energy,real s,real* denergyds,real* inweight,real* outweight,real* Z)
    boltzmannkernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>
      (plmd[0],-1,plmd->E_d,0,NULL,plmd->ensweight_d,plmd->weight_d,plmd->Z_d);

    cudaMemset(plmd->Zprofile_d,0,plmd->nprof*NBINS*sizeof(real));
// void profilekernel(struct_plmd plmd,real* lambda,real* inweight,real* weightprofile,real* outweight,real* Zprofile)
    profilekernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0],plmd->lambda_d,plmd->weight_d,NULL,NULL,plmd->Zprofile_d);
    freeenergykernel<<<plmd->nprof,NBINS>>>(plmd[0],plmd->Zprofile_d,plmd->G_d,plmd->L_d,NULL);
  }

  if (MOMENT) {
    energykernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0],plmd->x_d,plmd->mc_lambda_d,plmd->mc_E_d);
    cudaMemset(plmd->Esum_d,0,sizeof(real));
    dotenergykernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0],1,plmd->ensweight_d,plmd->E_d,plmd->Esum_d);
    cudaMemset(plmd->mc_Z_d,0,sizeof(real));
    boltzmannkernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>
      (plmd[0],1,plmd->mc_E_d,0,NULL,plmd->mc_ensweight_d,plmd->mc_weight_d,plmd->mc_Z_d);
    likelihoodkernel<<<1,1>>>(plmd[0],0,plmd->L_d,NULL);
  }

  cudaMemcpy(plmd->L,plmd->L_d,sizeof(real),cudaMemcpyDeviceToHost);

  fprintf(stdout,"New      L=%lg\n",plmd->L[0]);
}

void evaluateL_line(real s,struct_plmd *plmd)
{
  cudaMemset(plmd->L_d,0,sizeof(real));
  cudaMemset(plmd->dLds_d,0,sizeof(real));

  regularizelinekernel<<<(plmd->nx+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0],s);

  if (PROFILE) {
    cudaMemset(plmd->Z_d,0,sizeof(real));
// void boltzmannkernel(struct_plmd plmd,real sign,real* energy,real s,real* denergyds,real* inweight,real* outweight,real* Z)
    boltzmannkernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>
      (plmd[0],-1,plmd->E_d,s,plmd->dEds_d,plmd->ensweight_d,plmd->weight_d,plmd->Z_d);

    cudaMemset(plmd->Zprofile_d,0,plmd->nprof*NBINS*sizeof(real));
// void profilekernel(struct_plmd plmd,real* lambda,real* inweight,real* weightprofile,real* outweight,real* Zprofile)
    profilekernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0],plmd->lambda_d,plmd->weight_d,NULL,NULL,plmd->Zprofile_d);
    freeenergykernel<<<plmd->nprof,NBINS>>>(plmd[0],plmd->Zprofile_d,plmd->G_d,plmd->L_d,plmd->dLdZprofile_d);

    profilekernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0],plmd->lambda_d,plmd->weight_d,plmd->dLdZprofile_d,plmd->dLdE_d,NULL);
    dotenergykernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0],-1,plmd->dLdE_d,plmd->dEds_d,plmd->dLds_d);
  }

  if (MOMENT) {
    cudaMemset(plmd->mc_Z_d,0,sizeof(real));
    boltzmannkernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>
      (plmd[0],1,plmd->mc_E_d,s,plmd->mc_dEds_d,plmd->mc_ensweight_d,plmd->mc_weight_d,plmd->mc_Z_d);
    cudaMemset(plmd->mc_dEdssum_d,0,sizeof(real));
    dotenergykernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0],1,plmd->mc_dEds_d,plmd->mc_weight_d,plmd->mc_dEdssum_d);
    likelihoodkernel<<<1,1>>>(plmd[0],s,plmd->L_d,plmd->dLds_d);
  }

  /*if (s==64) { // DEBUG
  printarray(plmd->dLdE_d,plmd->B,"dLdE.dat");
  printarray(plmd->dLds_d,1,"dLds.dat");
  printarray(plmd->Zprofile_d,NBINS*plmd->nprof,"Zprofile.dat");
  printarray(plmd->dLdZprofile_d,NBINS*plmd->nprof,"dLdZprofile.dat");
  }*/

  cudaMemcpy(plmd->L,plmd->L_d,sizeof(real),cudaMemcpyDeviceToHost);
  cudaMemcpy(plmd->dLds,plmd->dLds_d,sizeof(real),cudaMemcpyDeviceToHost);
}

void evaluatedLdx(struct_plmd *plmd)
{
  cudaMemset(plmd->dLdx_d,0,plmd->nx*sizeof(real));

  regularizedLdxkernel<<<(plmd->nx+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0]);

  if (PROFILE) {
    freeenergykernel<<<plmd->nprof,NBINS>>>(plmd[0],plmd->Zprofile_d,plmd->G_d,NULL,plmd->dLdZprofile_d);
    profilekernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0],plmd->lambda_d,plmd->weight_d,plmd->dLdZprofile_d,plmd->dLdE_d,NULL);

// void weightedenergykernel(struct_plmd plmd,real* lambda,real* weight,real* dEdx)
    weightedenergykernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>
      (plmd[0],-1,plmd->lambda_d,plmd->dLdE_d,plmd->dLdx_d);
  }

  if (MOMENT) {
    cudaMemset(plmd->mc_moments_d,0,plmd->nbias*sizeof(real));
// void weightedenergykernel(struct_plmd plmd,real sign,real* lambda,real* weight,real* dEdx)
    weightedenergykernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0],1,plmd->mc_lambda_d,plmd->mc_weight_d,plmd->mc_moments_d);

// void gradientlikelihoodkernel(struct_plmd plmd,real* norm,real* dLdxin)
    gradientlikelihoodkernel<<<(plmd->nbias+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0],plmd->sumensweight_d,plmd->moments_d);
    gradientlikelihoodkernel<<<(plmd->nbias+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0],plmd->mc_Z_d,plmd->mc_moments_d);
  }
  
  cudaMemcpy(plmd->dLdx,plmd->dLdx_d,plmd->nx*sizeof(real),cudaMemcpyDeviceToHost);
}

void resetHinv(struct_plmd *plmd)
{
  int i;
  // N^2 Hinv
  // for (i=0; i<plmd->nx*plmd->nx; i++) {
  //   plmd->Hinv[i]=0.0L;
  // }
  // for (i=0; i<plmd->nx; i++) {
  //   plmd->Hinv[(plmd->nx+1)*i]=1.0L;
  // }
  // N^1 Hinv
  for (i=0; i<plmd->nx; i++) {
    plmd->x0[i]=plmd->x[i];
    plmd->dLdx0[i]=plmd->dLdx[i];
  }
}

void updateHinv(struct_plmd *plmd)
{
  /* // Begin N^2 Hinv
  int i,j;
  real DxDg,DgHinvDg;
  real c1,c2;

  DxDg=0.0L;
  for (i=0;i<plmd->nx;i++) {
    // Put Delta x and Delta dLdx in x0 and dLdx0, which hold previous values
    plmd->x0[i]=plmd->x[i]-plmd->x0[i];
    plmd->dLdx0[i]=plmd->dLdx[i]-plmd->dLdx0[i];
    DxDg+=plmd->x0[i]*plmd->dLdx0[i];
  }

  DgHinvDg=0.0L;
  for (i=0; i<plmd->nx; i++) {
    plmd->hi[i]=0.0L;
    for (j=0; j<plmd->nx; j++) {
      // put Hinv * Delta dLdx in hi (the search direction) as a buffer
      plmd->hi[i]+=plmd->Hinv[i*plmd->nx+j]*plmd->dLdx0[j];
    }
    DgHinvDg+=plmd->hi[i]*plmd->dLdx0[i];
  }
  c1=(1.0L+DgHinvDg/DxDg)/DxDg;
  c2=-1.0L/DxDg;

  for (i=0; i<plmd->nx; i++) {
    for (j=i; j<plmd->nx; j++) {
      plmd->Hinv[i*plmd->nx+j]+=c1*(plmd->x0[i]*plmd->x0[j])+c2*(plmd->hi[i]*plmd->x0[j]+plmd->x0[i]*plmd->hi[j]);
      plmd->Hinv[j*plmd->nx+i]=plmd->Hinv[i*plmd->nx+j];
    }
  }
  
  //   dx=xf-xi;
  //   dd=df-di;
  //   Hinv=Hinv+(1+(dd'*Hinv*dd)/(dx'*dd))*(dx*dx')/(dx'*dd)-((Hinv*dd*dx')+(Hinv*dd*dx')')/(dx'*dd);
  */ // End N^2 Hinv
  // Begin N^1 Hinv
  int i,j;

  if (plmd->Nmem<plmd->Nmemax) {
    plmd->Nmem++;
  }
  for (i=plmd->Nmem-1; i>0; i--) {
    for (j=0; j<plmd->nx; j++) {
      plmd->d_x[i*plmd->nx+j]=plmd->d_x[(i-1)*plmd->nx+j];
      plmd->d_dLdx[i*plmd->nx+j]=plmd->d_dLdx[(i-1)*plmd->nx+j];
    }
    plmd->rho[i]=plmd->rho[i-1];
  }

  plmd->rho[0]=0;
  for (i=0; i<plmd->nx; i++) {
    plmd->d_x[i]=plmd->x[i]-plmd->x0[i];
    plmd->d_dLdx[i]=plmd->dLdx[i]-plmd->dLdx0[i];
    plmd->rho[0]+=plmd->d_x[i]*plmd->d_dLdx[i];
  }
  plmd->rho[0]=1.0/plmd->rho[0];

  for (i=0; i<plmd->nx; i++) {
    plmd->x0[i]=plmd->x[i];
    plmd->dLdx0[i]=plmd->dLdx[i];
  }
  // End N^1 Hinv
}

void projectHinv(struct_plmd *plmd)
{
/* // Begin N^2 Hinv
  int i,j;
  real dLds;
  //   hi=-Hinv*df;
  for (i=0; i<plmd->nx; i++) {
    plmd->hi[i]=0;
    for (j=0; j<plmd->nx; j++) {
      plmd->hi[i]+=-plmd->Hinv[i*plmd->nx+j]*plmd->dLdx[j];
    }
  }

  dLds=0;
  for (i=0; i<plmd->nx; i++) {
    dLds+=plmd->hi[i]*plmd->dLdx[i];
  }

  if (dLds>0) {
    fprintf(stderr,"Bad direction, reset Hinv\n");
    for (i=0; i<plmd->nx*plmd->nx; i++) {
      plmd->Hinv[i]=0;
    }
    for (i=0; i<plmd->nx; i++) {
      plmd->Hinv[i*(plmd->nx+1)]=1;
      plmd->hi[i]=plmd->dLdx[i];
    }
  }
*/ // End N^2 Hinv
// Begin N^1 Hinv
  int i,j;

  for (i=0; i<plmd->nx; i++) {
    plmd->hi[i]=plmd->dLdx[i];
  }
  for (i=0; i<plmd->Nmem; i++) {
    plmd->alpha[i]=0;
    for (j=0; j<plmd->nx; j++) {
      plmd->alpha[i]+=plmd->d_x[i*plmd->nx+j]*plmd->hi[j];
    }
    plmd->alpha[i]*=plmd->rho[i];
    for (j=0; j<plmd->nx; j++) {
      plmd->hi[j]+=-plmd->alpha[i]*plmd->d_dLdx[i*plmd->nx+j];
    }
  }
  /*
  // According to wikipedia, this is to ensure the step length is always about unity
  // https://en.wikipedia.org/wiki/Limited-memory_BFGS
  if (plmd->Nmem>0) {
    numer=0.0L
    denom=0.0L;
    for (i=0; i<plmd->nx; i++) {
      numer+=plmd->d_x[i]*plmd->hi[i];
      denom+=plmd->d_dLdx[i]*plmd->d_dLdx[i];
    }
    numer/=denom;
    for (i=0; i<plmd->nx; i++) {
      plmd->hi[i]=numer*plmd->d_dLdx[i]; // This seems like a horrible idea, maybe wikipedia has a typo...
    }
  }
  */
  for (i=plmd->Nmem-1; i>=0; i--) {
    plmd->beta[i]=0;
    for (j=0; j<plmd->nx; j++) {
      plmd->beta[i]+=plmd->d_dLdx[i*plmd->nx+j]*plmd->hi[j];
    }
    plmd->beta[i]*=plmd->rho[i];
    for (j=0; j<plmd->nx; j++) {
      plmd->hi[j]+=(plmd->alpha[i]-plmd->beta[i])*plmd->d_x[i*plmd->nx+j];
    }
  }

  for (i=0; i<plmd->nx; i++) {
    plmd->hi[i]*=-1;
  }
// End N^1 Hinv

  /*{ // DEBUG
  FILE *fp;
  int i;
  fp=fopen("dxdsbad.dat","r");
  for (i=0;i<plmd->nx;i++) {
    double buffer;
    fscanf(fp,"%lf",&buffer);
    plmd->hi[i]=buffer;
  }
  fclose(fp);
  }*/

  cudaMemcpy(plmd->dxds_d,plmd->hi,plmd->nx*sizeof(real),cudaMemcpyHostToDevice);

  energykernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0],plmd->dxds_d,plmd->lambda_d,plmd->dEds_d);
  if (MOMENT) {
    energykernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0],plmd->dxds_d,plmd->mc_lambda_d,plmd->mc_dEds_d);
    cudaMemset(plmd->dEdssum_d,0,sizeof(real));
    dotenergykernel<<<(plmd->B+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0],1,plmd->ensweight_d,plmd->dEds_d,plmd->dEdssum_d);
  }
}

void update_line(int step,struct_plmd *plmd)
{
  int i;
  real a,b,c,s;
  real s1,s2,s3;
  real L1,L2,L3;
  real dLds1,dLds2,dLds3;
  real L0;

  for (i=0; i<plmd->nx; i++) {
    plmd->x0[i]=plmd->x[i];
    plmd->dLdx0[i]=plmd->dLdx[i];
  }

  L0=plmd->L[0];

  s1=0.0;
  evaluateL_line(s1,plmd);
  L1=plmd->L[0];
  dLds1=plmd->dLds[0];
  if (dLds1>0) {
    fprintf(stdout,"Error, hi is pointing wrong way - halting\n");
    plmd->done=true;
    return;
    // exit(1);
  }
  
  s3=1.0;
  evaluateL_line(s3,plmd);
  L3=plmd->L[0];
  dLds3=plmd->dLds[0];

  while (dLds3<0 && s3<1e+8) {
    fprintf(stdout,"Seek %4d s=%lg %lg\n          L=%lg %lg\n       dLds=%lg %lg\n",
            step,(double) s1,(double) s3,
            (double) L1,(double) L3,
            (double) dLds1,(double) dLds3);
    s2=s1-dLds1*(s3-s1)/(dLds3-dLds1);
    s3=((1.5*s2>8*s3 || 1.5*s2<=0) ? 8*s3 : 1.5*s2); // s2 is expected 0. Go past it by 50%, unless that's an increase of more than a factor of 8.
    evaluateL_line(s3,plmd);
    L3=plmd->L[0];
    dLds3=plmd->dLds[0];
  }

  while (!isfinite(dLds3) && s3>1e-8) {
    fprintf(stdout,"Warning, overshot bound\n");
    fprintf(stdout,"Seek %4d s=%lg %lg\n          L=%lg %lg\n       dLds=%lg %lg\n",
            step,(double) s1,(double) s3,
            (double) L1,(double) L3,
            (double) dLds1,(double) dLds3);
    s3=0.95*s3;
    evaluateL_line(s3,plmd);
    L3=plmd->L[0];
    dLds3=plmd->dLds[0];
  }

  if (!(dLds3>0)) {
    fprintf(stdout,"Warning: Step %4d unsuccessful, halting minimization\n",step);
    fprintf(stdout,"Seek %4d s=%lg %lg\n          L=%lg %lg\n       dLds=%lg %lg\n",
            step,(double) s1,(double) s3,
            (double) L1,(double) L3,
            (double) dLds1,(double) dLds3);
    plmd->done=true;
    return;
  }

  /*if (s3>6) { // DEBUG
    for (i=0; i<10; i++) {
      evaluateL_line(0.1*i*s3,plmd);
      fprintf(stderr,"Debug   s=%g L=%g dLds=%g\n",0.1*i*s3,plmd->L[0],plmd->dLds[0]);
    }
  }*/

  /*if (!isfinite(dLds3)) { // DEBUG
    printarray(plmd->x_d,plmd->nx,"xbad.dat");
    printarray(plmd->dxds_d,plmd->nx,"dxdsbad.dat");
    exit(1);
  }*/

  s2=s1-dLds1*(s3-s1)/(dLds3-dLds1);
  evaluateL_line(s2,plmd);
  L2=plmd->L[0];
  dLds2=plmd->dLds[0];

  fprintf(stdout,"Step %4d s=%lg %lg %lg\n          L=%lg %lg %lg\n       dLds=%lg %lg %lg\n",
          step,(double) s1,(double) s2,(double) s3,
          (double) L1,(double) L2,(double) L3,
          (double) dLds1,(double) dLds2,(double) dLds3);

  for (i=0; i<15; i++) {
    if ((s2-s1)/s2<5e-7 || (s3-s2)/s2<5e-7 || dLds2==0) break;

    // Quadratic interpolation
    a=dLds1/((s1-s2)*(s1-s3));
    a+=dLds2/((s2-s1)*(s2-s3));
    a+=dLds3/((s3-s1)*(s3-s2));
    b=-dLds1*(s2+s3)/((s1-s2)*(s1-s3));
    b+=-dLds2*(s1+s3)/((s2-s1)*(s2-s3));
    b+=-dLds3*(s1+s2)/((s3-s1)*(s3-s2));
    c=dLds1*s2*s3/((s1-s2)*(s1-s3));
    c+=dLds2*s1*s3/((s2-s1)*(s2-s3));
    c+=dLds3*s1*s2/((s3-s1)*(s3-s2));
    s=(-b+sqrt(b*b-4*a*c))/(2*a);

    if (dLds2<0) {
      s1=s2;
      L1=L2;
      dLds1=dLds2;
    } else { // dLds2==0 already addressed above
      s3=s2;
      L3=L2;
      dLds3=dLds2;
    }

    if (s>s1 && s<s3) {
      // Use the earlier quadratic interpolation
      s2=s;
    } else {
      // Linear interpolation (secant method)
      fprintf(stdout,"Warning, fell back on linear interpolation\n");
      fprintf(stdout,"a=%lg b=%lg c=%lg s-=%lg s+=%lg s=%lg\n",(double)a,(double)b,(double)c,(double)((-b-sqrt(b*b-4*a*c))/(2*a)),(double)((-b+sqrt(b*b-4*a*c))/(2*a)),(double)s);
      s2=s1-dLds1*(s3-s1)/(dLds3-dLds1);
    }

    evaluateL_line(s2,plmd);
    L2=plmd->L[0];
    dLds2=plmd->dLds[0];

    fprintf(stdout,"Step %4d s=%lg %lg %lg\n          L=%lg %lg %lg\n       dLds=%lg %lg %lg\n",
            step,(double) s1,(double) s2,(double) s3,
            (double) L1,(double) L2,(double) L3,
            (double) dLds1,(double) dLds2,(double) dLds3);
  }

  fprintf(stdout,"Step %4d s=%lg %lg %lg\n          L=%lg %lg %lg\n       dLds=%lg %lg %lg\n",
          step,(double) s1,(double) s2,(double) s3,
          (double) L1,(double) L2,(double) L3,
          (double) dLds1,(double) dLds2,(double) dLds3);

  real stepLength2=0;
  real initGrad2=0;

  for (i=0; i<plmd->nx; i++) {
    plmd->x[i]=plmd->x0[i]+s2*plmd->hi[i];
    stepLength2+=(s2*plmd->hi[i])*(s2*plmd->hi[i]);
    initGrad2+=plmd->dLdx[i]*plmd->dLdx[i];
  }

  // fprintf(stderr,"Step %d smid1 %lg smid2 %lg L=%lg -> L1=%lg L2=%lg\n",step,(double) smid1,(double) smid2,(double) L0,(double) Lmid1,(double) Lmid2);
  fprintf(stdout,"Step %4d L=%24.16lf -> L2=%24.16lf, dL=%lg, step length=%lg\n",step,(double)L0,(double)L2,(double)(L2-L0),(double) sqrt(stepLength2));

  fprintf(plmd->fplog,"%24.16lf %24.16lf %lg %lg\n",(double)L0,(double)L2,(double)sqrt(initGrad2),(double)sqrt(stepLength2));

  if (sqrt(stepLength2)<5e-7) plmd->done=true;
  if (sqrt(stepLength2/plmd->nx)<plmd->criteria) { // criteria was 1e-2
    plmd->doneCount+=1;
    if (plmd->doneCount==2) plmd->done=true;
  } else {
    plmd->doneCount=0;
  }
}

real lineL(real s,struct_plmd *plmd)
{
  int i;

  for (i=0; i<plmd->nx; i++) {
    plmd->x[i]=plmd->x0[i]+s*plmd->hi[i];
  }

  evaluateL(plmd);
  return plmd->L[0];
}

void itterate(int step,struct_plmd *plmd)
{
  evaluateL(plmd);
  evaluatedLdx(plmd);

  if (step==0) {
    resetHinv(plmd);
  } else {
    updateHinv(plmd);
  }

  projectHinv(plmd);
  update_line(step,plmd);
}

void run(struct_plmd *plmd)
{
  int s;

  evaluateGimp(plmd);
  plmd->done=false;
  plmd->doneCount=0;
  for (s=0; s<250; s++) {
    itterate(s,plmd);
    if (plmd->done) break;
  }
}

void finish(struct_plmd *plmd,int argc, char *argv[])
{
  int i,j;
  FILE *fp;

  fp=fopen(argv[6],"w");
  for (i=0; i<plmd->nx; i++) {
    fprintf(fp," %lg",(double) plmd->x[i]);
  }
  fclose(fp);

  plmd->mc_weight=(real*)calloc(plmd->B,sizeof(real));
  cudaMemcpy(plmd->mc_weight,plmd->mc_weight_d,plmd->B*sizeof(real),cudaMemcpyDeviceToHost);
  cudaMemcpy(plmd->mc_lambda,plmd->mc_lambda_d,plmd->B*plmd->nblocks*sizeof(real),cudaMemcpyDeviceToHost);

  fp=fopen("mc_weight.dat","w");
  for (i=0; i<plmd->B; i++) {
    fprintf(fp," %lg\n",(double) plmd->mc_weight[i]);
  }
  fclose(fp);

  fp=fopen("mc_Lambda.dat","w");
  for (i=0; i<plmd->B; i++) {
    for (j=0; j<plmd->nblocks; j++) {
      fprintf(fp," %lg",(double) plmd->mc_lambda[i*plmd->nblocks+j]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  free(plmd->mc_weight);

  free(plmd->nsubs);
  free(plmd->block0);
  cudaFree(plmd->block0_d);
  free(plmd->block2site);

  free(plmd->lambda);
  cudaFree(plmd->lambda_d);
  free(plmd->ensweight);
  cudaFree(plmd->ensweight_d);
  free(plmd->mc_lambda);
  cudaFree(plmd->mc_lambda_d);
  free(plmd->mc_ensweight);
  cudaFree(plmd->mc_ensweight_d);
  free(plmd->kx);
  cudaFree(plmd->kx_d);
  free(plmd->xr);
  cudaFree(plmd->xr_d);

  free(plmd->L);
  cudaFree(plmd->L_d);
  free(plmd->dLds);
  cudaFree(plmd->dLds_d);

  free(plmd->x);
  free(plmd->dLdx);
  cudaFree(plmd->dLdx_d);
  // free(plmd->Hinv);
  free(plmd->d_x);
  free(plmd->d_dLdx);
  free(plmd->rho);
  free(plmd->alpha);
  free(plmd->beta);
  free(plmd->hi);
  free(plmd->x0);
  free(plmd->dLdx0);

  cudaFree(plmd->E_d);
  cudaFree(plmd->dEds_d);
  cudaFree(plmd->mc_E_d);
  cudaFree(plmd->mc_dEds_d);
  cudaFree(plmd->weight_d);
  cudaFree(plmd->mc_weight_d);
  cudaFree(plmd->x_d);
  cudaFree(plmd->dxds_d);
  cudaFree(plmd->Z_d);
  cudaFree(plmd->mc_Z_d);
  cudaFree(plmd->Zprofile_d);
  cudaFree(plmd->mc_Zprofile_d);
  cudaFree(plmd->dLdZprofile_d);
  cudaFree(plmd->mc_dLdZprofile_d);
  cudaFree(plmd->dLdE_d);
  cudaFree(plmd->Gimp_d);
  cudaFree(plmd->G_d);

  cudaFree(plmd->Esum_d);
  cudaFree(plmd->dEdssum_d);
  cudaFree(plmd->mc_dEdssum_d);
  cudaFree(plmd->moments_d);
  cudaFree(plmd->mc_moments_d);
  cudaFree(plmd->sumensweight_d);

  fclose(plmd->fplog);

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
