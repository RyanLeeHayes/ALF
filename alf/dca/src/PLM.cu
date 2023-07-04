// Written by Ryan Hayes 2017-06-20
// plmDCA algorithm from R470 - DOI: 10.1016/j.jcp.2014.07.024
// Quasi newton equations from https://www.rose-hulman.edu/~bryan/lottamath/quasinewton.pdf

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <unistd.h>

//   ID=omp_get_thread_num();
//   NID=omp_get_max_threads();

#include "PLM.h"

#define MAXLENGTH 1024
#define PLMGPU

struct_plmd* setup(int argc, char *argv[])
{
  struct_plmd *plmd;
  int h,i,j,k,l;
  FILE *fp;
  char line[MAXLENGTH];
  char *linebuf;

  if (argc<5) {
    fprintf(stderr,"Error: not enough input arguments\n");
    exit(1);
  }

  plmd=(struct_plmd*) malloc(sizeof(struct_plmd));

  fp=fopen("../prep/nsubs","r");
  if (fp==NULL) {
    fprintf(stderr,"Error, ../prep/nsubs does not exist\n");
    exit(1);
  }
  plmd->nsites=0;
  while (fscanf(fp,"%d",&i)==1) {
    plmd->nsites++;
  }
  fclose(fp);

  fp=fopen("../prep/nsubs","r");
  i=0;
  plmd->nsubs=(int*) calloc(plmd->nsites,sizeof(int));
  plmd->nblocks=0;
  plmd->nsubsmax=0;
  for(i=0; i<plmd->nsites; i++) {
    fscanf(fp,"%d",&(plmd->nsubs[i]));
    plmd->nsubs[i]++;
    plmd->nblocks+=plmd->nsubs[i];
    plmd->nsubsmax=(plmd->nsubs[i]>plmd->nsubsmax?plmd->nsubs[i]:plmd->nsubsmax);
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

#ifdef PLMGPU
  cudaMalloc(&plmd->block0_d,(plmd->nsites+1)*sizeof(int));
  cudaMemcpy(plmd->block0_d,plmd->block0,(plmd->nsites+1)*sizeof(int),cudaMemcpyHostToDevice);
  cudaMalloc(&plmd->nsubs_d,plmd->nsites*sizeof(int));
  cudaMemcpy(plmd->nsubs_d,plmd->nsubs,plmd->nsites*sizeof(int),cudaMemcpyHostToDevice);
#endif

  plmd->h=(real*) calloc(plmd->nblocks,sizeof(real));
  plmd->J=(real*) calloc(plmd->nblocks*plmd->nblocks,sizeof(real));

  plmd->dhds=(real*) calloc(plmd->nblocks,sizeof(real));
  plmd->dJds=(real*) calloc(plmd->nblocks*plmd->nblocks,sizeof(real));

  plmd->dLdh=(reala*) calloc(plmd->nblocks,sizeof(reala));
  plmd->dLdJ=(reala*) calloc(plmd->nblocks*plmd->nblocks,sizeof(reala));

#ifdef PLMGPU
  cudaMalloc(&plmd->h_d,plmd->nblocks*sizeof(real));
  cudaMalloc(&plmd->J_d,plmd->nblocks*plmd->nblocks*sizeof(real));

  cudaMalloc(&plmd->dhds_d,plmd->nblocks*sizeof(real));
  cudaMalloc(&plmd->dJds_d,plmd->nblocks*plmd->nblocks*sizeof(real));

  cudaMalloc(&plmd->dLdh_d,plmd->nblocks*sizeof(reala));
  cudaMalloc(&plmd->dLdJ_d,plmd->nblocks*plmd->nblocks*sizeof(reala));
#endif

  plmd->B=0;
  for (h=0; h<argc-3; h++) {
  fp=fopen(argv[h+3],"r");
  for (; fgets(line,MAXLENGTH,fp) != NULL; plmd->B++) {
    ;
  }
  fclose(fp);
  }
  fprintf(stderr,"%d frames\n",plmd->B); // DEBUG

  plmd->BiBmax=(1<<20);
  if (plmd->B>plmd->BiBmax) { // Too much data, have to save to temporary files
    plmd->Seq_fp=tmpfile();
    plmd->BiB=plmd->BiBmax;
  } else { // All fits in memory
    plmd->BiB=plmd->B;
  }

  plmd->Seq=(int*) calloc(plmd->BiB*plmd->nsites,sizeof(int));
  i=0;
  for (h=0; h<argc-3; h++) {
  fp=fopen(argv[h+3],"r");
  for (; fgets(line,MAXLENGTH,fp) != NULL; i++) {
    linebuf=line;
    k=0;
    for (j=0;j<plmd->nsites;j++) {
      sscanf(linebuf,"%d%n",&(plmd->Seq[i*plmd->nsites+j]),&l);
      linebuf+=l;
      plmd->Seq[i*plmd->nsites+j]+=k;
      k+=plmd->nsubs[j];
    }
    if (plmd->B>plmd->BiBmax && i>=plmd->BiBmax) {
      fwrite(plmd->Seq,sizeof(int),i*plmd->nsites,plmd->Seq_fp);
      i=0;
    }
  }
  fclose(fp);
  }
  if (plmd->B>plmd->BiBmax) {
    fwrite(plmd->Seq,sizeof(int),i*plmd->nsites,plmd->Seq_fp);
  }

#ifdef PLMGPU
  cudaMalloc(&plmd->Seq_d,plmd->BiB*plmd->nsites*sizeof(int));
  cudaMemcpy(plmd->Seq_d,plmd->Seq,plmd->BiB*plmd->nsites*sizeof(int),cudaMemcpyHostToDevice);
#endif

  plmd->heff=(real*) calloc(plmd->BiB*plmd->nblocks,sizeof(real));
  plmd->dheffds=(real*) calloc(plmd->BiB*plmd->nblocks,sizeof(real));

#ifdef PLMGPU
  cudaMalloc(&plmd->heff_d,plmd->BiB*plmd->nblocks*sizeof(real));
  cudaMalloc(&plmd->dheffds_d,plmd->BiB*plmd->nblocks*sizeof(real));
#endif

  plmd->kh=1e-6;
  plmd->kJ=1e-6;

  plmd->L=(reala*) calloc(1,sizeof(reala));
  plmd->dLds=(reala*) calloc(1,sizeof(reala));

#ifdef PLMGPU
  cudaMalloc(&plmd->L_d,sizeof(reala));
  cudaMalloc(&plmd->dLds_d,sizeof(reala));
#endif


  k=plmd->nblocks;
  for (i=0; i<plmd->nblocks; i++) {
    for (j=plmd->block0[plmd->block2site[i]+1]; j<plmd->nblocks; j++) {
      k++;
    }
  }
  plmd->Jend=k;

  plmd->x=(real*) calloc(plmd->Jend,sizeof(real));
  plmd->dLdx=(reala*) calloc(plmd->Jend,sizeof(reala));
  // plmd->Hinv=(real*) calloc(plmd->Jend*plmd->Jend,sizeof(real));
  plmd->Nmemax=200;
  plmd->Nmem=0;
  plmd->d_x=(real*) calloc(plmd->Jend*plmd->Nmemax,sizeof(real));
  plmd->d_dLdx=(real*) calloc(plmd->Jend*plmd->Nmemax,sizeof(real));
  plmd->rho=(real*) calloc(plmd->Nmemax,sizeof(real));
  plmd->alpha=(real*) calloc(plmd->Jend*plmd->Nmemax,sizeof(real));
  plmd->beta=(real*) calloc(plmd->Jend*plmd->Nmemax,sizeof(real));

  plmd->hi=(real*) calloc(plmd->Jend,sizeof(real));
  plmd->x0=(real*) calloc(plmd->Jend,sizeof(real));
  plmd->dLdx0=(reala*) calloc(plmd->Jend,sizeof(reala));

  // for (i=0; i<plmd->Jend*plmd->Jend; i++) {
  //   plmd->Hinv[i]=0;
  // }
  // for (i=0; i<plmd->Jend; i++) {
  //   plmd->x[i]=0;
  //   plmd->Hinv[i*(plmd->Jend+1)]=1;
  // }

  /*{
    cudaResourceDesc resDesc;
    memset(&resDesc,0,sizeof(resDesc));
    resDesc.resType=cudaResourceTypeLinear;
    resDesc.res.linear.devPtr=plmd->h_d;
    resDesc.res.linear.desc=cudaCreateChannelDesc<real>();
    resDesc.res.linear.sizeInBytes=plmd->nblocks*sizeof(real);
    cudaTextureDesc texDesc;
    memset(&texDesc,0,sizeof(texDesc));
    texDesc.readMode=cudaReadModeElementType;
    cudaCreateTextureObject(&plmd->h_tex,&resDesc,&texDesc,NULL);
  }
  {
    cudaResourceDesc resDesc;
    memset(&resDesc,0,sizeof(resDesc));
    resDesc.resType=cudaResourceTypeLinear;
    resDesc.res.linear.devPtr=plmd->J_d;
    resDesc.res.linear.desc=cudaCreateChannelDesc<real>();
    resDesc.res.linear.sizeInBytes=plmd->nblocks*plmd->nblocks*sizeof(real);
    cudaTextureDesc texDesc;
    memset(&texDesc,0,sizeof(texDesc));
    texDesc.readMode=cudaReadModeElementType;
    cudaCreateTextureObject(&plmd->J_tex,&resDesc,&texDesc,NULL);
  }
  {
    cudaResourceDesc resDesc;
    memset(&resDesc,0,sizeof(resDesc));
    resDesc.resType=cudaResourceTypeLinear;
    resDesc.res.linear.devPtr=plmd->dhds_d;
    resDesc.res.linear.desc=cudaCreateChannelDesc<real>();
    resDesc.res.linear.sizeInBytes=plmd->nblocks*sizeof(real);
    cudaTextureDesc texDesc;
    memset(&texDesc,0,sizeof(texDesc));
    texDesc.readMode=cudaReadModeElementType;
    cudaCreateTextureObject(&plmd->dhds_tex,&resDesc,&texDesc,NULL);
  }
  {
    cudaResourceDesc resDesc;
    memset(&resDesc,0,sizeof(resDesc));
    resDesc.resType=cudaResourceTypeLinear;
    resDesc.res.linear.devPtr=plmd->dJds_d;
    resDesc.res.linear.desc=cudaCreateChannelDesc<real>();
    resDesc.res.linear.sizeInBytes=plmd->nblocks*plmd->nblocks*sizeof(real);
    cudaTextureDesc texDesc;
    memset(&texDesc,0,sizeof(texDesc));
    texDesc.readMode=cudaReadModeElementType;
    cudaCreateTextureObject(&plmd->dJds_tex,&resDesc,&texDesc,NULL);
  }*/

  return plmd;
}

void copyout_hJ(struct_plmd *plmd)
{
  int i,j,k;

  k=plmd->nblocks;
  for (i=0; i<plmd->nblocks; i++) {
    plmd->h[i]=plmd->x[i];
    for (j=plmd->block0[plmd->block2site[i]+1]; j<plmd->nblocks; j++) {
      plmd->J[i*plmd->nblocks+j]=plmd->x[k];
      plmd->J[j*plmd->nblocks+i]=plmd->x[k];
      k++;
    }
  }
}

void copyback_hJ(struct_plmd *plmd)
{
  int i,j,k;

  k=plmd->nblocks;
  for (i=0; i<plmd->nblocks; i++) {
    plmd->x[i]=plmd->h[i];
    for (j=plmd->block0[plmd->block2site[i]+1]; j<plmd->nblocks; j++) {
      plmd->x[k]=plmd->J[i*plmd->nblocks+j];
      k++;
    }
  }
}

void copyback_gradient(struct_plmd *plmd)
{
  int i,j,k;

  k=plmd->nblocks;
  for (i=0; i<plmd->nblocks; i++) {
    plmd->dLdx[i]=plmd->dLdh[i];
    for (j=plmd->block0[plmd->block2site[i]+1]; j<plmd->nblocks; j++) {
      plmd->dLdx[k]=(plmd->dLdJ[i*plmd->nblocks+j]+plmd->dLdJ[j*plmd->nblocks+i]);
      k++;
    }
  }
}

void copyout_gradient(struct_plmd *plmd)
{
  int i,j,k;

  k=plmd->nblocks;
  for (i=0; i<plmd->nblocks; i++) {
    plmd->dhds[i]=plmd->hi[i];
    for (j=plmd->block0[plmd->block2site[i]+1]; j<plmd->nblocks; j++) {
      plmd->dJds[i*plmd->nblocks+j]=plmd->hi[k];
      plmd->dJds[j*plmd->nblocks+i]=plmd->hi[k];
      k++;
    }
  }
}

#ifndef PLMGPU

void regularize_function(struct_plmd plmd)
{
  int i;
  reala L;

#ifdef PROFILE_CPU
  clock_t t1,t2;
  t1=clock();
#endif

  L=0.0;

  for (i=0; i<plmd.nblocks; i++) {
    L+=0.5*plmd.kh*plmd.h[i]*plmd.h[i];
  }

  for (i=0; i<plmd.nblocks*plmd.nblocks; i++) {
    L+=0.5*plmd.kJ*plmd.J[i]*plmd.J[i];
  }

  plmd.L[0]=L;

#ifdef PROFILE_CPU
  t2=clock();
  fprintf(stdout,"regularize_function time=%f\n",(double)(t2-t1)/(double)(CLOCKS_PER_SEC));
#endif
}

void regularize_function_line(real s,struct_plmd plmd)
{
  int i;
  reala L;
  reala dLds;

  L=0.0;
  dLds=0.0;

  for (i=0; i<plmd.nblocks; i++) {
    L+=0.5*plmd.kh*(plmd.h[i]+s*plmd.dhds[i])*(plmd.h[i]+s*plmd.dhds[i]);
    dLds+=plmd.kh*(plmd.h[i]+s*plmd.dhds[i])*plmd.dhds[i];
  }

  for (i=0; i<plmd.nblocks*plmd.nblocks; i++) {
    L+=0.5*plmd.kJ*(plmd.J[i]+s*plmd.dJds[i])*(plmd.J[i]+s*plmd.dJds[i]);
    dLds+=plmd.kJ*(plmd.J[i]+s*plmd.dJds[i])*plmd.dJds[i];
  }

  plmd.L[0]=L;
  plmd.dLds[0]=dLds;
}

void partition_function(int site,int block0,struct_plmd plmd)
{
  int b;
  int i,j;
  real heff;
  int nsubs;

  for (b=0; b<plmd.BiB; b++) {
    nsubs=plmd.nsubs[site];
    for (i=0; i<nsubs; i++) {
      heff=plmd.h[block0+i];
      for (j=0; j<plmd.nsites; j++) {
        heff+=plmd.J[plmd.nblocks*(block0+i) + plmd.Seq[b*plmd.nsites+j]];
      }
      plmd.heff[b*plmd.nblocks+block0+i]=heff;
    }
  }
}

void partition_function_line(int site,int block0,struct_plmd plmd)
{
  int b;
  int i,j;
  real heff;
  int nsubs;

  for (b=0; b<plmd.BiB; b++) {
    nsubs=plmd.nsubs[site];
    for (i=0; i<nsubs; i++) {
      heff=plmd.dhds[block0+i];
      for (j=0; j<plmd.nsites; j++) {
        heff+=plmd.dJds[plmd.nblocks*(block0+i) + plmd.Seq[b*plmd.nsites+j]];
      }
      plmd.dheffds[b*plmd.nblocks+block0+i]=heff;
    }
  }
}

void evaluate_function(int site,int block0,struct_plmd plmd)
{
  int b;
  int i;
  real Z;
  int nsubs;
  reala L;

  L=0.0;

  for (b=0; b<plmd.BiB; b++) {
    nsubs=plmd.nsubs[site];
    L+=-plmd.heff[b*plmd.nblocks+plmd.Seq[b*plmd.nsites+site]];
    Z=0.0;
    for (i=0; i<nsubs; i++) {
      Z+=exp(plmd.heff[b*plmd.nblocks+block0+i]);
    }
    L+=log(Z);
  }

  #pragma omp atomic
  plmd.L[0]+=L/plmd.B;
}

void evaluate_function_line(real s,int site,int block0,struct_plmd plmd)
{
  int b;
  int i;
  real m1,Z,w;
  int nsubs;
  reala dLds, L;

  L=0.0;
  dLds=0.0;

  for (b=0; b<plmd.BiB; b++) {
    nsubs=plmd.nsubs[site];
    L+=-plmd.heff[b*plmd.nblocks+plmd.Seq[b*plmd.nsites+site]];
    L+=-s*plmd.dheffds[b*plmd.nblocks+plmd.Seq[b*plmd.nsites+site]];
    dLds+=-plmd.dheffds[b*plmd.nblocks+plmd.Seq[b*plmd.nsites+site]];
    m1=0.0;
    Z=0.0;
    for (i=0; i<nsubs; i++) {
      w=exp(plmd.heff[b*plmd.nblocks+block0+i]+s*plmd.dheffds[b*plmd.nblocks+block0+i]);
      Z+=w;
      m1+=plmd.dheffds[b*plmd.nblocks+block0+i]*w;
    }
    L+=log(Z);
    dLds+=m1/Z;
  }

  #pragma omp atomic
  plmd.L[0]+=L/plmd.B;
  #pragma omp atomic
  plmd.dLds[0]+=dLds/plmd.B;
}

void regularize_gradient(struct_plmd plmd)
{
  int i;

  for (i=0; i<plmd.nblocks; i++) {
    plmd.dLdh[i]=plmd.kh*plmd.h[i];
  }
  for (i=0; i<plmd.nblocks*plmd.nblocks; i++) {
    plmd.dLdJ[i]=plmd.kJ*plmd.J[i];
  }
}

void evaluate_gradient_h(int site,int block0,struct_plmd plmd)
{
  int b;
  int i;
  int Seq;
  real Z;
  int nsubs;

  for (b=0; b<plmd.BiB; b++) {
    nsubs=plmd.nsubs[site];
    Seq=plmd.Seq[b*plmd.nsites+site]-block0;
    plmd.dLdh[block0+Seq]+=-1.0/plmd.B;
    Z=0;
    for (i=0; i<nsubs; i++) {
      Z+=exp(plmd.heff[b*plmd.nblocks+block0+i]);
    }
    for (i=0; i<nsubs; i++) {
      plmd.dLdh[block0+i]+=exp(plmd.heff[b*plmd.nblocks+block0+i])/Z/plmd.B;
    }
  }
}

void evaluate_gradient_J(int site1,int site2,int block01,int block02,struct_plmd plmd)
{
  int b;
  int i;
  int Seq1,Seq2;
  real Z;
  int nsubs1, nsubs2;

  nsubs1=plmd.nsubs[site1];
  nsubs2=plmd.nsubs[site2];
  for (b=0; b<plmd.BiB; b++) {
    Seq1=plmd.Seq[b*plmd.nsites+site1]-block01;
    Seq2=plmd.Seq[b*plmd.nsites+site2]-block02;
    plmd.dLdJ[plmd.nblocks*(block01+Seq1)+block02+Seq2]+=-1.0/plmd.B;
    Z=0;
    for (i=0; i<nsubs1; i++) {
      Z+=exp(plmd.heff[b*plmd.nblocks+block01+i]);
    }
    for (i=0; i<nsubs1; i++) {
      plmd.dLdJ[plmd.nblocks*(block01+i)+block02+Seq2]+=exp(plmd.heff[b*plmd.nblocks+block01+i])/Z/plmd.B;
    }
  }
}

#else

#define BLOCK 64
#define BATCH 16

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

__device__ static inline
int rectify_modulus(int a,int b)
{
  int c=a%b;
  c-=(c>=b?b:0);
  c+=(c<0?b:0);
  return c;
}

template <typename real_type>
__device__ inline
void reduce(real_type input,real_type *shared,real_type *global)
{
  real_type local=input;
  local+=__shfl_down_sync(0xFFFFFFFF,local,1);
  local+=__shfl_down_sync(0xFFFFFFFF,local,2);
  local+=__shfl_down_sync(0xFFFFFFFF,local,4);
  local+=__shfl_down_sync(0xFFFFFFFF,local,8);
  local+=__shfl_down_sync(0xFFFFFFFF,local,16);
  __syncthreads();
  if ((0x1F & threadIdx.x)==0) {
    shared[threadIdx.x>>5]=local;
  }
  __syncthreads();
  local=0;
  if (threadIdx.x < (blockDim.x>>5)) {
    local=shared[threadIdx.x];
  }
  if (threadIdx.x < 32) {
    if (blockDim.x>=64) local+=__shfl_down_sync(0xFFFFFFFF,local,1);
    if (blockDim.x>=128) local+=__shfl_down_sync(0xFFFFFFFF,local,2);
    if (blockDim.x>=256) local+=__shfl_down_sync(0xFFFFFFFF,local,4);
    if (blockDim.x>=512) local+=__shfl_down_sync(0xFFFFFFFF,local,8);
    if (blockDim.x>=1024) local+=__shfl_down_sync(0xFFFFFFFF,local,16);
  }
  if (threadIdx.x==0) {
    atomicAdd(global,local);
  }
}

__device__ inline
int sort(int i,int *sortBuffer)
{
  int i1,i2;
  int direction,otherThreadIdx,iother,bswitch;
  int itmp=i;
  int Ztmp=threadIdx.x;

  // Bitonic sort
  for (i1=1; i1<BLOCK; i1*=2) {
    direction=(((2*i1)&threadIdx.x)!=0); // 0 ascending, 1 descending
    for (i2=i1; i2>0; i2/=2) {
      otherThreadIdx=(threadIdx.x^i2);
      if (i2<32) {
        iother=__shfl_xor_sync(0xFFFFFFFF,itmp,i2);
        bswitch=(((otherThreadIdx>threadIdx.x)==(iother>itmp))==direction);
        bswitch=(iother==itmp?0:bswitch);
        itmp=__shfl_sync(0xFFFFFFFF,itmp,threadIdx.x^(i2*bswitch));
        Ztmp=__shfl_sync(0xFFFFFFFF,Ztmp,threadIdx.x^(i2*bswitch));
      } else {
        sortBuffer[threadIdx.x]=itmp;
        __syncthreads();
        iother=sortBuffer[otherThreadIdx];
        bswitch=(((otherThreadIdx>threadIdx.x)==(iother>itmp))==direction);
        bswitch=(iother==itmp?0:bswitch);
        itmp=sortBuffer[threadIdx.x^(i2*bswitch)];
        __syncthreads();
        sortBuffer[threadIdx.x]=Ztmp;
        __syncthreads();
        Ztmp=sortBuffer[threadIdx.x^(i2*bswitch)];
        __syncthreads();
      }
    }
  }
  // threadIdx.x is destination, Ztmp is source
  sortBuffer[Ztmp]=threadIdx.x;
  __syncthreads();
  Ztmp=sortBuffer[threadIdx.x];
  // Now threadIdx.x is source, and Ztmp is destination
  __syncthreads();
  sortBuffer[threadIdx.x]=itmp;
  __syncthreads();

  return Ztmp;
}

__device__ inline
void reduce_sorted(reala input,int *sortBuffer,int iDest,reala *sreduceBuffer,reala *shared)
{
  int i1;

  int itmp=sortBuffer[threadIdx.x];
  sreduceBuffer[iDest]=input;
  __syncthreads();

  // Reduction
  for (i1=1; i1<BLOCK; i1*=2) {
    if ((threadIdx.x&i1) && (threadIdx.x&(i1-1))==0) {
      if (itmp==sortBuffer[threadIdx.x-i1]) {
        sreduceBuffer[threadIdx.x-i1]+=sreduceBuffer[threadIdx.x];
      } else {
        shared[itmp]+=sreduceBuffer[threadIdx.x];
      }
    }
    __syncthreads();
  }
  if (threadIdx.x==0) {
    shared[itmp]+=sreduceBuffer[threadIdx.x];
  }
  __syncthreads();
}

__global__
void regularize_function(struct_plmd plmd)
{
  int t=blockIdx.x*blockDim.x+threadIdx.x;
  int i=t;
  reala L;
  __shared__ reala reduceBuffer[BLOCK>>5];

  L=0.0;

  if (i<plmd.nblocks) {
    L+=0.5*plmd.kh*plmd.h_d[i]*plmd.h_d[i];
    // real h=tex1Dfetch<real>(plmd.h_tex,i);
    // L+=0.5*plmd.kh*h*h;
  }

  if (i<plmd.nblocks*plmd.nblocks) {
    L+=0.5*plmd.kJ*plmd.J_d[i]*plmd.J_d[i];
    // real J=tex1Dfetch<real>(plmd.J_tex,i);
    // L+=0.5*plmd.kJ*J*J;
  }

  reduce(L,reduceBuffer,plmd.L_d);
}

__global__
void regularize_function_line(real s,struct_plmd plmd)
{
  int t=blockIdx.x*blockDim.x+threadIdx.x;
  int i=t;
  reala L;
  reala dLds;
  __shared__ reala reduceBuffer[BLOCK>>5];

  L=0.0;
  dLds=0.0;

  if (i<plmd.nblocks) {
    L+=0.5*plmd.kh*(plmd.h_d[i]+s*plmd.dhds_d[i])*(plmd.h_d[i]+s*plmd.dhds_d[i]);
    dLds+=plmd.kh*(plmd.h_d[i]+s*plmd.dhds_d[i])*plmd.dhds_d[i];
    // real h=tex1Dfetch<real>(plmd.h_tex,i);
    // real dhds=tex1Dfetch<real>(plmd.dhds_tex,i);
    // L+=0.5*plmd.kh*(h+s*dhds)*(h+s*dhds);
    // dLds+=plmd.kh*(h+s*dhds)*dhds;
  }

  if (i<plmd.nblocks*plmd.nblocks) {
    L+=0.5*plmd.kJ*(plmd.J_d[i]+s*plmd.dJds_d[i])*(plmd.J_d[i]+s*plmd.dJds_d[i]);
    dLds+=plmd.kJ*(plmd.J_d[i]+s*plmd.dJds_d[i])*plmd.dJds_d[i];
    // real J=tex1Dfetch<real>(plmd.J_tex,i);
    // real dJds=tex1Dfetch<real>(plmd.dJds_tex,i);
    // L+=0.5*plmd.kJ*(J+s*dJds)*(J+s*dJds);
    // dLds+=plmd.kJ*(J+s*dJds)*dJds;
  }

  reduce(L,reduceBuffer,plmd.L_d);
  reduce(dLds,reduceBuffer,plmd.dLds_d);
}

__global__
void partition_function(struct_plmd plmd)
{
  // int t=blockIdx.x*blockDim.x+threadIdx.x;
  int site=blockIdx.y;
  int ib,b;
  int i,j;
  real heff;
  int block0,nsubs;

  block0=plmd.block0_d[site];
  nsubs=plmd.nsubs_d[site];
  for (ib=0; ib<BATCH; ib++) {
    b=blockIdx.x*blockDim.x*BATCH+blockDim.x*ib+threadIdx.x;
    if (b<plmd.BiB) {
      for (i=0; i<nsubs; i++) {
        heff=plmd.h_d[block0+i];
        // heff=tex1Dfetch<real>(plmd.h_tex,block0+i);
        for (j=0; j<plmd.nsites; j++) {
          heff+=plmd.J_d[plmd.nblocks*(block0+i) + plmd.Seq_d[b*plmd.nsites+j]];
          // heff+=tex1Dfetch<real>(plmd.J_tex,plmd.nblocks*(block0+i) + plmd.Seq_d[b*plmd.nsites+j]);
        }
        plmd.heff_d[b*plmd.nblocks+block0+i]=heff;
      }
    }
  }
}

__global__
void partition_function_line(struct_plmd plmd)
{
  // int t=blockIdx.x*blockDim.x+threadIdx.x;
  int site=blockIdx.y;
  int ib,b;
  int i,j;
  real heff;
  int block0,nsubs;

  block0=plmd.block0_d[site];
  nsubs=plmd.nsubs_d[site];
  for (ib=0; ib<BATCH; ib++) {
    b=blockIdx.x*blockDim.x*BATCH+blockDim.x*ib+threadIdx.x;
    if (b<plmd.BiB) {
      for (i=0; i<nsubs; i++) {
        heff=plmd.dhds_d[block0+i];
        // heff=tex1Dfetch<real>(plmd.dhds_tex,block0+i);
        for (j=0; j<plmd.nsites; j++) {
          heff+=plmd.dJds_d[plmd.nblocks*(block0+i) + plmd.Seq_d[b*plmd.nsites+j]];
          // heff+=tex1Dfetch<real>(plmd.dJds_tex,plmd.nblocks*(block0+i) + plmd.Seq_d[b*plmd.nsites+j]);
        }
        plmd.dheffds_d[b*plmd.nblocks+block0+i]=heff;
      }
    }
  }
}

__global__
void evaluate_function(struct_plmd plmd)
{
  int site=blockIdx.y;
  int ib,b;
  int i;
  real Z;
  int block0,nsubs;
  reala L;
  __shared__ reala reduceBuffer[BLOCK>>5];

  L=0.0;

  block0=plmd.block0_d[site];
  nsubs=plmd.nsubs_d[site];
  for (ib=0; ib<BATCH; ib++) {
    b=blockIdx.x*blockDim.x*BATCH+blockDim.x*ib+threadIdx.x;
    if (b<plmd.BiB) {
      L+=-plmd.heff_d[b*plmd.nblocks+plmd.Seq_d[b*plmd.nsites+site]];
      Z=0.0;
      for (i=0; i<nsubs; i++) {
        Z+=exp(plmd.heff_d[b*plmd.nblocks+block0+i]);
      }
      L+=log(Z);
    }
  }

  L/=plmd.B;
  reduce(L,reduceBuffer,plmd.L_d);
}

__global__
void evaluate_function_line(real s,struct_plmd plmd)
{
  int site=blockIdx.y;
  int ib,b;
  int i;
  real m1,Z,w;
  int block0,nsubs;
  reala dLds, L;
  __shared__ reala reduceBuffer[BLOCK>>5];

  L=0.0;
  dLds=0.0;

  block0=plmd.block0_d[site];
  nsubs=plmd.nsubs_d[site];
  for (ib=0; ib<BATCH; ib++) {
    b=blockIdx.x*blockDim.x*BATCH+blockDim.x*ib+threadIdx.x;
    if (b<plmd.BiB) {
      L+=-plmd.heff_d[b*plmd.nblocks+plmd.Seq_d[b*plmd.nsites+site]];
      L+=-s*plmd.dheffds_d[b*plmd.nblocks+plmd.Seq_d[b*plmd.nsites+site]];
      dLds+=-plmd.dheffds_d[b*plmd.nblocks+plmd.Seq_d[b*plmd.nsites+site]];
      m1=0.0;
      Z=0.0;
      for (i=0; i<nsubs; i++) {
        w=exp(plmd.heff_d[b*plmd.nblocks+block0+i]+s*plmd.dheffds_d[b*plmd.nblocks+block0+i]);
        Z+=w;
        m1+=plmd.dheffds_d[b*plmd.nblocks+block0+i]*w;
      }
      L+=log(Z);
      dLds+=m1/Z;
    }
  }

  L/=plmd.B;
  dLds/=plmd.B;
  reduce(L,reduceBuffer,plmd.L_d);
  reduce(dLds,reduceBuffer,plmd.dLds_d);
}

__global__
void regularize_gradient(struct_plmd plmd)
{
  int i=blockIdx.x*blockDim.x+threadIdx.x;

  if (i<plmd.nblocks) {
    plmd.dLdh_d[i]=plmd.kh*plmd.h_d[i];
    // real h=tex1Dfetch<real>(plmd.h_tex,i);
    // plmd.dLdh_d[i]=plmd.kh*h;
  }
  if (i<plmd.nblocks*plmd.nblocks) {
    plmd.dLdJ_d[i]=plmd.kJ*plmd.J_d[i];
    // real J=tex1Dfetch<real>(plmd.J_tex,i);
    // plmd.dLdJ_d[i]=plmd.kJ*J;
  }
}

__global__
void evaluate_gradient_h(struct_plmd plmd)
{
  int site=blockIdx.y;
  int ib,b;
  int i;
  int Seq;
  real Z;
  int block0,nsubs;
  __shared__ reala reduceBuffer[BLOCK>>5];
  extern __shared__ reala dLdh_s[];
  reala dLdh;

  block0=plmd.block0_d[site];
  nsubs=plmd.nsubs_d[site];
  if (threadIdx.x<nsubs) {
    dLdh_s[threadIdx.x]=0;
  }
  for (ib=0; ib<BATCH; ib++) {
    b=blockIdx.x*blockDim.x*BATCH+blockDim.x*ib+threadIdx.x;
    if (b<plmd.BiB) Seq=plmd.Seq_d[b*plmd.nsites+site]-block0;
    Z=0;
    for (i=0; i<nsubs; i++) {
      if (b<plmd.BiB) Z+=exp(plmd.heff_d[b*plmd.nblocks+block0+i]);
    }
    for (i=0; i<nsubs; i++) {
      dLdh=0;
      if (b<plmd.BiB) dLdh=exp(plmd.heff_d[b*plmd.nblocks+block0+i])/Z-(i==Seq);
      reduce(dLdh,reduceBuffer,&dLdh_s[i]);
    }
  }

  __syncthreads();
  if (threadIdx.x<nsubs) {
    dLdh_s[threadIdx.x]/=plmd.B;
    atomicAdd(&plmd.dLdh_d[block0+threadIdx.x],dLdh_s[threadIdx.x]);
  }
}

__global__
void evaluate_gradient_J(struct_plmd plmd)
{
  int site1=blockIdx.y;
  int site2=blockIdx.z;
  int ib,b;
  int i;
  int Seq1,Seq2;
  int iDest;
  real Z;
  int block01,block02,nsubs1,nsubs2;
  __shared__ reala sreduceBuffer[BLOCK];
  __shared__ int sortBuffer[BLOCK];
  extern __shared__ reala dLdJ_s[];
  reala dLdJ;

  if (site1==site2) return;

  block01=plmd.block0_d[site1];
  block02=plmd.block0_d[site2];
  nsubs1=plmd.nsubs_d[site1];
  nsubs2=plmd.nsubs_d[site2];
  if (threadIdx.x<nsubs1*nsubs2) {
    dLdJ_s[threadIdx.x]=0;
  }
  for (ib=0; ib<BATCH; ib++) {
    b=blockIdx.x*blockDim.x*BATCH+blockDim.x*ib+threadIdx.x;
    Seq1=0;
    Seq2=0;
    if (b<plmd.BiB) Seq1=plmd.Seq_d[b*plmd.nsites+site1]-block01;
    if (b<plmd.BiB) Seq2=plmd.Seq_d[b*plmd.nsites+site2]-block02;
    // iSource=threadIdx.x;
    iDest=sort(Seq2,sortBuffer);
    Z=0;
    for (i=0; i<nsubs1; i++) {
      if (b<plmd.BiB) Z+=exp(plmd.heff_d[b*plmd.nblocks+block01+i]);
    }
    for (i=0; i<nsubs1; i++) {
      dLdJ=0;
      if (b<plmd.BiB) dLdJ=exp(plmd.heff_d[b*plmd.nblocks+block01+i])/Z-(i==Seq1);
      reduce_sorted(dLdJ,sortBuffer,iDest,sreduceBuffer,&dLdJ_s[i*nsubs2]);
    }
  }

  __syncthreads();
  if (threadIdx.x<nsubs1*nsubs2) {
    dLdJ_s[threadIdx.x]/=plmd.B;
    Seq1=threadIdx.x/nsubs2;
    Seq2=rectify_modulus(threadIdx.x,nsubs2);
    atomicAdd(&plmd.dLdJ_d[plmd.nblocks*(block01+Seq1)+block02+Seq2],dLdJ_s[threadIdx.x]);
  }
}
#endif

void evaluateL(struct_plmd *plmd)
{
  copyout_hJ(plmd);
#ifdef PLMGPU
  cudaMemcpy(plmd->h_d,plmd->h,plmd->nblocks*sizeof(real),cudaMemcpyHostToDevice);
  cudaMemcpy(plmd->J_d,plmd->J,plmd->nblocks*plmd->nblocks*sizeof(real),cudaMemcpyHostToDevice);

  cudaMemset(plmd->L_d,0,sizeof(reala));
  regularize_function<<<(plmd->nblocks*plmd->nblocks+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0]);
#else
  regularize_function(plmd[0]);
#endif

  if (plmd->B>plmd->BiBmax) {
    rewind(plmd->Seq_fp);
  }
  for (plmd->iB=0; plmd->iB<(plmd->B+plmd->BiBmax-1)/plmd->BiBmax; plmd->iB++) {
    plmd->BiB=plmd->B-plmd->iB*plmd->BiBmax;
    plmd->BiB=(plmd->BiB>plmd->BiBmax?plmd->BiBmax:plmd->BiB);
    if (plmd->B>plmd->BiBmax) {
      fread(plmd->Seq,sizeof(int),plmd->BiB*plmd->nsites,plmd->Seq_fp);
      cudaMemcpy(plmd->Seq_d,plmd->Seq,plmd->BiB*plmd->nsites*sizeof(int),cudaMemcpyHostToDevice);
    }
#ifdef PLMGPU
    {
      dim3 block_dim(BLOCK);
      dim3 grid_dim((plmd->BiB+BLOCK*BATCH+1)/(BLOCK*BATCH),plmd->nsites,1);
      partition_function<<<grid_dim,block_dim>>>(plmd[0]);
      evaluate_function<<<grid_dim,block_dim>>>(plmd[0]);
    }
#else
    #pragma omp parallel
    {
      int ID=omp_get_thread_num();
      int NID=omp_get_max_threads();
      for (int site=ID; site<plmd->nsites; site+=NID) {
        partition_function(site,plmd->block0[site],plmd[0]);
        evaluate_function(site,plmd->block0[site],plmd[0]);
      }
    }
#endif
  }

#ifdef PLMGPU
  cudaMemcpy(plmd->L,plmd->L_d,sizeof(reala),cudaMemcpyDeviceToHost);
#endif
}

void evaluateL_line(real s,struct_plmd *plmd)
{
#ifdef PLMGPU
  cudaMemset(plmd->L_d,0,sizeof(reala));
  cudaMemset(plmd->dLds_d,0,sizeof(reala));
  regularize_function_line<<<(plmd->nblocks*plmd->nblocks+BLOCK-1)/BLOCK,BLOCK>>>(s,plmd[0]);
#else
  regularize_function_line(s,plmd[0]);
#endif

  if (plmd->B>plmd->BiBmax) {
    rewind(plmd->Seq_fp);
  }
  for (plmd->iB=0; plmd->iB<(plmd->B+plmd->BiBmax-1)/plmd->BiBmax; plmd->iB++) {
    plmd->BiB=plmd->B-plmd->iB*plmd->BiBmax;
    plmd->BiB=(plmd->BiB>plmd->BiBmax?plmd->BiBmax:plmd->BiB);
    if (plmd->B>plmd->BiBmax) {
      fread(plmd->Seq,sizeof(int),plmd->BiB*plmd->nsites,plmd->Seq_fp);
#ifdef PLMGPU
      cudaMemcpy(plmd->Seq_d,plmd->Seq,plmd->BiB*plmd->nsites*sizeof(int),cudaMemcpyHostToDevice);
      dim3 block_dim(BLOCK);
      dim3 grid_dim((plmd->BiB+BLOCK*BATCH+1)/(BLOCK*BATCH),plmd->nsites,1);
      partition_function<<<grid_dim,block_dim>>>(plmd[0]);
      partition_function_line<<<grid_dim,block_dim>>>(plmd[0]);
#else
      #pragma omp parallel
      {
        int ID=omp_get_thread_num();
        int NID=omp_get_max_threads();
        for (int site=ID; site<plmd->nsites; site+=NID) {
          partition_function(site,plmd->block0[site],plmd[0]);
          partition_function_line(site,plmd->block0[site],plmd[0]);
        }
      }
#endif
    }
#ifdef PLMGPU
    {
      dim3 block_dim(BLOCK);
      dim3 grid_dim((plmd->BiB+BLOCK*BATCH+1)/(BLOCK*BATCH),plmd->nsites,1);
      evaluate_function_line<<<grid_dim,block_dim>>>(s,plmd[0]);
    }
#else
    #pragma omp parallel
    {
      int ID=omp_get_thread_num();
      int NID=omp_get_max_threads();
      for (int site=ID; site<plmd->nsites; site+=NID) {
        evaluate_function_line(s,site,plmd->block0[site],plmd[0]);
      }
    }
#endif
  }

#ifdef PLMGPU
  cudaMemcpy(plmd->L,plmd->L_d,sizeof(reala),cudaMemcpyDeviceToHost);
  cudaMemcpy(plmd->dLds,plmd->dLds_d,sizeof(reala),cudaMemcpyDeviceToHost);
#endif
}

void evaluatedLdx(struct_plmd *plmd)
{
#ifdef PLMGPU
  regularize_gradient<<<(plmd->nblocks*plmd->nblocks+BLOCK-1)/BLOCK,BLOCK>>>(plmd[0]);
#else
  regularize_gradient(plmd[0]);
#endif

  if (plmd->B>plmd->BiBmax) {
    rewind(plmd->Seq_fp);
  }
  for (plmd->iB=0; plmd->iB<(plmd->B+plmd->BiBmax-1)/plmd->BiBmax; plmd->iB++) {
    plmd->BiB=plmd->B-plmd->iB*plmd->BiBmax;
    plmd->BiB=(plmd->BiB>plmd->BiBmax?plmd->BiBmax:plmd->BiB);
    if (plmd->B>plmd->BiBmax) {
      fread(plmd->Seq,sizeof(int),plmd->BiB*plmd->nsites,plmd->Seq_fp);
#ifdef PLMGPU
      cudaMemcpy(plmd->Seq_d,plmd->Seq,plmd->BiB*plmd->nsites*sizeof(int),cudaMemcpyHostToDevice);
      dim3 block_dim(BLOCK);
      dim3 grid_dim((plmd->BiB+BLOCK*BATCH+1)/(BLOCK*BATCH),plmd->nsites,1);
      partition_function<<<grid_dim,block_dim>>>(plmd[0]);
#else
      #pragma omp parallel
      {
        int ID=omp_get_thread_num();
        int NID=omp_get_max_threads();
        for (int site=ID; site<plmd->nsites; site+=NID) {
          partition_function(site,plmd->block0[site],plmd[0]);
        }
      }
#endif
    }
#ifdef PLMGPU
    {
      dim3 block_dim(BLOCK);
      dim3 grid_dim((plmd->BiB+BLOCK*BATCH+1)/(BLOCK*BATCH),plmd->nsites,1);
      evaluate_gradient_h<<<grid_dim,block_dim>>>(plmd[0]);
    }
    {
      dim3 block_dim(BLOCK);
      dim3 grid_dim((plmd->BiB+BLOCK*BATCH+1)/(BLOCK*BATCH),plmd->nsites,plmd->nsites);
      int shmem=plmd->nsubsmax*plmd->nsubsmax*sizeof(reala);
      evaluate_gradient_J<<<grid_dim,block_dim,shmem>>>(plmd[0]);
    }
#else
    #pragma omp parallel
    {
      int ID=omp_get_thread_num();
      int NID=omp_get_max_threads();
      for (int site1=ID; site1<plmd->nsites; site1+=NID) {
        int block01=plmd->block0[site1];
        evaluate_gradient_h(site1,block01,plmd[0]);
        for (int site2=0; site2<plmd->nsites; site2++) {
          int block02=plmd->block0[site2];
          if (site1 != site2) {
            evaluate_gradient_J(site1,site2,block01,block02,plmd[0]);
          }
        }
      }
    }
#endif
  }

#ifdef PLMGPU
  cudaMemcpy(plmd->dLdh,plmd->dLdh_d,plmd->nblocks*sizeof(reala),cudaMemcpyDeviceToHost);
  cudaMemcpy(plmd->dLdJ,plmd->dLdJ_d,plmd->nblocks*plmd->nblocks*sizeof(reala),cudaMemcpyDeviceToHost);
#endif

  copyback_gradient(plmd);
}

void resetHinv(struct_plmd *plmd)
{
  int i;
  // N^2 Hinv
  // for (i=0; i<plmd->Jend*plmd->Jend; i++) {
  //   plmd->Hinv[i]=0.0;
  // }
  // for (i=0; i<plmd->Jend; i++) {
  //   plmd->Hinv[(plmd->Jend+1)*i]=1.0;
  // }
  // N^1 Hinv
  for (i=0; i<plmd->Jend; i++) {
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

  DxDg=0.0;
  for (i=0;i<plmd->Jend;i++) {
    // Put Delta x and Delta dLdx in x0 and dLdx0, which hold previous values
    plmd->x0[i]=plmd->x[i]-plmd->x0[i];
    plmd->dLdx0[i]=plmd->dLdx[i]-plmd->dLdx0[i];
    DxDg+=plmd->x0[i]*plmd->dLdx0[i];
  }

  DgHinvDg=0.0;
  for (i=0; i<plmd->Jend; i++) {
    plmd->hi[i]=0.0;
    for (j=0; j<plmd->Jend; j++) {
      // put Hinv * Delta dLdx in hi (the search direction) as a buffer
      plmd->hi[i]+=plmd->Hinv[i*plmd->Jend+j]*plmd->dLdx0[j];
    }
    DgHinvDg+=plmd->hi[i]*plmd->dLdx0[i];
  }
  c1=(1.0+DgHinvDg/DxDg)/DxDg;
  c2=-1.0/DxDg;

  for (i=0; i<plmd->Jend; i++) {
    for (j=i; j<plmd->Jend; j++) {
      plmd->Hinv[i*plmd->Jend+j]+=c1*(plmd->x0[i]*plmd->x0[j])+c2*(plmd->hi[i]*plmd->x0[j]+plmd->x0[i]*plmd->hi[j]);
      plmd->Hinv[j*plmd->Jend+i]=plmd->Hinv[i*plmd->Jend+j];
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
    for (j=0; j<plmd->Jend; j++) {
      plmd->d_x[i*plmd->Jend+j]=plmd->d_x[(i-1)*plmd->Jend+j];
      plmd->d_dLdx[i*plmd->Jend+j]=plmd->d_dLdx[(i-1)*plmd->Jend+j];
    }
    plmd->rho[i]=plmd->rho[i-1];
  }

  plmd->rho[0]=0;
  for (i=0; i<plmd->Jend; i++) {
    plmd->d_x[i]=plmd->x[i]-plmd->x0[i];
    plmd->d_dLdx[i]=plmd->dLdx[i]-plmd->dLdx0[i];
    plmd->rho[0]+=plmd->d_x[i]*plmd->d_dLdx[i];
  }
  plmd->rho[0]=1.0/plmd->rho[0];

  for (i=0; i<plmd->Jend; i++) {
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
  for (i=0; i<plmd->Jend; i++) {
    plmd->hi[i]=0;
    for (j=0; j<plmd->Jend; j++) {
      plmd->hi[i]+=-plmd->Hinv[i*plmd->Jend+j]*plmd->dLdx[j];
    }
  }

  dLds=0;
  for (i=0; i<plmd->Jend; i++) {
    dLds+=plmd->hi[i]*plmd->dLdx[i];
  }

  if (dLds>0) {
    fprintf(stderr,"Bad direction, reset Hinv\n");
    for (i=0; i<plmd->Jend*plmd->Jend; i++) {
      plmd->Hinv[i]=0;
    }
    for (i=0; i<plmd->Jend; i++) {
      plmd->Hinv[i*(plmd->Jend+1)]=1;
      plmd->hi[i]=plmd->dLdx[i];
    }
  }
*/ // End N^2 Hinv
// Begin N^1 Hinv
  int i,j;

  for (i=0; i<plmd->Jend; i++) {
    plmd->hi[i]=plmd->dLdx[i];
  }
  for (i=0; i<plmd->Nmem; i++) {
    plmd->alpha[i]=0;
    for (j=0; j<plmd->Jend; j++) {
      plmd->alpha[i]+=plmd->d_x[i*plmd->Jend+j]*plmd->hi[j];
    }
    plmd->alpha[i]*=plmd->rho[i];
    for (j=0; j<plmd->Jend; j++) {
      plmd->hi[j]+=-plmd->alpha[i]*plmd->d_dLdx[i*plmd->Jend+j];
    }
  }
  /*
  // According to wikipedia, this is to ensure the step length is always about unity
  // https://en.wikipedia.org/wiki/Limited-memory_BFGS
  if (plmd->Nmem>0) {
    numer=0.0L
    denom=0.0L;
    for (i=0; i<plmd->Jend; i++) {
      numer+=plmd->d_x[i]*plmd->hi[i];
      denom+=plmd->d_dLdx[i]*plmd->d_dLdx[i];
    }
    numer/=denom;
    for (i=0; i<plmd->Jend; i++) {
      plmd->hi[i]=numer*plmd->d_dLdx[i]; // This seems like a horrible idea, maybe wikipedia has a typo...
    }
  }
  */
  for (i=plmd->Nmem-1; i>=0; i--) {
    plmd->beta[i]=0;
    for (j=0; j<plmd->Jend; j++) {
      plmd->beta[i]+=plmd->d_dLdx[i*plmd->Jend+j]*plmd->hi[j];
    }
    plmd->beta[i]*=plmd->rho[i];
    for (j=0; j<plmd->Jend; j++) {
      plmd->hi[j]+=(plmd->alpha[i]-plmd->beta[i])*plmd->d_x[i*plmd->Jend+j];
    }
  }

  for (i=0; i<plmd->Jend; i++) {
    plmd->hi[i]*=-1;
  }
// End N^1 Hinv

  // New stuff
  copyout_gradient(plmd);
#ifdef PLMGPU
  cudaMemcpy(plmd->dhds_d,plmd->dhds,plmd->nblocks*sizeof(real),cudaMemcpyHostToDevice);
  cudaMemcpy(plmd->dJds_d,plmd->dJds,plmd->nblocks*plmd->nblocks*sizeof(real),cudaMemcpyHostToDevice);
#endif

  if (plmd->B>plmd->BiBmax) {
    rewind(plmd->Seq_fp);
  }
  for (plmd->iB=0; plmd->iB<(plmd->B+plmd->BiBmax-1)/plmd->BiBmax; plmd->iB++) {
    plmd->BiB=plmd->B-plmd->iB*plmd->BiBmax;
    plmd->BiB=(plmd->BiB>plmd->BiBmax?plmd->BiBmax:plmd->BiB);
    if (!(plmd->B>plmd->BiBmax)) {
#ifdef PLMGPU
      dim3 block_dim(BLOCK);
      dim3 grid_dim((plmd->BiB+BLOCK*BATCH+1)/(BLOCK*BATCH),plmd->nsites,1);
      partition_function_line<<<grid_dim,block_dim>>>(plmd[0]);
#else
      #pragma omp parallel
      {
        int ID=omp_get_thread_num();
        int NID=omp_get_max_threads();
        for (int site=ID; site<plmd->nsites; site+=NID) {
          partition_function_line(site,plmd->block0[site],plmd[0]);
        }
      }
#endif
    }
  }
}

void update_line(int step,struct_plmd *plmd)
{
  int i;
  reala a,b,c,s;
  real s1,s2,s3;
  reala L1,L2,L3;
  reala dLds1,dLds2,dLds3;
  reala L0;

  for (i=0; i<plmd->Jend; i++) {
    plmd->x0[i]=plmd->x[i];
    plmd->dLdx0[i]=plmd->dLdx[i];
  }

  L0=plmd->L[0];

  s1=0.0;
  evaluateL_line(s1,plmd);
  L1=plmd->L[0];
  dLds1=plmd->dLds[0];
  if (dLds1>0) {
    fprintf(stderr,"Error, hi is pointing wrong way\n");
    exit(1);
  }
  
  s3=1.0;
  evaluateL_line(s3,plmd);
  L3=plmd->L[0];
  dLds3=plmd->dLds[0];

  while (dLds3<0 && s3<100000000L) {
    s2=s1-dLds1*(s3-s1)/(dLds3-dLds1);
    s3=(1.5*s2>8*s3 ? 8*s3 : 1.5*s2); // s2 is expected 0. Go past it by 50%, unless that's an increase of more than a factor of 8.
    evaluateL_line(s3,plmd);
    L3=plmd->L[0];
    dLds3=plmd->dLds[0];
  }

  s2=s1-dLds1*(s3-s1)/(dLds3-dLds1);
  evaluateL_line(s2,plmd);
  L2=plmd->L[0];
  dLds2=plmd->dLds[0];

  fprintf(stderr,"Step %4d s=%lg %lg %lg\n          L=%lg %lg %lg\n       dLds=%lg %lg %lg\n",
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

    fprintf(stderr,"Step %4d s=%lg %lg %lg\n          L=%lg %lg %lg\n       dLds=%lg %lg %lg\n",
            step,(double) s1,(double) s2,(double) s3,
            (double) L1,(double) L2,(double) L3,
            (double) dLds1,(double) dLds2,(double) dLds3);
  }

  fprintf(stderr,"Step %4d s=%lg %lg %lg\n          L=%lg %lg %lg\n       dLds=%lg %lg %lg\n",
          step,(double) s1,(double) s2,(double) s3,
          (double) L1,(double) L2,(double) L3,
          (double) dLds1,(double) dLds2,(double) dLds3);

  real stepLength2=0;

  for (i=0; i<plmd->Jend; i++) {
    plmd->x[i]=plmd->x0[i]+s2*plmd->hi[i];
    stepLength2+=(s2*plmd->hi[i])*(s2*plmd->hi[i]);
  }

  // fprintf(stderr,"Step %d smid1 %lg smid2 %lg L=%lg -> L1=%lg L2=%lg\n",step,(double) smid1,(double) smid2,(double) L0,(double) Lmid1,(double) Lmid2);
  fprintf(stderr,"Step %4d L=%24.16lf -> L2=%24.16lf, dL=%lg, step length=%lg\n",step,(double)L0,(double)L2,(double)(L2-L0),(double) sqrt(stepLength2));

  if (sqrt(stepLength2)<5e-7) plmd->done=true;
}

real lineL(real s,struct_plmd *plmd)
{
  int i;

  for (i=0; i<plmd->Jend; i++) {
    plmd->x[i]=plmd->x0[i]+s*plmd->hi[i];
  }

  evaluateL(plmd);
  return plmd->L[0];
}

void update(int step,struct_plmd *plmd)
{
  real smin,smid1,smid2,smax;
  real Lmin,Lmid1,Lmid2,Lmax;
  real L0;
  int i,ss;
  real phi=(1.0L+sqrtl(5.0L))/2.0L;

  for (i=0; i<plmd->Jend; i++) {
    plmd->x0[i]=plmd->x[i];
    plmd->dLdx0[i]=plmd->dLdx[i];
  }
  smin=0.0L;
  smax=1.0L;

  // evaluateL(plmd);
  Lmin=plmd->L[0];
  L0=Lmin;

  Lmax=lineL(smax,plmd);
  // fprintf(stderr,"Step %d smax %g L=%g -> L=%g\n",step,smax,L0,Lmax);

  while (Lmax<Lmin && smax<100000000L) {
    smax*=2.0L;

    Lmax=lineL(smax,plmd);
    // fprintf(stderr,"Step %d smax %g L=%g -> L=%g\n",step,smax,L0,Lmax);
  }

  smid1=(smax-smin)/(phi*phi)+smin;
  Lmid1=lineL(smid1,plmd);

  smid2=-1.0L;
  
  for (ss=0; ss<25; ss++) {
    if (smid1<0.0L) {
      smid1=(smid2-smin)/phi+smin;
      Lmid1=lineL(smid1,plmd);
    } else {
      smid2=smax-(smax-smid1)/phi;
      Lmid2=lineL(smid2,plmd);
    }
    // fprintf(stderr,"Step %d smid1 %g smid2 %g L=%g -> L1=%g L2=%g\n",step,smid1,smid2,L0,Lmid1,Lmid2);
    fprintf(stderr,"Step %d smid1 %lg Lmin=%lg dL=%lg %lg %lg %lg\n",step,(double) smid1,(double) Lmin,(double) (Lmin-Lmin),(double) (Lmid1-Lmin),(double) (Lmid2-Lmin),(double) (Lmax-Lmin));

    if (Lmid1<=Lmid2) {
      smax=smid2;
      Lmax=Lmid2;
      smid2=smid1;
      Lmid2=Lmid1;
      smid1=-1.0L;
    } else {
      smin=smid1;
      Lmin=Lmid1;
      smid1=smid2;
      Lmid1=Lmid2;
      smid2=-1.0L;
    }
  }
  fprintf(stderr,"Step %d smid1 %lg smid2 %lg L=%lg -> L1=%lg L2=%lg\n",step,(double) smid1,(double) smid2,(double) L0,(double) Lmid1,(double) Lmid2);
}

/*
void lineL(double s,double *L,double *dLds,struct_plmd *plmd)
{
  int i;

  for (i=0; i<plmd->Jend; i++) {
    plmd->x[i]=plmd->x0[i]+s*plmd->hi[i];
  }

  evaluateL(plmd);
  evaluatedLdx(plmd);

  L[0]=plmd->L[0];

  dLds[0]=0;
  for (i=0; i<plmd->Jend; i++) {
    dLds[0]+=plmd->dLdx[i]*plmd->hi[i];
  }
}

void update(int step,struct_plmd *plmd)
{
  double smin,smid,smax;
  double Lmin,Lmid,Lmax;
  double dLdsmin,dLdsmid,dLdsmax;
  double L0;
  double d2Lds2;
  int i,ss;

  for (i=0; i<plmd->Jend; i++) {
    plmd->x0[i]=plmd->x[i];
    plmd->dLdx0[i]=plmd->dLdx[i];
  }
  smin=0;
  smax=1;

  // evaluateL(plmd);
  Lmin=plmd->L[0];
  L0=Lmin;

  dLdsmin=0;
  for (i=0; i<plmd->Jend; i++) {
    dLdsmin+=plmd->dLdx[i]*plmd->hi[i];
  }

  lineL(smax,&Lmax,&dLdsmax,plmd);
  fprintf(stderr,"Step %d srange %g %g Lrange %g %g dLdsrange %g %g\n",step,smin,smax,Lmin,Lmax,dLdsmin,dLdsmax);

  for (ss=0; ss<10; ss++) {
    d2Lds2=(dLdsmax-dLdsmin)/(smax-smin);
    smid=smin-dLdsmin/d2Lds2;
    lineL(smid,&Lmid,&dLdsmid,plmd);

    if (dLdsmid>0) {
      smax=smid;
      Lmax=Lmid;
      dLdsmax=dLdsmid;
    } else if (dLdsmid<0) {
      smin=smid;
      Lmin=Lmid;
      dLdsmin=dLdsmid;
    } else {
      break;
    }
    fprintf(stderr,"Step %d srange %g %g Lrange %g %g dLdsrange %g %g\n",step,smin,smax,Lmin,Lmax,dLdsmin,dLdsmax);
  }
  fprintf(stderr,"Step %d smid %g L=%g -> L=%g\n",step,smid,L0,Lmid);
}
*/
/*
void initialize(struct_plmd *plmd)
{
  evaluateL(plmd);
  evaluatedLdx(plmd);

  projectHinv(plmd);
  update(-1,plmd);
}
*/

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

  //initialize(plmd);
  plmd->done=false;
  for (s=0; s<1000; s++) {
    itterate(s,plmd);
    if (plmd->done) break;
  }
}

void finish(struct_plmd *plmd,int argc, char *argv[])
{
  int i,j;
  FILE *fp;

  fp=fopen(argv[1],"w");
  for (i=0; i<plmd->nblocks; i++) {
    fprintf(fp," %lg",(double) plmd->h[i]);
  }
  fclose(fp);

  fp=fopen(argv[2],"w");
  for (i=0; i<plmd->nblocks; i++) {
    for (j=0; j<plmd->nblocks; j++) {
      fprintf(fp," %lg",(double) plmd->J[i*plmd->nblocks+j]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

  free(plmd->nsubs);
  free(plmd->block2site);
  free(plmd->Seq);
  free(plmd->heff);

  free(plmd->x);
  free(plmd->dLdx);
  // free(plmd->Hinv);
  free(plmd->d_x);
  free(plmd->d_dLdx);
  free(plmd->rho);
  free(plmd->alpha);
  free(plmd->beta);
  free(plmd->hi);
  free(plmd->x0);
  free(plmd->dLdx0);

  // fclose(plmd->Seq_fp);

  /*cudaDestroyTextureObject(plmd->h_tex);
  cudaDestroyTextureObject(plmd->J_tex);
  cudaDestroyTextureObject(plmd->dhds_tex);
  cudaDestroyTextureObject(plmd->dJds_tex);*/
}

int main(int argc, char *argv[])
{
  struct_plmd *plmd;

  plmd = setup(argc,argv);
 
  run(plmd);

  finish(plmd,argc,argv);

  return 0;
}
