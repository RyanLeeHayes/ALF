// Written by Ryan Hayes 2017-06-20
// plmDCA algorithm from R470 - DOI: 10.1016/j.jcp.2014.07.024
// Quasi newton equations from https://www.rose-hulman.edu/~bryan/lottamath/quasinewton.pdf

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

//   ID=omp_get_thread_num();
//   NID=omp_get_max_threads();

#include "LM.h"

#define MAXLENGTH 1024
#define BLOCK 256
#define MPI_REALTYPE MPI_DOUBLE

struct_plmd* setup(int argc, char *argv[])
{
  struct_plmd *plmd;
  int i,j,k,l;
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
  for(i=0; i<plmd->nsites; i++) {
    fscanf(fp,"%d",&(plmd->nsubs[i]));
    plmd->nsubs[i]++;
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

  plmd->h=(real*) calloc(plmd->nblocks,sizeof(real));
  plmd->J=(real*) calloc(plmd->nblocks*plmd->nblocks,sizeof(real));

  plmd->dLdh=(real*) calloc(plmd->nblocks,sizeof(real));
  plmd->dLdJ=(real*) calloc(plmd->nblocks*plmd->nblocks,sizeof(real));

  plmd->kh=1e-6;
  plmd->kJ=1e-6;

  plmd->L=(real*) calloc(1,sizeof(real));


  k=plmd->nblocks;
  for (i=0; i<plmd->nblocks; i++) {
    for (j=plmd->block0[plmd->block2site[i]+1]; j<plmd->nblocks; j++) {
      k++;
    }
  }
  plmd->Jend=k;

  plmd->x=(real*) calloc(plmd->Jend,sizeof(real));
  plmd->dLdx=(real*) calloc(plmd->Jend,sizeof(real));
  plmd->Hinv=(real*) calloc(plmd->Jend*plmd->Jend,sizeof(real));

  plmd->hi=(real*) calloc(plmd->Jend,sizeof(real));
  plmd->x0=(real*) calloc(plmd->Jend,sizeof(real));
  plmd->dLdx0=(real*) calloc(plmd->Jend,sizeof(real));

  for (i=0; i<plmd->Jend*plmd->Jend; i++) {
    plmd->Hinv[i]=0;
  }
  for (i=0; i<plmd->Jend; i++) {
    plmd->x[i]=0;
    plmd->Hinv[i*(plmd->Jend+1)]=1;
  }

  for(i=0; i<plmd->nblocks; i++) {
    plmd->h[i]=0;
  }

  for(i=0; i<plmd->nblocks; i++) {
    for(j=0; j<plmd->nblocks; j++) {
      plmd->J[i*plmd->nblocks+j]=0;
    }
  }

  plmd->Z=(real**) calloc(omp_get_max_threads()+1,sizeof(real*));
  plmd->m1=(real**) calloc(omp_get_max_threads()+1,sizeof(real*));
  plmd->m2=(real**) calloc(omp_get_max_threads()+1,sizeof(real*));
  for (i=0; i<omp_get_max_threads()+1; i++) {
    plmd->Z[i]=(real*) calloc(1,sizeof(real));
    plmd->m1[i]=(real*) calloc(plmd->nblocks,sizeof(real));
    plmd->m2[i]=(real*) calloc(plmd->nblocks*plmd->nblocks,sizeof(real));
  }

  plmd->m1obs=(real*) calloc(plmd->nblocks,sizeof(real));
  plmd->m2obs=(real*) calloc(plmd->nblocks*plmd->nblocks,sizeof(real));

  fp=fopen(argv[3],"r");
  if (fp==NULL) {
    fprintf(stderr,"Error: %s does not exist\n",argv[3]);
    exit(1);
  }
  for(i=0; i<plmd->nblocks; i++) {
    fscanf(fp,"%lg\n",&plmd->m1obs[i]);
  }
  fclose(fp);

  fp=fopen(argv[4],"r");
  if (fp==NULL) {
    fprintf(stderr,"Error: %s does not exist\n",argv[4]);
    exit(1);
  }
  for(i=0; i<plmd->nblocks; i++) {
    for(j=0; j<plmd->nblocks; j++) {
      fscanf(fp,"%lg\n",&plmd->m2obs[i*plmd->nblocks+j]);
    }
  }
  fclose(fp);

  return plmd;
}

void regularize_function(struct_plmd plmd)
{
  int i;
  real L;
  int mID;

  MPI_Comm_rank(MPI_COMM_WORLD, &mID);

  if (mID==0) {
  L=0.0L;

  for (i=0; i<plmd.nblocks; i++) {
    L+=0.5L*plmd.kh*plmd.h[i]*plmd.h[i];
  }

  for (i=0; i<plmd.nblocks*plmd.nblocks; i++) {
    L+=0.5L*plmd.kJ*plmd.J[i]*plmd.J[i];
  }

  plmd.L[0]=L;
  } else {
  plmd.L[0]=0.0L;
  }
}

void runZ(struct_plmd *plmd) {
  long long int s;
  long long int prod;
  int Seq[plmd->nsites];
  int i,j;
  real H,w;
  real H0,HC;
  int ID,NID,mID,mNID,oID,oNID;
  long long int smin,smax;
  char line[MAXLENGTH];
  int offset;

  oID=omp_get_thread_num();
  oNID=omp_get_max_threads();

  MPI_Comm_rank(MPI_COMM_WORLD, &mID);
  MPI_Comm_size(MPI_COMM_WORLD, &mNID);

  ID=oNID*mID+oID;
  NID=mNID*oNID;

  prod=1;
  for (i=0; i<plmd->nsites; i++) {
    prod*=plmd->nsubs[i];
  }

  smin=(prod * ID)/NID;
  smax=(prod * (ID+1))/NID;

  s=smin;
  offset=0;
  for (i=0; i<plmd->nsites; i++) {
    Seq[i]=plmd->block0[i]+(s%plmd->nsubs[i]);
    s/=plmd->nsubs[i];
    offset+=sprintf(line+offset," %d",Seq[i]);
  }
  // fprintf(stderr,"%12ld to %12ld on %d\nStarting sequence %s on %d\n",smin,smax,ID,line,ID);

  plmd->Z[oID][0]=0.0L;
  for (i=0; i<plmd->nblocks; i++) {
    plmd->m1[oID][i]=0.0L;
    for (j=0; j<plmd->nblocks; j++) {
      plmd->m2[oID][i*plmd->nblocks+j]=0.0L;
    }
  }

  // for (s=0; s<prod; s++)
  for (s=smin; s<smax; s++) {
    #pragma omp master
    {
    if (s%1000000==0 && mID==0) {
      // fprintf(stderr,"%12ld sequences scanned by master of %12ld\n",s,smax);
    }
    }

    HC=plmd->h[Seq[0]];
    for (i=1; i<plmd->nsites; i++) {
      HC+=0.5*plmd->J[Seq[0]*plmd->nblocks+Seq[i]];
      HC+=0.5*plmd->J[Seq[i]*plmd->nblocks+Seq[0]];
    }

    if (Seq[0]==0 || s==smin) {
      H0=0;
      for (i=1; i<plmd->nsites; i++) {
        H0+=plmd->h[Seq[i]];
        for (j=1; j<plmd->nsites; j++) {
          H0+=0.5*plmd->J[Seq[i]*plmd->nblocks+Seq[j]];
        }
      }
    }
    H=H0+HC;
    w=exp(H);

    plmd->Z[oID][0]+=w;
    for (i=0; i<plmd->nsites; i++) {
      plmd->m1[oID][Seq[i]]+=w;
      for (j=0; j<plmd->nsites; j++) {
        plmd->m2[oID][Seq[i]*plmd->nblocks+Seq[j]]+=w;
      }
    }

    Seq[0]+=1;
    for (i=0; i<plmd->nsites-1; i++) {
      if (Seq[i]==plmd->block0[i+1]) {
        Seq[i]=plmd->block0[i];
        Seq[i+1]+=1;
      } else {
        break;
      }
    }
  }
}

void reduceZ(struct_plmd *plmd)
{
  int ID,NID,i,j;

  NID=omp_get_max_threads();

  for (ID=1; ID<NID; ID++) {
    plmd->Z[0][0]+=plmd->Z[ID][0];
    for (i=0; i<plmd->nblocks; i++) {
      plmd->m1[0][i]+=plmd->m1[ID][i];
      for (j=0; j<plmd->nblocks; j++) {
        plmd->m2[0][i*plmd->nblocks+j]+=plmd->m2[ID][i*plmd->nblocks+j];
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Reduce(plmd->Z[0],plmd->Z[NID],1,MPI_REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(plmd->m1[0],plmd->m1[NID],plmd->nblocks,MPI_REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(plmd->m2[0],plmd->m2[NID],plmd->nblocks*plmd->nblocks,MPI_REALTYPE,MPI_SUM,0,MPI_COMM_WORLD);
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

  MPI_Bcast(plmd->h,plmd->nblocks,MPI_REALTYPE,0,MPI_COMM_WORLD);
  MPI_Bcast(plmd->J,plmd->nblocks*plmd->nblocks,MPI_REALTYPE,0,MPI_COMM_WORLD);
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

  MPI_Bcast(plmd->x,k,MPI_REALTYPE,0,MPI_COMM_WORLD);
}

void regularize_gradient(struct_plmd plmd)
{
  int i;
  int mID;

  MPI_Comm_rank(MPI_COMM_WORLD, &mID);

  if (mID==0) {
  for (i=0; i<plmd.nblocks; i++) {
    plmd.dLdh[i]=plmd.kh*plmd.h[i];
  }
  for (i=0; i<plmd.nblocks*plmd.nblocks; i++) {
    plmd.dLdJ[i]=plmd.kJ*plmd.J[i];
  }
  } else {
  for (i=0; i<plmd.nblocks; i++) {
    plmd.dLdh[i]=0.0L;
  }
  for (i=0; i<plmd.nblocks*plmd.nblocks; i++) {
    plmd.dLdJ[i]=0.0L;
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

  MPI_Bcast(plmd->dLdx,k,MPI_REALTYPE,0,MPI_COMM_WORLD);
}

void evaluateL(struct_plmd *plmd)
{
  int i,j;
  int mID,oNID;

  oNID=omp_get_max_threads();
  MPI_Comm_rank(MPI_COMM_WORLD, &mID);

  copyout_hJ(plmd);

  #pragma omp parallel
  {
  runZ(plmd);
  }
  reduceZ(plmd);

  regularize_function(plmd[0]);

  if (mID==0) {
  plmd->L[0]+=logl(plmd->Z[oNID][0]);
  for (i=0; i<plmd->nblocks; i++) {
    plmd->L[0]+=-plmd->m1obs[i]*plmd->h[i];
    for (j=plmd->block0[plmd->block2site[i]+1]; j<plmd->nblocks; j++) {
      plmd->L[0]+=-plmd->m2obs[i*plmd->nblocks+j]*plmd->J[i*plmd->nblocks+j];
    }
  }
  }

  MPI_Bcast(plmd->L,1,MPI_REALTYPE,0,MPI_COMM_WORLD);
}

void evaluatedLdx(struct_plmd *plmd)
{
  int i,j,ij;
  int mID,oNID;

  oNID=omp_get_max_threads();
  MPI_Comm_rank(MPI_COMM_WORLD, &mID);

  regularize_gradient(plmd[0]);

  if (mID==0) {
  for (i=0; i<plmd->nblocks; i++) {
    plmd->dLdh[i]+=plmd->m1[oNID][i]/plmd->Z[oNID][0]-plmd->m1obs[i];
    for (j=0; j<plmd->nblocks; j++) {
      if (plmd->block2site[i] != plmd->block2site[j]) {
        ij=plmd->nblocks*i+j;
        plmd->dLdJ[ij]+=0.5L*(plmd->m2[oNID][ij]/plmd->Z[oNID][0]-plmd->m2obs[ij]);
      }
    }
  }
  }

  copyback_gradient(plmd);
}

void resetHinv(struct_plmd *plmd)
{
  int i;
  int mID;

  MPI_Comm_rank(MPI_COMM_WORLD, &mID);

  if (mID==0) {
  for (i=0; i<plmd->Jend*plmd->Jend; i++) {
    plmd->Hinv[i]=0.0L;
  }
  for (i=0; i<plmd->Jend; i++) {
    plmd->Hinv[(plmd->Jend+1)*i]=1.0L;
  }
  }
}

void updateHinv(struct_plmd *plmd)
{
  int i,j;
  real DxDg,DgHinvDg;
  real c1,c2;
  int mID;

  MPI_Comm_rank(MPI_COMM_WORLD, &mID);

  if (mID==0) {
  DxDg=0.0L;
  for (i=0;i<plmd->Jend;i++) {
    // Put Delta x and Delta dLdx in x0 and dLdx0, which hold previous values
    plmd->x0[i]=plmd->x[i]-plmd->x0[i];
    plmd->dLdx0[i]=plmd->dLdx[i]-plmd->dLdx0[i];
    DxDg+=plmd->x0[i]*plmd->dLdx0[i];
  }

  DgHinvDg=0.0L;
  for (i=0; i<plmd->Jend; i++) {
    plmd->hi[i]=0.0L;
    for (j=0; j<plmd->Jend; j++) {
      // put Hinv * Delta dLdx in hi (the search direction) as a buffer
      plmd->hi[i]+=plmd->Hinv[i*plmd->Jend+j]*plmd->dLdx0[j];
    }
    DgHinvDg+=plmd->hi[i]*plmd->dLdx0[i];
  }
  c1=(1.0L+DgHinvDg/DxDg)/DxDg;
  c2=-1.0L/DxDg;

  for (i=0; i<plmd->Jend; i++) {
    for (j=i; j<plmd->Jend; j++) {
      plmd->Hinv[i*plmd->Jend+j]+=c1*(plmd->x0[i]*plmd->x0[j])+c2*(plmd->hi[i]*plmd->x0[j]+plmd->x0[i]*plmd->hi[j]);
      plmd->Hinv[j*plmd->Jend+i]=plmd->Hinv[i*plmd->Jend+j];
    }
  }
  
  //   dx=xf-xi;
  //   dd=df-di;
  //   Hinv=Hinv+(1+(dd'*Hinv*dd)/(dx'*dd))*(dx*dx')/(dx'*dd)-((Hinv*dd*dx')+(Hinv*dd*dx')')/(dx'*dd);
  }
}

void projectHinv(struct_plmd *plmd)
{
  int i,j;
  real dLds;
  int mID;

  MPI_Comm_rank(MPI_COMM_WORLD, &mID);

  if (mID==0) {
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
  }
}

real lineL(real s,struct_plmd *plmd)
{
  int i;
  int mID;

  MPI_Comm_rank(MPI_COMM_WORLD, &mID);

  if (mID==0) {
  for (i=0; i<plmd->Jend; i++) {
    plmd->x[i]=plmd->x0[i]+s*plmd->hi[i];
  }
  }

  MPI_Bcast(plmd->x,plmd->Jend,MPI_REALTYPE,0,MPI_COMM_WORLD);

  evaluateL(plmd);
  return plmd->L[0];
}

int update(int step,struct_plmd *plmd)
{
  real smin,smid1,smid2,smax;
  real Lmin,Lmid1,Lmid2,Lmax;
  real L0;
  int i,ss;
  real phi=(1.0L+sqrtl(5.0L))/2.0L;
  int mID;

  MPI_Comm_rank(MPI_COMM_WORLD, &mID);

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
    MPI_Bcast(&smax,1,MPI_REALTYPE,0,MPI_COMM_WORLD);

    Lmax=lineL(smax,plmd);
    // fprintf(stderr,"Step %d smax %g L=%g -> L=%g\n",step,smax,L0,Lmax);
  }

  smid1=(smax-smin)/(phi*phi)+smin;
  Lmid1=lineL(smid1,plmd);

  smid2=-1.0L;
  
  for (ss=0; ss<20; ss++) {
    if (smid1<0.0L) {
      smid1=(smid2-smin)/phi+smin;
      MPI_Bcast(&smid1,1,MPI_REALTYPE,0,MPI_COMM_WORLD);
      Lmid1=lineL(smid1,plmd);
    } else {
      smid2=smax-(smax-smid1)/phi;
      MPI_Bcast(&smid2,1,MPI_REALTYPE,0,MPI_COMM_WORLD);
      Lmid2=lineL(smid2,plmd);
    }
    // fprintf(stderr,"Step %d smid1 %g smid2 %g L=%g -> L1=%g L2=%g\n",step,smid1,smid2,L0,Lmid1,Lmid2);
    if (mID==0) {
    fprintf(stderr,"Step %d smid1 %lg Lmin=%lg dL=%lg %lg %lg %lg\n",step,(double) smid1,(double) Lmin,(double) (Lmin-Lmin),(double) (Lmid1-Lmin),(double) (Lmid2-Lmin),(double) (Lmax-Lmin));
    }

    if (Lmid1<=Lmid2) {
      smax=smid2;
      MPI_Bcast(&smax,1,MPI_REALTYPE,0,MPI_COMM_WORLD);
      Lmax=Lmid2;
      smid2=smid1;
      MPI_Bcast(&smid2,1,MPI_REALTYPE,0,MPI_COMM_WORLD);
      Lmid2=Lmid1;
      smid1=-1.0L;
      MPI_Bcast(&smid1,1,MPI_REALTYPE,0,MPI_COMM_WORLD);
    } else {
      smin=smid1;
      MPI_Bcast(&smin,1,MPI_REALTYPE,0,MPI_COMM_WORLD);
      Lmin=Lmid1;
      smid1=smid2;
      MPI_Bcast(&smid1,1,MPI_REALTYPE,0,MPI_COMM_WORLD);
      Lmid1=Lmid2;
      smid2=-1.0L;
      MPI_Bcast(&smid2,1,MPI_REALTYPE,0,MPI_COMM_WORLD);
    }
  }
  if (mID==0) {
  fprintf(stderr,"Step %d smid1 %lg smid2 %lg L=%lg -> L1=%lg L2=%lg\n",step,(double) smid1,(double) smid2,(double) L0,(double) Lmid1,(double) Lmid2);
  }

  return (L0-plmd->L[0]<1e-15);
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

int itterate(int step,struct_plmd *plmd)
{
  evaluateL(plmd);
  evaluatedLdx(plmd);

  if (step==0) {
    resetHinv(plmd);
  } else {
    updateHinv(plmd);
  }

  projectHinv(plmd);
  return update(step,plmd);
}

void run(struct_plmd *plmd)
{
  int b_conv;
  int s;

  copyback_hJ(plmd);
  //initialize(plmd);
  b_conv=0;
  for (s=0; s<1000; s++) {
    b_conv=itterate(s,plmd);
    if (b_conv) {
      break;
    }
  }
}

void finish(struct_plmd *plmd,int argc, char *argv[])
{
  int i,j;
  FILE *fp;

  fp=fopen(argv[1],"w");
  if (fp==NULL) {
    fprintf(stderr,"Error: %s does not exist\n",argv[6]);
    exit(1);
  }
  for (i=0; i<plmd->nblocks; i++) {
    fprintf(fp," %lg",(double) plmd->h[i]);
  }
  fclose(fp);

  fp=fopen(argv[2],"w");
  if (fp==NULL) {
    fprintf(stderr,"Error: %s does not exist\n",argv[7]);
    exit(1);
  }
  for (i=0; i<plmd->nblocks; i++) {
    for (j=0; j<plmd->nblocks; j++) {
      fprintf(fp," %lg",(double) plmd->J[i*plmd->nblocks+j]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

  free(plmd->nsubs);
  free(plmd->block2site);
  // free(plmd->Seq);

  free(plmd->x);
  free(plmd->dLdx);
  free(plmd->Hinv);
  free(plmd->hi);
  free(plmd->x0);
  free(plmd->dLdx0);
}

int main(int argc, char *argv[])
{
  struct_plmd *plmd;

  MPI_Init(&argc,&argv);

  plmd = setup(argc,argv);
 
  run(plmd);

  finish(plmd,argc,argv);

  MPI_Finalize();

  return 0;
}
