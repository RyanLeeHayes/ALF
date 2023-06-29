// Written by Ryan Hayes 2017-06-23
// Compute moments and generate a random sequence given h and J

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

#include "GenerateMoments.h"

#define MAXLENGTH 1024

struct_plmd* setup(int argc, char *argv[])
{
  struct_plmd *plmd;
  int i,j,k;
  FILE *fp;

  plmd=(struct_plmd*) malloc(sizeof(struct_plmd));

  fp=fopen("../../nblocks","r");
  if (fp==NULL) {
    fprintf(stderr,"Error: %s does not exist\n","../../nblocks");
    exit(1);
  }
  fscanf(fp,"%d",&(plmd->nblocks));
  fclose(fp);

  fp=fopen("../../nsubs","r");
  if (fp==NULL) {
    fprintf(stderr,"Error: %s does not exist\n","../../nsubs");
    exit(1);
  }
  i=0;
  plmd->nsites=0;
  while(i<plmd->nblocks) {
    fscanf(fp,"%d",&j);
    i+=j;
    plmd->nsites++;
  }
  fclose(fp);

  fp=fopen("../../nsubs","r");
  i=0;
  plmd->nsubs=(int*) calloc(plmd->nsites,sizeof(int));
  for(i=0; i<plmd->nsites; i++) {
    fscanf(fp,"%d",&(plmd->nsubs[i]));
    plmd->nsubs[i]++;
  }
  plmd->nblocks+=plmd->nsites;
  fclose(fp);

  plmd->block0=(int*) calloc(plmd->nsites,sizeof(int));
  plmd->block2site=(int*) calloc(plmd->nblocks,sizeof(int));
  k=0;
  for(i=0; i<plmd->nsites; i++) {
    plmd->block0[i]=k;
    for(j=0; j<plmd->nsubs[i]; j++) {
      plmd->block2site[k]=i;
      k++;
    }
  }

  plmd->h=(double*) calloc(plmd->nblocks,sizeof(double));
  plmd->J=(double*) calloc(plmd->nblocks*plmd->nblocks,sizeof(double));

  fp=fopen(argv[1],"r");
  if (fp==NULL) {
    fprintf(stderr,"Error: %s does not exist\n",argv[1]);
    exit(1);
  }
  for(i=0; i<plmd->nblocks; i++) {
    fscanf(fp,"%lg\n",&plmd->h[i]);
  }
  fclose(fp);

  fp=fopen(argv[2],"r");
  if (fp==NULL) {
    fprintf(stderr,"Error: %s does not exist\n",argv[2]);
    exit(1);
  }
  for(i=0; i<plmd->nblocks; i++) {
    for(j=0; j<plmd->nblocks; j++) {
      fscanf(fp,"%lg\n",&plmd->J[i*plmd->nblocks+j]);
    }
  }
  fclose(fp);

  plmd->Z=(double**) calloc(omp_get_max_threads()+1,sizeof(double*));
  plmd->m1=(double**) calloc(omp_get_max_threads()+1,sizeof(double*));
  plmd->m2=(double**) calloc(omp_get_max_threads()+1,sizeof(double*));
  for (i=0; i<omp_get_max_threads()+1; i++) {
    plmd->Z[i]=(double*) calloc(1,sizeof(double));
    plmd->m1[i]=(double*) calloc(plmd->nblocks,sizeof(double));
    plmd->m2[i]=(double*) calloc(plmd->nblocks*plmd->nblocks,sizeof(double));
  }

  return plmd;
}

void run(struct_plmd *plmd) {
  long long int s;
  long long int prod;
  int Seq[plmd->nsites];
  int i,j;
  double H,w;
  double H0,HC;
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
  fprintf(stderr,"%12ld to %12ld on %d\nStarting sequence %s on %d\n",smin,smax,ID,line,ID);

  // for (s=0; s<prod; s++)
  for (s=smin; s<smax; s++) {
    #pragma omp master
    {
    if (s%1000000==0 && mID==0) {
      fprintf(stderr,"%12ld sequences scanned by master of %12ld\n",s,smax);
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

void reduce(struct_plmd *plmd)
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

  MPI_Reduce(plmd->Z[0],plmd->Z[NID],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(plmd->m1[0],plmd->m1[NID],plmd->nblocks,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(plmd->m2[0],plmd->m2[NID],plmd->nblocks*plmd->nblocks,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
}

void finish(struct_plmd *plmd,int argc, char *argv[])
{
  int i,j;
  FILE *fp;
  int oNID,mID;

  oNID=omp_get_max_threads();

  MPI_Comm_rank(MPI_COMM_WORLD, &mID);

  if (mID==0) {

  fp=fopen(argv[3],"w");
  if (fp==NULL) {
    fprintf(stderr,"Error: %s does not exist\n",argv[3]);
    exit(1);
  }
  for (i=0; i<plmd->nblocks; i++) {
    fprintf(fp," %g",plmd->m1[oNID][i]/plmd->Z[oNID][0]);
  }
  fprintf(fp,"\n");
  fclose(fp);

  fp=fopen(argv[4],"w");
  if (fp==NULL) {
    fprintf(stderr,"Error: %s does not exist\n",argv[4]);
    exit(1);
  }
  for (i=0; i<plmd->nblocks; i++) {
    for (j=0; j<plmd->nblocks; j++) {
      fprintf(fp," %g",plmd->m2[oNID][i*plmd->nblocks+j]/plmd->Z[oNID][0]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

  }
}

int main(int argc, char *argv[])
{
  struct_plmd *plmd;

  MPI_Init(&argc,&argv);

  plmd = setup(argc,argv);

  #pragma omp parallel
  {
  run(plmd);
  }
  reduce(plmd);

  finish(plmd,argc,argv);

  MPI_Finalize();

  return 0;
}
