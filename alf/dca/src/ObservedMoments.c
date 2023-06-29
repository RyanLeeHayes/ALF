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

#include "ObservedMoments.h"

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

  plmd->filter=(int*) calloc(plmd->nsites,sizeof(int));
  plmd->m0=0;
  plmd->m1=(int*) calloc(plmd->nblocks,sizeof(int));
  plmd->m2=(int*) calloc(plmd->nblocks*plmd->nblocks,sizeof(int));

  plmd->fp=fopen(argv[1],"r");
  if (plmd->fp==NULL) {
    fprintf(stderr,"Error: %s does not exist\n",argv[1]);
    exit(1);
  }

  plmd->fpout1=fopen(argv[2],"w");
  if (plmd->fpout1==NULL) {
    fprintf(stderr,"Error: %s does not exist\n",argv[2]);
    exit(1);
  }

  plmd->fpout2=fopen(argv[3],"w");
  if (plmd->fpout2==NULL) {
    fprintf(stderr,"Error: %s does not exist\n",argv[3]);
    exit(1);
  }

  return plmd;
}

int read_filter(struct_plmd *plmd)
{
  int i,c;

  for (i=0; i<plmd->nsites; i++) {
    c=fscanf(plmd->fp,"%d",&plmd->filter[i]);
  }
  return c;
}

void run(struct_plmd *plmd)
{
  int i,j;
  int isite,jsite;

  while (read_filter(plmd) == 1) {
    plmd->m0++;
    for (i=0; i<plmd->nblocks; i++) {
      isite=plmd->block2site[i];
      if (plmd->filter[isite]==i-plmd->block0[isite]) {
        plmd->m1[i]++;
        for (j=0; j<plmd->nblocks; j++) {
          jsite=plmd->block2site[j];
          if (plmd->filter[jsite]==j-plmd->block0[jsite]) {
            plmd->m2[i*plmd->nblocks+j]++;
          }
        }
      }
    }
  }
}

void finish(struct_plmd *plmd)
{
  int i,j;

  free(plmd->nsubs);
  free(plmd->block2site);
  // free(plmd->Seq);

  for (i=0; i<plmd->nblocks; i++) {
    fprintf(plmd->fpout1," %12.10f",(1.0*plmd->m1[i])/plmd->m0);
    for (j=0; j<plmd->nblocks; j++) {
      fprintf(plmd->fpout2," %12.10f",(1.0*plmd->m2[i*plmd->nblocks+j])/plmd->m0);
    }
    fprintf(plmd->fpout2,"\n");
  }
  fprintf(plmd->fpout1,"\n");

  free(plmd->filter);
  free(plmd->m1);
  free(plmd->m2);
}

int main(int argc, char *argv[])
{
  struct_plmd *plmd;

  MPI_Init(&argc,&argv);

  plmd = setup(argc,argv);
 
  run(plmd);

  finish(plmd);

  MPI_Finalize();

  return 0;
}
