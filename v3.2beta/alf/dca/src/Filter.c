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

#include "Filter.h"

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

  plmd->lambda=(real*) calloc(plmd->nblocks,sizeof(real));
  plmd->filter=(int*) calloc(plmd->nsites,sizeof(int));

  if (sscanf(argv[1],"%d",&plmd->freq) != 1) {
    fprintf(stderr,"Error: %s is not a valid frequency for reading frames\n",argv[1]);
    exit(1);
  }

  plmd->fp=fopen(argv[2],"r");
  if (plmd->fp==NULL) {
    fprintf(stderr,"Error: %s does not exist\n",argv[2]);
    exit(1);
  }

  plmd->fpout=fopen(argv[3],"w");
  if (plmd->fpout==NULL) {
    fprintf(stderr,"Error: %s does not exist\n",argv[3]);
    exit(1);
  }

  return plmd;
}

int read_lambda(struct_plmd *plmd)
{
  int i,c;

  for (i=0; i<plmd->nblocks; i++) {
    c=fscanf(plmd->fp,"%lg",&plmd->lambda[i]);
  }
  return c;
}

void write_filter(struct_plmd *plmd)
{
  int i;
  for (i=0; i<plmd->nsites; i++) {
    fprintf(plmd->fpout," %d",plmd->filter[i]);
  }
  fprintf(plmd->fpout,"\n");
}

void run(struct_plmd *plmd)
{
  int i,j;
  int count=0;

  while (read_lambda(plmd) == 1) {
    if ((count%plmd->freq)==0) {
      for (i=0; i<plmd->nsites; i++) {
        plmd->filter[i]=0;
        for (j=0; j<plmd->nsubs[i]; j++) {
          if (plmd->lambda[plmd->block0[i]+j]>0.99) {
            plmd->filter[i]=j+1;
          }
        }
      }
      write_filter(plmd);
    }
    count++;
  }
}

void finish(struct_plmd *plmd)
{
  free(plmd->nsubs);
  free(plmd->block2site);
  // free(plmd->Seq);

  free(plmd->lambda);
  free(plmd->filter);
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
