typedef double real;

typedef struct struct_plmd {
  int nblocks;
  int nsites;
  int *nsubs;
  int *block0;
  int *block2site;
  FILE *fp;
  FILE *fpout1;
  FILE *fpout2;
  int *filter;
  int m0;
  int *m1;
  int *m2;
} struct_plmd;
