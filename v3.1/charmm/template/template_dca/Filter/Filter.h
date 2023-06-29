typedef double real;

typedef struct struct_plmd {
  int nblocks;
  int nsites;
  int *nsubs;
  int *block0;
  int *block2site;
  FILE *fp;
  FILE *fpout;
  real *lambda;
  int *filter;
  int freq;
} struct_plmd;
