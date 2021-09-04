
typedef struct struct_plmd {
  int nblocks;
  int nsites;
  int *nsubs;
  int *block0;
  int *block2site;
  double *h;
  double *J;
  double **Z;
  double **m1;
  double **m2;
} struct_plmd;

