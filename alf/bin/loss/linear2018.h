typedef struct struct_hist {
  double max;
  double min;
  int N;
  double dx;
} struct_hist;

typedef struct struct_data {
  int NL;
  int NF;
  int ms; // whether to calculate multisite interaction parameters
  int msprof; // whether to calculate intersite profiles
  double *T_h;
  double *beta_h;
  double *beta_d;
  double beta_t;
  struct_hist B[2];
  struct_hist B2d[3];
  int ND;
  int NDmax;
  int Nsim;
  int Ndim;
  double *D_h; // data
  double *D_d;
  int *i_h; // simulation index
  int *i_d;
  double *lnw_h; // simulation weight
  double *lnw_d;
  double *lnDenom_h;
  double *lnDenom_d;
  int *n_h;
  int *n_d;
  double *f_h;
  double *f_d;
  double *invf_h;
  double *invf_d;
  double *lnZ_h;
  double *lnZ_d;
  double **dlnZ_hN;
  double *dlnZ_dN;
  double *dlnZ_d;
  double *Gimp_h;
  double *Gimp_d;
  double *CC_h;
  double *CC_d;
  double *C_h;
  double *C_d;
  double *CV_h;
  double *CV_d;
  int Nblocks;
  int Nsites;
  int *Nsubs;
  int *block0;
  int iN;
  int jN;
  // int *jNij;
} struct_data;
