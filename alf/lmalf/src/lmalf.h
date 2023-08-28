// typedef long double real;
// typedef double real;
typedef float real;

typedef struct struct_plmd {
  int nblocks;
  int nsites;
  int *nsubs;
  int *block0;
  int *block0_d;
  int *block2site;
  int B; // number of frames
  real kT;
  real *kx; // Regularization constants
  real *kx_d;
  real *kprofile;
  real *kprofile_d;
  int ms; // multisite coupling flag
  int msprof;
  real *lambda;
  real *lambda_d;
  real *mc_lambda;
  real *mc_lambda_d;
  real *ensweight;
  real *ensweight_d;
  real *mc_ensweight;
  real *mc_ensweight_d;
  real *mc_weight;
  int nbias;
  int nprof;
// Calculation variables
  real *weight_d;
  real *mc_weight_d;
  real *E_d;
  real *dEds_d;
  real *mc_E_d;
  real *mc_dEds_d;
  real *Z_d;
  real *mc_Z_d;
  real *Zprofile_d;
  real *mc_Zprofile_d;
  real *dLdZprofile_d;
  real *mc_dLdZprofile_d;
  real *dLdE_d;
  real *Gimp_d;
  real *G_d;
  real *Esum_d;
  real *dEdssum_d;
  real *mc_dEdssum_d;
  real *moments_d;
  real *mc_moments_d;
  real *sumensweight_d;
// Minimization variables
  int nx;
  real *L; // log likelihood
  real *L_d;
  real *x;
  real *x_d;
  real *x0;
  real *dLdx;
  real *dLdx_d;
  real *dLdx0;
  real *dLds;
  real *dLds_d;
  real *hi;
  real *dxds_d;
  int doneCount;
  bool done;
  real criteria;
  // real *Hinv;
  int Nmem;
  int Nmemax;
  real *d_x; // Nmemax * Jend
  real *d_dLdx; // Nmemax * Jend
  real *rho; // Nmemax
  real *alpha; // Nmemax
  real *beta; // Nmemax
  FILE *fplog;
} struct_plmd;

