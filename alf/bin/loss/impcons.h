// double is comparable on data center gpus and 4:3 slower on gaming gpus.
// float causes crashes on multisite systems and weird drift on single site
// systems with many samples. Teh stability of double is worth the mild extra
// cost.
// typedef long double real;
typedef double real;
// typedef float real;

typedef enum emode {
  gimp,
  sample,
  emodeend}
emode;

typedef enum eimpcons {
  fnex2011,
  fnexdozen2024,
  fnpwise2026,
  eimpconsend}
eimpcons;

typedef struct struct_plmd {
  int nblocks;
  int nsites;
  int *nsubs;
  int *block0;
  int *block0_d;
  int B; // number of frames
  real *mc_lambda;
  real *mc_lambda_d;
  real *mc_ensweight;
  real *mc_ensweight_d;
  int nprof;
  real *mc_Zprofile_d;
  emode mode;
  eimpcons impcons;
  real fnexp_c;
  real dozen_b;
  real pwise_w;
  real pwise_k;
  int NBINS2;
  int NBINS;
} struct_plmd;

