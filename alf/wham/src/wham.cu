#include <stdio.h>
#include <math.h>

#include "wham.h"

// #define kB 0.008314
#define kB 0.00198614
#define MAXLENGTH 4096
#define BLOCK 256
#define SBLOCK 16

__device__ __host__ inline
double logadd(double lnA,double lnB)
{
  if (lnA>lnB) {
    return lnA+log(1+exp(lnB-lnA));
  } else if (lnA==(-INFINITY)) {
    return lnB;
  } else {
    return lnB+log(1+exp(lnA-lnB));
  }
}

__device__ inline
void atomic_logadd(double *p_lnA,double lnB)
{
  double lnA, lnC;
  double tmp_lnA;
  tmp_lnA=p_lnA[0];
  do {
    lnA=tmp_lnA;
    if (lnA>lnB) {
      lnC=lnA+log(1+exp(lnB-lnA));
    } else if (lnA==(-INFINITY)) {
      lnC=lnB;
    } else {
      lnC=lnB+log(1+exp(lnA-lnB));
    }
    // tmp_lnA=atomicCAS(p_lnA,lnA,lnC);
    // tmp_lnA=__int_as_float(atomicCAS((int*) p_lnA,__float_as_int(lnA),__float_as_int(lnC)));
    tmp_lnA=__longlong_as_double(atomicCAS((unsigned long long int*) p_lnA,__double_as_longlong(lnA),__double_as_longlong(lnC)));
  } while (lnA != tmp_lnA);
}

struct_data* readdata(int arg1, double arg2, int arg3, int arg4)
{
  const char *Edir,*Ddir;
  FILE *fp,*fpE,*fpQ;
  int i;
  int s1,s2;
  int j,jN;
  int iN;
  // int iB0,iB1;
  int B_N;
  int iB0;
  char fnm[MAXLENGTH];
  char line[MAXLENGTH];
  char *linebuffer;
  int ibuffer;
  int n;
  double E,q;
  int B0_MAX=INT_MIN;
  int B0_MIN=INT_MAX;
  // int B1_MAX=INT_MIN;
  // int B1_MIN=INT_MAX;
  struct_data *data;

  data=(struct_data*) malloc(sizeof(struct_data));

  data->Nsim=arg1;

  fp=fopen("../prep/nsubs","r");
  if (fp==NULL) {
    fprintf(stderr,"Error, ../prep/nsubs does not exist\n");
    exit(1);
  }
  data->Nsites=0;
  while (fscanf(fp,"%d",&i)==1) {
    data->Nsites++;
  }
  fclose(fp);
  data->Nsubs=(int*) malloc(data->Nsites*sizeof(int));
  data->block0=(int*) malloc((data->Nsites+1)*sizeof(int));
  fp=fopen("../prep/nsubs","r");
  data->Nblocks=0;
  data->block0[0]=0;
  for (i=0; i<data->Nsites; i++) {
    fscanf(fp,"%d",&(data->Nsubs[i]));
    data->Nblocks+=data->Nsubs[i];
    data->block0[i+1]=data->block0[i]+data->Nsubs[i];
  }
  fclose(fp);

  /*if (argc>3) {
    Edir=argv[2];
    Ddir=argv[3];
  } else {
    fprintf(stderr,"Error, missing input and output directory.\n");
    exit(1);
  }*/
  Edir="Energy";
  Ddir="Lambda";

  /*if (argc>4) {
    sscanf(argv[4],"%d",&(data->ms));
  } else {
    fprintf(stderr,"Error, 4th argument should indicate whether to use multisite parameters.\n");
    exit(1);
  }*/
  data->ms=arg3;

  /*if (argc>5) {
    sscanf(argv[5],"%d",&(data->msprof));
  } else {
    fprintf(stderr,"Error, 5th argument should indicate whether to use multisite profiles.\n");
    exit(1);
  }*/
  data->msprof=arg4;

  data->NL=data->Nblocks;
  data->NF=data->Nsim;

  data->Ndim=data->Nsim+data->NL+1+2; // E(1) Lambda(6) Energies(Nsim) BinInd(1) ReactCoord(1)

  data->T_h=(double*) malloc(data->NF*sizeof(double));
  data->beta_h=(double*) malloc(data->NF*sizeof(double));
  cudaMalloc(&(data->beta_d),data->NF*sizeof(double));

  for (i=0; i<data->NF; i++) {
    data->T_h[i]=arg2;
    data->beta_h[i]=1.0/(kB*data->T_h[i]);
  }
  cudaMemcpy(data->beta_d,data->beta_h,data->NF*sizeof(double),cudaMemcpyHostToDevice);

  data->beta_t=1.0/(kB*arg2);

  data->B[0].dx=0.1;
  data->B[1].dx=0.002500025;

  data->B2d[0].dx=0.1;
  data->B2d[1].dx=0.0500005;
  data->B2d[2].dx=0.0500005;

  data->n_h=(int*) malloc(data->NF*sizeof(int));
  cudaMalloc(&(data->n_d),data->NF*sizeof(int));
  data->ND=0;
  data->NDmax=0;
  data->NDmax+=MAXLENGTH;
  data->D_h=(double*) malloc(data->NDmax*data->Ndim*(sizeof(double)));
  data->i_h=(int*) malloc(data->NDmax*(sizeof(int)));
  data->lnw_h=(double*) malloc(data->NF*(sizeof(double)));

  for (i=0; i<data->NF; i++) {
    sprintf(fnm,"%s/ESim%d.dat",Edir,i+1);
    fpE=fopen(fnm,"r");
    if (fpE==NULL) {
      fprintf(stderr,"Error, energy file %s does not exist\n",fnm);
      exit(1);
    }
    sprintf(fnm,"%s/Lambda%d.dat",Ddir,i+1);
    fpQ=fopen(fnm,"r");
    if (fpQ==NULL) {
      fprintf(stderr,"Error, contact file %s does not exist\n",fnm);
      exit(1);
    }

    data->lnw_h[i]=(data->NF-i-1)*log(1.0);
    // data->lnw_h[i]=(data->NF-i-1)*(-1.0);

    n=0;

    while (fgets(line,MAXLENGTH,fpE) != NULL) {
      if (data->ND>=data->NDmax) {
        data->NDmax+=MAXLENGTH;
        data->D_h=(double*) realloc(data->D_h,data->NDmax*data->Ndim*sizeof(double));
        data->i_h=(int*) realloc(data->i_h,data->NDmax*sizeof(int));
      }
      n++;

      linebuffer=line;
      sscanf(linebuffer,"%lf%n",&E,&ibuffer);
      linebuffer+=ibuffer;
      data->D_h[data->ND*data->Ndim]=E;
      for (j=0; j<data->NF; j++) {
        sscanf(linebuffer,"%lf%n",&E,&ibuffer);
        linebuffer+=ibuffer;
        data->D_h[data->ND*data->Ndim+data->NL+1+j]=E;
      }
      // data->D_h[data->ND*data->Ndim]=E;

      fgets(line,MAXLENGTH,fpQ);
      linebuffer=line;
      for (j=0; j<data->Nblocks; j++) {
        sscanf(linebuffer,"%lf%n",&q,&ibuffer);
        linebuffer+=ibuffer;
        data->D_h[data->ND*data->Ndim+1+j]=q;
      }

      data->i_h[data->ND]=i;

      iB0=(int) floor(E/data->B[0].dx);
      // iB1=(int) floor(q/data->B[1].dx);
      if (iB0<B0_MIN) {
        B0_MIN=iB0;
      }
      if (iB0>B0_MAX) {
        B0_MAX=iB0;
      }
      // if (iB1<B1_MIN) {
      //   B1_MIN=iB1;
      // }
      // if (iB1>B1_MAX) {
      //   B1_MAX=iB1;
      // }
      
      data->ND++;
    }

    data->n_h[i]=n;

    fclose(fpE);
    fclose(fpQ);
  }

  data->B[0].min=B0_MIN*data->B[0].dx;
  data->B[0].max=(B0_MAX+1)*data->B[0].dx;
  data->B[0].N=(B0_MAX-B0_MIN)+1;

  data->B2d[0].min=data->B[0].min;
  data->B2d[0].max=data->B[0].max;
  data->B2d[0].N=data->B[0].N;

  data->B[1].min=0;
  data->B[1].max=1;
  data->B[1].N=400;
  B_N=data->B[1].N;

  data->B2d[1].min=0;
  data->B2d[1].max=1;
  data->B2d[1].N=20;
  data->B2d[2].min=0;
  data->B2d[2].max=1;
  data->B2d[2].N=20;
  if (data->B2d[1].N*data->B2d[2].N>B_N) {
    B_N=data->B2d[1].N*data->B2d[2].N;
  }

  cudaMemcpy(data->n_d,data->n_h,data->NF*sizeof(int),cudaMemcpyHostToDevice);
  cudaMalloc(&(data->D_d),data->NDmax*data->Ndim*sizeof(double));
  cudaMemcpy(data->D_d,data->D_h,data->NDmax*data->Ndim*sizeof(double),cudaMemcpyHostToDevice);

  cudaMalloc(&(data->i_d),data->NDmax*sizeof(int));
  cudaMemcpy(data->i_d,data->i_h,data->NDmax*sizeof(int),cudaMemcpyHostToDevice);

  cudaMalloc(&(data->lnw_d),data->NF*sizeof(double));
  cudaMemcpy(data->lnw_d,data->lnw_h,data->NF*sizeof(double),cudaMemcpyHostToDevice);

  data->lnDenom_h=(double*) malloc(data->NDmax*sizeof(double));
  cudaMalloc(&(data->lnDenom_d),data->NDmax*sizeof(double));

  /*
  V=[None]*2

  for iB in range(0,2):
    V[iB]=[None]*B[iB].N
    for i in range(0,B[iB].N):
      V[iB][i]=(B[iB].max-B[iB].min)/B[iB].N*i+B[iB].min
  */

  data->f_h=(double*) malloc(data->NF*sizeof(double));
  cudaMalloc(&(data->f_d),data->NF*sizeof(double));
  for (i=0; i<data->NF; i++) {
    data->f_h[i]=0.0;
  }
  cudaMemcpy(data->f_d,data->f_h,data->NF*sizeof(double),cudaMemcpyHostToDevice);
  data->invf_h=(double*) malloc(data->NF*sizeof(double));
  cudaMalloc(&(data->invf_d),data->NF*sizeof(double));

  fprintf(stderr,"Warning, DOS allocation is not sparse, requesting %d doubles\n",data->B[0].N*B_N);
  data->lnZ_h=(double*) malloc(B_N*sizeof(double));
  cudaMalloc(&(data->lnZ_d),B_N*sizeof(double));

  iN=0;
  for (s1=0; s1<data->Nsites; s1++) {
    for (s2=s1; s2<data->Nsites; s2++) {
      if (s1==s2) {
        if (data->Nsubs[s1]==2) {
          iN+=data->Nsubs[s1]+data->Nsubs[s1]*(data->Nsubs[s1]-1)/2;
        } else {
          iN+=data->Nsubs[s1]+2*data->Nsubs[s1]*(data->Nsubs[s1]-1)/2;
        }
      } else if (data->msprof) {
        iN+=data->Nsubs[s1]*data->Nsubs[s2];
      }
    }
  }
  data->iN=iN;

  jN=0;
  // data->jNij=(int*) malloc((data->Nsites*(data->Nsites+1))/2+1,sizeof(int));
  // i=0;
  for (s1=0; s1<data->Nsites; s1++) {
    for (s2=s1; s2<data->Nsites; s2++) {
      // data->jNij[i]=jN;
      if (s1==s2) {
        jN+=data->Nsubs[s1]+5*data->Nsubs[s1]*(data->Nsubs[s1]-1)/2;
      } else if (data->ms==1) {
        jN+=5*data->Nsubs[s1]*data->Nsubs[s2];
      } else if (data->ms==2) {
        jN+=data->Nsubs[s1]*data->Nsubs[s2];
      }
      // i++;
    }
  }
  // data->jNij[i]=jN;
  data->jN=jN;

  data->dlnZ_hN=(double**) malloc(jN*sizeof(double*));
  for (j=0; j<jN; j++) {
    data->dlnZ_hN[j]=(double*) malloc(B_N*sizeof(double));
  }
  cudaMalloc(&(data->dlnZ_d),B_N*sizeof(double));
  cudaMalloc(&(data->dlnZ_dN),jN*B_N*sizeof(double));
  data->Gimp_h=(double*) malloc(B_N*sizeof(double));
  cudaMalloc(&(data->Gimp_d),B_N*sizeof(double));

  data->C_h=(double*) malloc(jN*sizeof(double));
  cudaMalloc(&(data->C_d),jN*sizeof(double));
  data->CV_h=(double*) malloc(jN*sizeof(double));
  cudaMalloc(&(data->CV_d),jN*sizeof(double));
  data->CC_h=(double*) malloc(jN*jN*sizeof(double));
  cudaMalloc(&(data->CC_d),jN*jN*sizeof(double));
  
  // data->one_h=(double*) malloc(B_N*sizeof(double));
  // data->E_h=(double*) malloc(B_N*sizeof(double));
  // data->E2_h=(double*) malloc(B_N*sizeof(double));

  return data;
}

__global__ void resetlogdata(double *d,int N)
{
  int i=blockIdx.x*blockDim.x+threadIdx.x;
  if (i<N) {
    d[i]=-INFINITY;
  }
}

__global__ void resetdata(double *d,int N)
{
  int i=blockIdx.x*blockDim.x+threadIdx.x;
  if (i<N) {
    d[i]=0.0;
  }
}

__global__ void sumdenom(struct_data data)
{
  int t=blockIdx.x*blockDim.x+threadIdx.x;
  int i;
  double lnwn;
  double E;
  double lnDenom=-INFINITY;

  if (t<data.ND) {
    for (i=0; i<data.NF; i++) {
      E=data.D_d[t*data.Ndim+data.NL+1+i];
      lnwn=data.lnw_d[i]+log((double) data.n_d[i]);
      lnDenom=logadd(lnDenom,lnwn+data.f_d[i]-data.beta_d[i]*E);
    }
    data.lnDenom_d[t]=lnDenom;
  }
}

__global__
void getf(struct_data data)
{
  int t,tmin,tmax;
  int i;
  double lnw;
  double beta;
  double E;
  __shared__ double invf[BLOCK];

  tmin=(data.ND*threadIdx.x)/blockDim.x;
  tmax=(data.ND*(threadIdx.x+1))/blockDim.x;

  i=blockIdx.x;

  beta=data.beta_d[i];

  invf[threadIdx.x]=-INFINITY;
  for (t=tmin; t<tmax; t++) {
    // E=data.D_d[t][0];
    lnw=data.lnw_d[data.i_d[t]];
    E=data.D_d[t*data.Ndim+data.NL+1+i];
    invf[threadIdx.x]=logadd(invf[threadIdx.x],lnw-beta*E-data.lnDenom_d[t]);
  }

  __syncthreads();

  for (t=1; t<blockDim.x; t*=2) {
    if ((threadIdx.x % (2*t)) == 0) {
      invf[threadIdx.x]=logadd(invf[threadIdx.x],invf[threadIdx.x+t]);
    }
    __syncthreads();
  }

  if (threadIdx.x==0) {
    data.invf_d[i]=invf[0];
    data.f_d[i]=-invf[0];
  }
}

__global__ void normf(double *f,double f_avg,int N)
{
  int i=blockIdx.x*blockDim.x+threadIdx.x;
  if (i<N) {
    f[i]-=f_avg;
  }
}

void iteratedata(struct_data *data)
{
  int escape_flag=0;
  int itt;
  int mitt=1000;
  int i;
  double f_sum;
  FILE *fp;

  for (itt=0; itt<mitt; itt++) {
    escape_flag=1;
    // resetlogdata <<< (data->ND+BLOCK-1)/BLOCK, BLOCK >>> (data->lnDenom_d,data->ND);
    sumdenom <<< (data->ND+BLOCK-1)/BLOCK, BLOCK >>> (data[0]);
    // resetlogdata <<< (data->NF+BLOCK-1)/BLOCK, BLOCK >>> (data->invf_d,data->NF);
    getf <<< data->NF, BLOCK >>> (data[0]);

    cudaMemcpy(data->invf_h,data->invf_d,data->NF*sizeof(double),cudaMemcpyDeviceToHost);
    f_sum=0;
    for (i=0; i<data->NF; i++) {
      fprintf(stdout," %f",data->f_h[i]);
      // fprintf(stderr,"%d %d %f %f\n",itt,i,data->f_h[i],-data->invf_h[i]);
      if (data->f_h[i]+data->invf_h[i]<-0.00005 || data->f_h[i]+data->invf_h[i]>0.00005){
        escape_flag=0;
      }
      data->f_h[i]=-data->invf_h[i];
      f_sum+=data->f_h[i];
    }
    fprintf(stdout,"\n");
    f_sum/=data->NF;
    normf <<< (data->NF+BLOCK-1)/BLOCK,BLOCK >>> (data->f_d,f_sum,data->NF);
    if (escape_flag==1) {
      break;
    }
  }

  fp=fopen("f.dat","w");
  for (i=0; i<data->NF; i++) {
    fprintf(fp," %12.5f",data->f_h[i]);
  }
  fclose(fp);
}

__global__ void bin1(struct_data data,int i1)
{
  double q1;
  int t;
  int iB1;

  t=blockIdx.x*blockDim.x+threadIdx.x;

  if (t<data.ND) {
    q1=data.D_d[t*data.Ndim+1+i1];
    iB1=(int) floor((q1-data.B[1].min)/data.B[1].dx);
    data.D_d[t*data.Ndim+1+data.NL+data.Nsim]=iB1;
  }
}

__global__ void bin12(struct_data data,int i1,int i2)
{
  double q1,q2,q;
  int t;
  int iB12;

  t=blockIdx.x*blockDim.x+threadIdx.x;

  if (t<data.ND) {
    q1=data.D_d[t*data.Ndim+1+i1];
    q2=data.D_d[t*data.Ndim+1+i2];
    if (q1+q2>0.8) {
      q=q1/(q1+q2);
      iB12=(int) floor((q-data.B[1].min)/data.B[1].dx);
    } else {
      iB12=-1;
    }
    data.D_d[t*data.Ndim+1+data.NL+data.Nsim]=iB12;
  }
}

__global__ void bin2(struct_data data,int i1,int i2)
{
  double q1,q2;
  int t;
  int iB1,iB2;

  t=blockIdx.x*blockDim.x+threadIdx.x;

  if (t<data.ND) {
    q1=data.D_d[t*data.Ndim+1+i1];
    q2=data.D_d[t*data.Ndim+1+i2];
    iB1=(int) floor((q1-data.B2d[1].min)/data.B2d[1].dx);
    iB2=(int) floor((q2-data.B2d[2].min)/data.B2d[2].dx);
    data.D_d[t*data.Ndim+1+data.NL+data.Nsim]=iB1*data.B2d[1].N+iB2;
  }
}

double bin_all(struct_data *data,int *ptype,int i)
{
  double wnorm;
  int s1,s2;
  int i1,i2;
  int B_N;
  char fnm[MAXLENGTH];
  FILE *fp;

  for (s1=0; s1<data->Nsites; s1++) {
    for (s2=s1; s2<data->Nsites; s2++) {

      if (s1==s2) {

        for (i1=data->block0[s1]; i1<data->block0[s1+1]; i1++) {
          if (i==0) {
            fprintf(stderr,"1D Profile %d\n",i1);
            wnorm=1.0;
            bin1 <<< (data->ND+BLOCK-1)/BLOCK, BLOCK >>> (data[0],i1);
            sprintf(fnm,"G_imp/G1_%d.dat",data->Nsubs[s1]);
            B_N=data->B[1].N;
            ptype[0]=0;
          }
          i--;
        }

        for (i1=data->block0[s1]; i1<data->block0[s1+1]; i1++) {
          for (i2=i1+1; i2<data->block0[s1+1]; i2++) {
            if (i==0) {
              fprintf(stderr,"1D Profile %d,%d\n",i1,i2);
              wnorm=1.0/((data->Nsubs[s1]-1)/2.0);
              bin12 <<< (data->ND+BLOCK-1)/BLOCK, BLOCK >>> (data[0],i1,i2);
              sprintf(fnm,"G_imp/G12_%d.dat",data->Nsubs[s1]);
              B_N=data->B2d[1].N*data->B2d[2].N;
              ptype[0]=1;
            }
            i--;
          }
        }

        if (data->Nsubs[s1]>2) {
          for (i1=data->block0[s1]; i1<data->block0[s1+1]; i1++) {
            for (i2=i1+1; i2<data->block0[s1+1]; i2++) {
              if (i==0) {
                fprintf(stderr,"2D Profile %d,%d\n",i1,i2);
                wnorm=1.0/((data->Nsubs[s1]-1)/2.0);
                bin2 <<< (data->ND+BLOCK-1)/BLOCK, BLOCK >>> (data[0],i1,i2);
                sprintf(fnm,"G_imp/G2_%d.dat",data->Nsubs[s1]);
                B_N=data->B2d[1].N*data->B2d[2].N;
                ptype[0]=2;
              }
              i--;
            }
          }
        }

      } else if (data->msprof) { // Site-site interaction

        for (i1=data->block0[s1]; i1<data->block0[s1+1]; i1++) {
          for (i2=data->block0[s2]; i2<data->block0[s2+1]; i2++) {
            if (i==0) {
              fprintf(stderr,"2D SS Profile %d,%d\n",i1,i2);
              wnorm=1.0/(data->Nsubs[s1]*data->Nsubs[s2]);
              bin2 <<< (data->ND+BLOCK-1)/BLOCK, BLOCK >>> (data[0],i1,i2);
              sprintf(fnm,"G_imp/G1_%d_%d.dat",data->Nsubs[s1],data->Nsubs[s2]);
              B_N=data->B2d[1].N*data->B2d[2].N;
              ptype[0]=3;
            }
            i--;
          }
        }

      }
    }
  }

  fp=fopen(fnm,"r");
  if (fp==NULL) {
    fprintf(stderr,"Error, %s does not exist\n",fnm);
    exit(1);
  }
  for (i=0; i<B_N; i++) {
    fscanf(fp,"%lf",&data->Gimp_h[i]);
  }
  fclose(fp);

  cudaMemcpy(data->Gimp_d,data->Gimp_h,B_N*sizeof(double),cudaMemcpyHostToDevice);

  return wnorm;
}

// slope
__global__ void reactioncoord_phi(struct_data data,int i1)
{
  double q1;
  int t;

  t=blockIdx.x*blockDim.x+threadIdx.x;

  if (t<data.ND) {
    q1=data.D_d[t*data.Ndim+1+i1];
    data.D_d[t*data.Ndim+1+data.NL+data.Nsim+1]=q1;
  }
}

// quadratic
__global__ void reactioncoord_psi(struct_data data,int i1,int i2)
{
  double q1,q2;
  int t;

  t=blockIdx.x*blockDim.x+threadIdx.x;

  if (t<data.ND) {
    q1=data.D_d[t*data.Ndim+1+i1];
    q2=data.D_d[t*data.Ndim+1+i2];
    data.D_d[t*data.Ndim+1+data.NL+data.Nsim+1]=q1*q2;
  }
}

// omega - sharp endpoint
__global__ void reactioncoord_omega(struct_data data,int i1,int i2)
{
  double q1,q2;
  int t;

  t=blockIdx.x*blockDim.x+threadIdx.x;

  if (t<data.ND) {
    q1=data.D_d[t*data.Ndim+1+i1];
    q2=data.D_d[t*data.Ndim+1+i2];
    data.D_d[t*data.Ndim+1+data.NL+data.Nsim+1]=q2*(1-1/(q1/0.017+1));
  }
}

// slope
__global__ void reactioncoord_chi(struct_data data,int i1,int i2)
{
  double q1,q2;
  int t;

  t=blockIdx.x*blockDim.x+threadIdx.x;

  if (t<data.ND) {
    q1=data.D_d[t*data.Ndim+1+i1];
    q2=data.D_d[t*data.Ndim+1+i2];
    data.D_d[t*data.Ndim+1+data.NL+data.Nsim+1]=q2*(1-exp(-q1/0.18));
  }
}

void reactioncoord_all(struct_data *data,int i)
{
  int s1,s2;
  int j1,j2;

  for (s1=0; s1<data->Nsites; s1++) {
    for (s2=s1; s2<data->Nsites; s2++) {

      if (s1==s2) {

        for (j1=data->block0[s1]; j1<data->block0[s1+1]; j1++) {
          if (i==0) {
            reactioncoord_phi <<< (data->ND+BLOCK-1)/BLOCK, BLOCK >>> (data[0],j1);
          }
          i--;
        }

        for (j1=data->block0[s1]; j1<data->block0[s1+1]; j1++) {
          for (j2=j1+1; j2<data->block0[s1+1]; j2++) {
            if (i==0) {
              reactioncoord_psi <<< (data->ND+BLOCK-1)/BLOCK, BLOCK >>> (data[0],j1,j2);
            }
            i--;
          }
        }

        for (j1=data->block0[s1]; j1<data->block0[s1+1]; j1++) {
          for (j2=data->block0[s1]; j2<data->block0[s1+1]; j2++) {
            if (j1 != j2) {
              if (i==0) {
                reactioncoord_chi <<< (data->ND+BLOCK-1)/BLOCK, BLOCK >>> (data[0],j1,j2);
              }
              i--;
            }
          }
        }

        for (j1=data->block0[s1]; j1<data->block0[s1+1]; j1++) {
          for (j2=data->block0[s1]; j2<data->block0[s1+1]; j2++) {
            if (j1 != j2) {
              if (i==0) {
                reactioncoord_omega <<< (data->ND+BLOCK-1)/BLOCK, BLOCK >>> (data[0],j1,j2);
              }
              i--;
            }
          }
        }

      } else if (data->ms) { // Different sites

        for (j1=data->block0[s1]; j1<data->block0[s1+1]; j1++) {
          for (j2=data->block0[s2]; j2<data->block0[s2+1]; j2++) {
            if (i==0) {
              reactioncoord_psi <<< (data->ND+BLOCK-1)/BLOCK, BLOCK >>> (data[0],j1,j2);
            }
            i--;
          }
        }
        if (data->ms==1) {
        for (j1=data->block0[s1]; j1<data->block0[s1+1]; j1++) {
          for (j2=data->block0[s2]; j2<data->block0[s2+1]; j2++) {
            if (i==0) {
              reactioncoord_chi <<< (data->ND+BLOCK-1)/BLOCK, BLOCK >>> (data[0],j1,j2);
            }
            i--;
            if (i==0) {
              reactioncoord_chi <<< (data->ND+BLOCK-1)/BLOCK, BLOCK >>> (data[0],j2,j1);
            }
            i--;
          }
        }

        for (j1=data->block0[s1]; j1<data->block0[s1+1]; j1++) {
          for (j2=data->block0[s2]; j2<data->block0[s2+1]; j2++) {
            if (i==0) {
              reactioncoord_omega <<< (data->ND+BLOCK-1)/BLOCK, BLOCK >>> (data[0],j1,j2);
            }
            i--;
            if (i==0) {
              reactioncoord_omega <<< (data->ND+BLOCK-1)/BLOCK, BLOCK >>> (data[0],j2,j1);
            }
            i--;
          }
        }
        }
      }
    }
  }
}

__global__ void get_lnZ(struct_data data,double beta)
{
  double E;
  double lnw;
  int t,i;
  double *p_lnZ;
  __shared__ double loc_lnZ[400];
  int iB;

  t=blockIdx.x*blockDim.x+threadIdx.x;

  for (i=threadIdx.x; i<400; i+=blockDim.x) {
    loc_lnZ[i]=-INFINITY;
  }

  __syncthreads();

  for (i=SBLOCK*t; i<SBLOCK*(t+1); i++) {
    if (i<data.ND) {
      E=data.D_d[i*data.Ndim+0];
      iB=(int) data.D_d[i*data.Ndim+1+data.NL+data.Nsim];
      if (iB>=0) {
        lnw=data.lnw_d[data.i_d[i]];
        p_lnZ=&loc_lnZ[iB];
        atomic_logadd(p_lnZ,lnw-data.lnDenom_d[i]-beta*E);
      }
    }
  }

  __syncthreads();

  for (i=threadIdx.x; i<400; i+=blockDim.x) {
    if (loc_lnZ[i]>-INFINITY) {
      p_lnZ=&data.lnZ_d[i];
      atomic_logadd(p_lnZ,loc_lnZ[i]);
    }
  }
}

__global__ void get_dlnZ(struct_data data,int j1,double beta)
{
  double E,q;
  double lnw;
  int t,i;
  double *p_dlnZ;
  __shared__ double loc_dlnZ[400];
  int iB;

  t=blockIdx.x*blockDim.x+threadIdx.x;

  for (i=threadIdx.x; i<400; i+=blockDim.x) {
    loc_dlnZ[i]=-INFINITY;
  }

  __syncthreads();

  for (i=SBLOCK*t; i<SBLOCK*(t+1); i++) {
    if (i<data.ND) {
      E=data.D_d[i*data.Ndim+0];
      iB=(int) data.D_d[i*data.Ndim+1+data.NL+data.Nsim];
      q=data.D_d[i*data.Ndim+1+data.NL+data.Nsim+1];
      if (iB>=0) {
        lnw=data.lnw_d[data.i_d[i]];
        p_dlnZ=&loc_dlnZ[iB];
        atomic_logadd(p_dlnZ,lnw-data.lnDenom_d[i]-beta*E+log(q));
      }
    }
  }

  __syncthreads();

  for (i=threadIdx.x; i<400; i+=blockDim.x) {
    if (loc_dlnZ[i]>-INFINITY) {
      p_dlnZ=&data.dlnZ_dN[400*j1+i];
      atomic_logadd(p_dlnZ,loc_dlnZ[i]);
    }
  }
}

__global__ void get_CC(struct_data data,int i,double beta,double wnorm,int ptype)
{
  int j1=blockIdx.x;
  int j2=blockIdx.y;
  int k;
  double weight;
  double lnZ, dlnZ1, dlnZ2;
  double myCC;
  __shared__ double loc_CC[100];

  myCC=0;

  for (k=threadIdx.x;k<400;k+=100) {
    weight=wnorm;
    if ((ptype==0 || ptype==3) && k==400-1) {
      weight*=100.0;
    }
    lnZ=data.lnZ_d[k];
    if (isfinite(lnZ)) {
      dlnZ1=data.dlnZ_dN[j1*400+k];
      dlnZ2=data.dlnZ_dN[j2*400+k];
      myCC+=weight*exp(dlnZ1-lnZ)*exp(dlnZ2-lnZ);
    }
  }

  loc_CC[threadIdx.x]=myCC;

  __syncthreads();

  for (k=1; k<100; k*=2) {
    if (threadIdx.x%(2*k) == 0) {
      if (threadIdx.x+k < 100) {
        loc_CC[threadIdx.x]+=loc_CC[threadIdx.x+k];
      }
    }
    __syncthreads();
  }

  if (threadIdx.x==0) {
    // jN=data->NL+5*data->NL*(data->NL-1)/2;
    // jN=gridDim.x;
    data.CC_d[gridDim.x*j1+j2]=loc_CC[0];
  }
}

void getfofq(struct_data *data,double beta)
{
  int B_N;
  int i,iN;
  int j1,j2,jN;
  int k;
  double *C, *V;
  double wnorm;
  int ptype; // profile type, affects bin weight
  double weight;
  char fnm[MAXLENGTH];
  FILE *fpC,*fpV,*fp;

  B_N=data->B[1].N;
  if (data->B2d[1].N*data->B2d[2].N>B_N) {
    B_N=data->B2d[1].N*data->B2d[2].N;
  }

  iN=data->iN;
  jN=data->jN;

  C=(double*) malloc((jN+iN)*(jN+iN)*sizeof(double));
  V=(double*) malloc((jN+iN)*sizeof(double));
  for (j1=0;j1<(jN+iN);j1++) {
    for (j2=0;j2<(jN+iN);j2++) {
      C[j1*(jN+iN)+j2]=0;
    }
    V[j1]=0;
  }

  sumdenom <<< (data->ND+BLOCK-1)/BLOCK, BLOCK >>> (data[0]);

  for (i=0; i<iN; i++) {
    wnorm=bin_all(data,&ptype,i);

    resetlogdata <<< (B_N+BLOCK-1)/BLOCK, BLOCK >>> (data->lnZ_d,B_N);
    get_lnZ <<< (data->ND+(100*SBLOCK)-1)/(100*SBLOCK), 100 >>> (data[0],data->beta_t);
    cudaMemcpy(data->lnZ_h,data->lnZ_d,B_N*sizeof(double),cudaMemcpyDeviceToHost);

    sprintf(fnm,"multisite/G%d.dat",i+1);
    fp=fopen(fnm,"w");
    for (k=0; k<B_N; k++) {
      fprintf(fp,"%g\n",(-data->lnZ_h[k]-data->Gimp_h[k])/data->beta_t);
    }
    fclose(fp);

    for (k=0;k<B_N;k++) {
      weight=wnorm;
      if ((ptype==0 || ptype==3) && k==400-1) {
        weight*=100.0;
      }
      if (isfinite(data->lnZ_h[k])) {
        V[jN+i]+=weight*(-data->lnZ_h[k]-data->Gimp_h[k])/data->beta_t;
        C[(jN+i)*(jN+iN)+jN+i]+=weight;
      }
      // if (C[(jN+i)*(jN+iN)+jN+i]==0) {
      //   C[(jN+i)*(jN+iN)+jN+i]=1.0;
      // }
    }

    for (j1=0; j1<jN; j1++) {
      reactioncoord_all(data,j1);
      resetlogdata <<< (B_N+BLOCK-1)/BLOCK, BLOCK >>> (&(data->dlnZ_dN[B_N*j1]),B_N);
      get_dlnZ <<< (data->ND+(100*SBLOCK)-1)/(100*SBLOCK), 100 >>> (data[0],j1,data->beta_t);
      cudaMemcpy(data->dlnZ_hN[j1],&(data->dlnZ_dN[B_N*j1]),B_N*sizeof(double),cudaMemcpyDeviceToHost);

      /*
      sprintf(fnm,"dG%d_d%d.dat",i+1,j1+1);
      fp=fopen(fnm,"w");
      for (k=0; k<B_N; k++) {
        fprintf(fp,"%g\n",exp(data->dlnZ_hN[j1][k]-data->lnZ_h[k]));
      }
      fclose(fp);
      */
    }

    get_CC <<< make_uint3(jN,jN,1), 100 >>> (data[0],i,data->beta_t,wnorm,ptype);
    cudaMemcpy(data->CC_h,data->CC_d,jN*jN*sizeof(double),cudaMemcpyDeviceToHost);

    for (j1=0; j1<jN; j1++) {
      for (k=0;k<B_N;k++) {
        weight=wnorm;
        if ((ptype==0 || ptype==3) && k==400-1) {
          weight*=100.0;
        }
        if (isfinite(data->lnZ_h[k])) {
          if (isfinite(data->Gimp_h[k])==0) {
            fprintf(stderr,"Fatal error, implicit constraint entropy is undefined at bin %d\n",k);
            exit(1);
          }
          V[j1]+=weight*exp(data->dlnZ_hN[j1][k]-data->lnZ_h[k])*(-data->lnZ_h[k]-data->Gimp_h[k])/data->beta_t;
          C[j1*(jN+iN)+jN+i]+=weight*exp(data->dlnZ_hN[j1][k]-data->lnZ_h[k]);
          C[(jN+i)*(jN+iN)+j1]+=weight*exp(data->dlnZ_hN[j1][k]-data->lnZ_h[k]);
        }
      }

      for (j2=0; j2<jN; j2++) {
        /*
        for (k=0;k<B_N;k++) {
          if (wnorm==1.0 && k==B_N-1) {
            weight=100.0;
          } else {
            weight=wnorm;
          }
          if (isfinite(data->lnZ_h[k])) {
            C[j1*(jN+iN)+j2]+=weight*exp(data->dlnZ_hN[j1][k]-data->lnZ_h[k])*exp(data->dlnZ_hN[j2][k]-data->lnZ_h[k]);
          }
        }
        */
        C[j1*(jN+iN)+j2]+=data->CC_h[j1*jN+j2];
      }
    }
  }

  fpC=fopen("multisite/C.dat","w");
  fpV=fopen("multisite/V.dat","w");
  for (j1=0; j1<(jN+iN); j1++) {
    for (j2=0; j2<(jN+iN); j2++) {
      fprintf(fpC," %f",C[j1*(jN+iN)+j2]);
    }
    fprintf(fpC,"\n");
    fprintf(fpV," %f\n",V[j1]);
  }
  fclose(fpC);
  fclose(fpV);
}

// Changed main to wham to prevent compiler warnings, since it is just being used as a library, not an executable
extern "C" int wham(int arg1, double arg2, int arg3, int arg4)
{
  struct_data *data;
  
  data=readdata(arg1,arg2,arg3,arg4);

  iteratedata(data);

  getfofq(data,data->beta_t);

  return 0;
}

