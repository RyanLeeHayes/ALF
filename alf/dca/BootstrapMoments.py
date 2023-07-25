#! /usr/bin/env python

def BootstrapMomentsDCA(alf_info,NF,Path,NBS=50):
  """
  Bootstraps moments for Potts model analysis

  Used in alf.BSMomentDCA, see further documentation there. 

  Parameters
  ----------
  alf_info : dict
      Dictionary of variables alf needs to run
  NF : int
      Number of independent trials during production
  Path : str
      Path to a data directory for Potts model estimation, typically 'data'
  NBS : int, optional
      Number of bootstrap samples to take. (Default is 50)
  """

  import numpy as np
  import sys, os

  nsubs=alf_info['nsubs']+0 # add 0 to copy
  nblocks=np.sum(nsubs)
  nsites=np.size(nsubs,0)

  nsubs+=1
  nblocks+=nsites

  m1=np.zeros((1,nblocks,NF))
  m2=np.zeros((nblocks,nblocks,NF))

  for ifile in range(NF):
    m1[:,:,ifile]=np.loadtxt(Path+'/m1.'+str(ifile)+'.obs.dat')
    m2[:,:,ifile]=np.loadtxt(Path+'/m2.'+str(ifile)+'.obs.dat')

  np.savetxt(Path+'/m1.obs.dat',np.mean(m1,axis=2))
  np.savetxt(Path+'/m2.obs.dat',np.mean(m2,axis=2))

  np.random.seed(2401)
  for i in range(NBS):
    bs=np.random.randint(0,NF,(NF,))
    m1mean=np.mean(m1[:,:,bs],axis=2)
    m2mean=np.mean(m2[:,:,bs],axis=2)

    np.savetxt(Path+'/bs'+str(i)+'.dat',bs,fmt='%d')
    np.savetxt(Path+'/m1.bs'+str(i)+'.obs.dat',m1mean)
    np.savetxt(Path+'/m2.bs'+str(i)+'.obs.dat',m2mean)
