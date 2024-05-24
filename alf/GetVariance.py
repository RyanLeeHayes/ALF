#! /usr/bin/env python

def GetVariance(alf_info,NF,NBS=50,lc=0.99):
  """
  Calculates the free energy change upon alchemical mutation

  Estimates free energy changes using the histogram based estimator and
  alchemical trajectories in analysis[i]/data/Lambda.[idupl].[irep].dat,
  where [i] is the cycle of alf, [idupl] is the numerical index for each
  of the indepdendent trials, and [irep] is the replica index for replica
  exchange. The histogram based estimator uses a hard-coded lambda cutoff
  lc of 0.99. The free energies are computed as -kT*log(P), where P is the
  number of times an alchemical state occurs in the trajectory. The energy
  of the biases at the endpoints is added back in. Uncertainty is
  estimated by bootstrapping from the independent trials [NBS] times.
  The default/recommended value for NBS is 50, which strikes a balance
  between computational expense and accuracy. Bootstrapped uncertainties
  tend to be slightly low becaause rerunning the entire alf workflow will
  result in different biases, which allow different levels of convergence,
  whereas all samples bootstrapped by GetVariance were run with the same
  bias. Results are written to analysis[i]/Result.txt. If the are multiple
  replicas, each replica will contribute to analysis[i]/Result.txt, but
  the result for each isolated replica is also written to
  analysis[i]/Result.[irep].txt.

  This routine is called by the routine alf.postprocess, and should be run
  from the analysis[i] directory.

  Parameters
  ----------
  alf_info : dict
      Dictionary of variables alf needs to run
  NF : int
      The number of independent trials run by runprod
  NBS : int, optional
      The number of bootstrap samples to take (default is 50)
  lc : float, optional
      The lambda cutoff for the histogram estimator (default is 0.99)
  """

  import sys, os, os.path
  import numpy as np
  import copy

  nblocks=alf_info['nblocks']
  nsubs=alf_info['nsubs']
  nreps=alf_info['nreps']
  nlig=np.prod(nsubs)

  b=np.loadtxt('b_prev.dat')
  if os.path.isfile('b_corr.dat'):
    b=b+np.loadtxt('b_corr.dat')
  c=np.loadtxt('c_prev.dat')
  x=np.loadtxt('x_prev.dat')
  s=np.loadtxt('s_prev.dat')

  b_shift=np.loadtxt('../nbshift/b_shift.dat')
  c_shift=np.loadtxt('../nbshift/c_shift.dat')
  x_shift=np.loadtxt('../nbshift/x_shift.dat')
  s_shift=np.loadtxt('../nbshift/s_shift.dat')
  ncentral=alf_info['ncentral']

  f=np.loadtxt('f.dat')

  kT=0.001987*alf_info['temp']

  G=np.zeros((NF,nlig))
  ind=np.zeros((nlig,len(nsubs)),dtype='int')
  for i in range(1,nlig):
    ind[i,:]=ind[i-1,:]
    for j in range(len(nsubs)-1,-1,-1):
      if (ind[i,j]+1)<nsubs[j]:
        ind[i,j]+=1
        break
      else:
        ind[i,j]=0
  blk=copy.deepcopy(ind)
  for i in range(1,len(nsubs)):
    blk[:,i:]+=nsubs[i-1]

  Eall=np.zeros((nreps,NF,nlig))
  Eshift=np.zeros((nreps,NF,nlig))
  lndenom=np.zeros((nreps,NF,nlig))
  nframes=np.zeros((nreps,NF))
  for irep in range(0,nreps):
    for i in range(0,NF):
      isim=i*nreps+irep
      L=np.loadtxt('data/Lambda.'+str(i)+'.'+str(irep)+'.dat')
      nframes[irep,i]=L.shape[0]
      for j in range(0,nlig):
        LList=np.zeros((nblocks,))
        LList[blk[j,:]]=1
        # E=np.sum(b[blk[j,:]])
        Eall[irep,i,j]=np.dot(LList,-b)+np.dot(np.dot(LList,-c),LList)+np.dot(np.dot(1-np.exp(-5.56*LList),-x),LList)+np.dot(np.dot(LList/(LList+0.017),-s),LList)
        Eshift[irep,i,j]=(irep-ncentral)*(np.dot(LList,-b_shift)+np.dot(np.dot(LList,-c_shift),LList)+np.dot(np.dot(1-np.exp(-5.56*LList),-x_shift),LList)+np.dot(np.dot(LList/(LList+0.017),-s_shift),LList))
        # lndenom[irep,i,j]=np.log(nframes[irep,i])+f[isim]-(Eall[irep,i,j]+Eshift[irep,i,j])/kT
        lndenom[irep,i,j]=np.log(nframes[irep,i])+f[isim]-(Eshift[irep,i,j])/kT

  PkeepA=np.zeros((nreps,NF,nlig))
  for irep in range(0,nreps):
    Pkeep=np.zeros((NF,nlig))
    G=np.zeros((NF,nlig))
    for i in range(0,NF):
      # print(i) # Was useful for keeping track of progress on slow analysis runs
      L=np.loadtxt('data/Lambda.'+str(i)+'.'+str(irep)+'.dat')
      for j in range(0,nlig):
        P=np.sum(np.all(L[:,blk[j,:]]>lc,axis=1))
        Pkeep[i,j]=P
        G[i,j]=-Eall[irep,i,j]-Eshift[irep,i,j]-kT*np.log(P)
    PkeepA[irep,:,:]=Pkeep

    np.savetxt('G.'+str(irep)+'.dat',G)

    Gmin=np.min(G,axis=0)
    Value=Gmin-kT*np.log(np.mean(np.exp(-(G-Gmin)/kT),axis=0))
    Value-=Value[0]

    np.random.seed(2401)
    GS=np.zeros((NBS,nlig))
    for i in range(0,NBS):
      GS[i,:]=Gmin-kT*np.log(np.mean(np.exp(-(G[np.random.randint(0,NF,(NF,1)),:]-Gmin)/kT),axis=0))
    Error=np.std(GS,axis=0)

    fp=open('Result.'+str(irep)+'.txt','w')

    for i in range(0,nlig):
      for j in range(0,len(nsubs)):
        fp.write('%2d ' % ind[i,j])
      fp.write('%8.3f +/- %5.3f\n' % (Value[i],Error[i]))

    fp.close()

  # Pkeep=np.sum(PkeepA,axis=0)
  for i in range(0,NF):
    for j in range(0,nlig):
      # Actual expressions
      # G[i,j]=Eall[ncentral,i,j]-kT*np.log(np.sum(PkeepA[:,i,j]*np.exp(-Eshift[ncentral,i,j]/kT)/np.sum(np.exp(lndenom[:,i,j]))))
      # More efficient? - exp(-Eshift(central)/kT)=1`
      G[i,j]=-Eall[ncentral,i,j]-kT*np.log(np.sum(PkeepA[:,i,j]))+kT*np.log(np.sum(np.exp(lndenom[:,i,j])))

  np.savetxt('G.dat',G)

  Gmin=np.min(G,axis=0)
  Value=Gmin-kT*np.log(np.mean(np.exp(-(G-Gmin)/kT),axis=0))
  Value-=Value[0]

  np.random.seed(2401)
  GS=np.zeros((NBS,nlig))
  for i in range(0,NBS):
    GS[i,:]=Gmin-kT*np.log(np.mean(np.exp(-(G[np.random.randint(0,NF,(NF,1)),:]-Gmin)/kT),axis=0))
  Error=np.std(GS,axis=0)

  fp=open('Result.txt','w')

  for i in range(0,nlig):
    for j in range(0,len(nsubs)):
      fp.write('%2d ' % ind[i,j])
    fp.write('%8.3f +/- %5.3f\n' % (Value[i],Error[i]))

  fp.close()
