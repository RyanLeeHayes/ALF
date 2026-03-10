#! /usr/bin/env python

def GetVarianceDCA(alf_info,NF,Path,NBS=50):
  """
  Computes free energies and places them in Results.txt

  Will bail out if there are more than 2^20 alchemical states. Used in
  alf.FinishDCA. See further documentation there. Must be run from inside
  dca analysis directory.

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

  import sys, os, os.path
  import numpy as np
  import subprocess
  import copy

  if os.path.exists(Path+'/h.LM.dat'):
    tag='LM'
  elif os.path.exists(Path+'/h.PLM.dat'):
    tag='PLM'
  else:
    print('Error: Neither h.LM.dat or h.PLM.dat found.')
    quit(1)

  nblocks=alf_info['nblocks']+0 # Add 0 so python makes a copy
  nsubs=alf_info['nsubs']+0 # Add 0 so python makes a copy
  # nreps=np.loadtxt('../nreps',dtype='int')
  # ncentral=np.loadtxt('../ncentral',dtype='int')

  if (np.prod(nsubs)>1024*1024):
    print("Too many states")
    quit(1)

  nblocks+=len(nsubs)
  nsubs+=1
  nlig=np.prod(nsubs)
  nlig_ng=np.prod(nsubs-1)
  # _ng stands for no gaps

  b=np.loadtxt('b_prev.dat')
  if os.path.isfile('b_corr.dat'):
    b=b+np.loadtxt('b_corr.dat')
  c=np.loadtxt('c_prev.dat')
  x=np.loadtxt('x_prev.dat')
  s=np.loadtxt('s_prev.dat')
  if alf_info['bias']=='bcxstu2026':
    t=np.loadtxt('t_prev.dat')
    u=np.loadtxt('u_prev.dat')
  kT=0.001987*alf_info['temp']

  if alf_info['bias']=='bcxs2018':
    alpha=0.017
  elif alf_info['bias']=='bcxstu2026':
    alpha=0.012

  G=np.zeros((nlig_ng,))
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

  ind_ng=np.zeros((nlig_ng,len(nsubs)),dtype='int') # index, no gaps
  for i in range(1,nlig_ng):
    ind_ng[i,:]=ind_ng[i-1,:]
    for j in range(len(nsubs)-1,-1,-1):
      if (ind_ng[i,j]+1)<(nsubs[j]-1):
        ind_ng[i,j]+=1
        break
      else:
        ind_ng[i,j]=0
  blk_ng=copy.deepcopy(ind_ng) # block, no gaps
  for i in range(1,len(nsubs)):
    blk_ng[:,i:]+=nsubs[i-1]-1

  h_fnm=Path+'/h.'+tag+'.dat'
  J_fnm=Path+'/J.'+tag+'.dat'
  h=np.loadtxt(h_fnm)
  J=np.loadtxt(J_fnm)
  jno0=0
  for j in range(nlig):
    # P=np.exp(np.sum(h[blk[j,:]])+0.5*np.sum(np.sum(J[blk[j,:],blk[j,:]])))
    Epotts=np.sum(h[blk[j,:]])+0.5*np.sum(np.sum(J[blk[j,:]][:,blk[j,:]]))
    if np.all(ind[j,:]>0):
      LList=np.zeros((nblocks-len(nsubs),))
      LList[blk_ng[jno0,:]]=1
      E=np.dot(LList,b)+np.dot(np.dot(LList,c),LList)+np.dot(np.dot(1-np.exp(-5.56*LList),x),LList)+np.dot(np.dot(LList/(LList+alpha),s),LList)
      if alf_info['bias']=='bcxstu2026':
        E=E+np.dot(np.dot(LList/(LList-(1+alpha)),t),LList)+np.dot(np.dot(LList/(LList+alpha),u),LList**2)
      G[jno0]=E-kT*Epotts
      jno0+=1

  GS=np.zeros((NBS,nlig_ng))
  # GS=np.zeros((NBS,))
  Value=np.zeros((nlig_ng,))
  Error=np.zeros((nlig_ng,))
  h=np.zeros((1,nblocks,NBS))
  J=np.zeros((nblocks,nblocks,NBS))
  for i in range(NBS):
    print(i)
    h_fnm=Path+'/h.bs'+str(i)+'.'+tag+'.dat'
    J_fnm=Path+'/J.bs'+str(i)+'.'+tag+'.dat'
    h[:,:,i]=np.loadtxt(h_fnm)
    J[:,:,i]=np.loadtxt(J_fnm)

  jno0=0
  for j in range(nlig):
    if np.all(ind[j,:]>0):
      LList=np.zeros((nblocks-len(nsubs),))
      LList[blk_ng[jno0,:]]=1
      E=np.dot(LList,b)+np.dot(np.dot(LList,c),LList)+np.dot(np.dot(1-np.exp(-5.56*LList),x),LList)+np.dot(np.dot(LList/(LList+alpha),s),LList)
      if alf_info['bias']=='bcxstu2026':
        E=E+np.dot(np.dot(LList/(LList-(1+alpha)),t),LList)+np.dot(np.dot(LList/(LList+alpha),u),LList**2)
      for i in range(NBS):
        # P=np.exp(np.sum(h[blk[j,:]])+0.5*np.sum(np.sum(J[blk[j,:],blk[j,:]])))
        Epotts=np.sum(h[0,blk[j,:],i])+0.5*np.sum(np.sum(J[blk[j,:]][:,blk[j,:],i]))
        GS[i,jno0]=E-kT*Epotts
      # Value[jno0]=G[jno0]-G[0]
      # Error[jno0]=np.sqrt(np.mean((GS-G[jno0])**2))
      Error[jno0]=np.sqrt(np.mean((GS[:,jno0]-G[jno0])**2,axis=0))
      jno0+=1
  Value=G-G[0]

  np.savetxt('G.dat',G)
  np.savetxt('GS.dat',GS)

  fp=open('ResultGuess.txt','w')

  for i in range(0,nlig_ng):
    for j in range(0,len(nsubs)):
      fp.write('%2d ' % ind_ng[i,j])
    fp.write('%8.3f +/- %5.3f\n' % (Value[i],Error[i]))

  fp.close()

def GetVarianceMaskedDCA(alf_info,Path):
  """
  Takes DCA free energies in Results.txt and masks them with nan if a
  combination was unsampled

  Parameters
  ----------
  alf_info : dict
      Dictionary of variables alf needs to run
  Path : str
      Path to a data directory for Potts model estimation, typically 'data'
  """

  import numpy as np

  nsubs=alf_info['nsubs']
  m1=np.loadtxt(Path+'/m1.obs.dat')
  m2=np.loadtxt(Path+'/m2.obs.dat')

  fpin=open('ResultGuess.txt','r')
  fpout=open('Result.txt','w')

  ind=np.zeros((len(nsubs),),dtype=int)
  for line in fpin:
    for i in range(len(nsubs)):
      ind[i]=int(line.split()[i])
    shift=0
    iblock0=0
    for i in range(len(nsubs)):
      iind=ind[i]+iblock0+1
      if m1[iind]==0:
        shift=np.nan
      jblock0=iblock0+nsubs[i]+1
      for j in range(i+1,len(nsubs)):
        jind=ind[j]+jblock0+1
        if m2[iind,jind]==0:
          shift=np.nan
        jblock0=jblock0+nsubs[j]+1
      iblock0=iblock0+nsubs[i]+1
    for i in range(len(nsubs)):
      fpout.write('%2d ' % ind[i])
    Value=float(line.split()[len(nsubs)])+shift
    Error=float(line.split()[len(nsubs)+2])+shift
    fpout.write('%8.3f +/- %5.3f\n' % (Value,Error))
  fpin.close()
  fpout.close()

