#! /usr/bin/env python

def GetVarianceDCA(alf_info,NF,Path):
  import sys, os, os.path
  import numpy as np
  import subprocess
  import copy

  NB=50

  if os.path.exists(Path+'/h.LM.dat'):
    tag='LM'
  elif os.path.exists(Path+'/h.PLM.dat'):
    tag='PLM'
  else:
    print('Error: Neither h.LM.dat or h.PLM.dat found.')
    quit()

  nblocks=alf_info['nblocks']+0 # Add 0 so python makes a copy
  nsubs=alf_info['nsubs']+0 # Add 0 so python makes a copy
  # nreps=np.loadtxt('../nreps',dtype='int')
  # ncentral=np.loadtxt('../ncentral',dtype='int')

  if (np.prod(nsubs)>1024*1024):
    print("Too many states")
    quit()

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
  kT=0.001987*alf_info['temp']

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
      E=np.dot(LList,b)+np.dot(np.dot(LList,c),LList)+np.dot(np.dot(1-np.exp(-5.56*LList),x),LList)+np.dot(np.dot(LList/(LList+0.017),s),LList)
      G[jno0]=E-kT*Epotts
      jno0+=1

  GS=np.zeros((NB,nlig_ng))
  # GS=np.zeros((NB,))
  Value=np.zeros((nlig_ng,))
  Error=np.zeros((nlig_ng,))
  h=np.zeros((1,nblocks,NB))
  J=np.zeros((nblocks,nblocks,NB))
  for i in range(NB):
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
      E=np.dot(LList,b)+np.dot(np.dot(LList,c),LList)+np.dot(np.dot(1-np.exp(-5.56*LList),x),LList)+np.dot(np.dot(LList/(LList+0.017),s),LList)
      for i in range(NB):
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

  fp=open('Result.txt','w')

  for i in range(0,nlig_ng):
    for j in range(0,len(nsubs)):
      fp.write('%2d ' % ind_ng[i,j])
    fp.write('%8.3f +/- %5.3f\n' % (Value[i],Error[i]))

  fp.close()
