#! /usr/bin/env python

def GetModelDCA(alf_info,NF,Path):
  import sys, os
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

  nsubs=alf_info['nsubs']
  nblocks=alf_info['nblocks']
  # nreps=np.loadtxt('../nreps',dtype='int')
  # ncentral=np.loadtxt('../ncentral',dtype='int')

  nsites=len(nsubs)
  nblocks=np.sum(nsubs)

  block2site=np.zeros((nblocks,),dtype='int')
  k=0;
  for i in range(nsites):
    for j in range(nsubs[i]):
      block2site[k]=i
      k+=1

  b=np.loadtxt('b_prev.dat')
  if os.path.isfile('b_corr.dat'):
    b=b+np.loadtxt('b_corr.dat')
  c=np.loadtxt('c_prev.dat')
  x=np.loadtxt('x_prev.dat')
  s=np.loadtxt('s_prev.dat')
  kT=0.001987*alf_info['temp']

  c0=0
  c1=1
  x0=1-np.exp(-5.56*0)
  x1=1-np.exp(-5.56*1)
  s0=0/(0+0.017)
  s1=1/(1+0.017)

  h_bias=np.zeros((nblocks,))
  J_bias=np.zeros((nblocks,nblocks))
  # E=np.dot(LList,b)+np.dot(np.dot(LList,c),LList)+np.dot(np.dot(1-np.exp(-5.56*LList),x),LList)+np.dot(np.dot(LList/(LList+0.017),s),LList)
  for i in range(nblocks):
    h_bias[i]=b[i]
    for j in range(nblocks):
      if (block2site[i]!=block2site[j]):
        J_bias[i,j]=c1*(c[i,j]+c[j,i])+x1*(x[i,j]+x[j,i])+s1*(s[i,j]+s[j,i])
  np.savetxt(Path+'/h.bias.dat',h_bias,fmt=' %10.6f')
  np.savetxt(Path+'/J.bias.dat',J_bias,fmt=' %10.6f')

  h_fnm=Path+'/h.'+tag+'.dat'
  J_fnm=Path+'/J.'+tag+'.dat'
  h=np.loadtxt(h_fnm)
  J=np.loadtxt(J_fnm)
  block0=0
  for i in range(nsites):
    h=np.delete(h,block0,axis=0)
    J=np.delete(J,block0,axis=0)
    J=np.delete(J,block0,axis=1)
    block0+=nsubs[i]
  h*=-kT
  J*=-kT
  np.savetxt(Path+'/h.model.dat',h+h_bias,fmt=' %10.6f')
  np.savetxt(Path+'/J.model.dat',J+J_bias,fmt=' %10.6f')

  for iB in range(NB):
    print(iB)
    h_fnm=Path+'/h.bs'+str(iB)+'.'+tag+'.dat'
    J_fnm=Path+'/J.bs'+str(iB)+'.'+tag+'.dat'
    h=np.loadtxt(h_fnm)
    J=np.loadtxt(J_fnm)
    block0=0
    for i in range(nsites):
      h=np.delete(h,block0,axis=0)
      J=np.delete(J,block0,axis=0)
      J=np.delete(J,block0,axis=1)
      block0+=nsubs[i]
    h*=-kT
    J*=-kT
    np.savetxt(Path+'/h.bs'+str(iB)+'.model.dat',h+h_bias,fmt=' %10.6f')
    np.savetxt(Path+'/J.bs'+str(iB)+'.model.dat',J+J_bias,fmt=' %10.6f')
