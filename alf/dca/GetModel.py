#! /usr/bin/env python

def GetModelDCA(alf_info,NF,Path,NBS=50):
  """
  Computes the h and J fields and couplings from Potts model estimator

  Used in alf.FinishDCA. See further documentation there.

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

  import sys, os
  import numpy as np
  import subprocess
  import copy

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
  if alf_info['bias']=='bcxstu2026':
    t=np.loadtxt('t_prev.dat')
    u=np.loadtxt('u_prev.dat')
  kT=0.001987*alf_info['temp']

  if alf_info['bias']=='bcxs2018':
    alpha=0.017
  elif alf_info['bias']=='bcxstu2026':
    alpha=0.012
  c1=1
  x1=1-np.exp(-5.56*1)
  s1=1/(1+alpha)
  if alf_info['bias']=='bcxstu2026':
    t1=1/(1-(1+alpha))
    u1=1/(1+alpha)

  h_bias=np.zeros((nblocks,))
  J_bias=np.zeros((nblocks,nblocks))
  # E=np.dot(LList,b)+np.dot(np.dot(LList,c),LList)+np.dot(np.dot(1-np.exp(-5.56*LList),x),LList)+np.dot(np.dot(LList/(LList+0.017),s),LList)
  for i in range(nblocks):
    h_bias[i]=b[i]
    for j in range(nblocks):
      if (block2site[i]!=block2site[j]):
        J_bias[i,j]=c1*(c[i,j]+c[j,i])+x1*(x[i,j]+x[j,i])+s1*(s[i,j]+s[j,i])
        if alf_info['bias']=='bcxstu2026':
          J_bias[i,j]=J_bias[i,j]+t1*(t[i,j]+t[j,i])+u1*(u[i,j]+u[j,i])
  np.savetxt(Path+'/h.bias.dat',h_bias,fmt=' %10.6f')
  np.savetxt(Path+'/J.bias.dat',J_bias,fmt=' %10.6f')

  iBStrList=['']
  for iB in range(NBS):
    iBStrList.append('bs'+str(iB)+'.')

  for iBStr in iBStrList:
    h_fnm=Path+'/h.'+iBStr+tag+'.dat'
    J_fnm=Path+'/J.'+iBStr+tag+'.dat'
    m1_fnm=Path+'/m1.'+iBStr+'obs.dat'
    m2_fnm=Path+'/m2.'+iBStr+'obs.dat'
    h=np.loadtxt(h_fnm)
    J=np.loadtxt(J_fnm)
    m1=np.loadtxt(m1_fnm)
    m2=np.loadtxt(m2_fnm)
    hMask=np.loadtxt(h_fnm)
    hMask[m1==0]=np.nan
    JMask=np.loadtxt(J_fnm)
    JMask[m2==0]=np.nan
    block0=0
    for i in range(nsites):
      h=np.delete(h,block0,axis=0)
      J=np.delete(J,block0,axis=0)
      J=np.delete(J,block0,axis=1)
      hMask=np.delete(hMask,block0,axis=0)
      JMask=np.delete(JMask,block0,axis=0)
      JMask=np.delete(JMask,block0,axis=1)
      JMask[block0:block0+nsubs[i],block0:block0+nsubs[i]]=0
      block0+=nsubs[i]
    h*=-kT
    J*=-kT
    hMask*=-kT
    JMask*=-kT
    np.savetxt(Path+'/hGuess.'+iBStr+'model.dat',h+h_bias,fmt=' %10.6f')
    np.savetxt(Path+'/JGuess.'+iBStr+'model.dat',J+J_bias,fmt=' %10.6f')
    np.savetxt(Path+'/h.'+iBStr+'model.dat',hMask+h_bias,fmt=' %10.6f')
    np.savetxt(Path+'/J.'+iBStr+'model.dat',JMask+J_bias,fmt=' %10.6f')

