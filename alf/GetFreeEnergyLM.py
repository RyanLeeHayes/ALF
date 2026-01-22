#! /usr/bin/env python

def GetFreeEnergyLM(alf_info,ms,msprof):
  """
  Unpack the output of lmalf to determine optimal bias changes

  Used with alf_info['loss']=='nonlinear2024'

  RunLMALF.sh produces an output file OUT.dat in analysis[i] of all the
  optimal changes in the bias parameters. This program should be run from
  analysis[i] and upacks that list into matrices for the changes in each
  bias term b.dat, c.dat, x.dat, s.dat, etc. The program also caps changes
  in the bias potentials to a specific limit, and scales all bias changes
  down if changes exceed this limit. SetVars.py then takes the output and
  adds it to the current biasing potentials and translates it into a
  format readable by the molecular dynamics engine.

  Parameters
  ----------
  alf_info : dict
      Dictionary of variables alf needs to run
  ms : int
      Flag for whether to include intersite biases. 0 for no, 1 for c, x,
      and s biases, 2 for just c biases. Typically taken from the first
      element of ntersite list.
  msprof : int
      Flag for whether to include intersite profiles. 0 for no, 1 for yes.
      Typically taken from the second element of ntersite list.
  """

  import sys
  import numpy as np

  kT=0.001987*alf_info['temp']

  cutb=2
  cutc=8
  cutx=2
  cuts=1
  if alf_info['bias']=='bcxstu2026':
    cutt=1
    cutu=1
  cutc2=2
  cutx2=0.5
  cuts2=0.5
  if alf_info['bias']=='bcxstu2026':
    cutt2=0.05
    cutu2=0.05

  nsubs=alf_info['nsubs']
  nblocks=alf_info['nblocks']

  b_prev=np.loadtxt('b_prev.dat')
  c_prev=np.loadtxt('c_prev.dat')
  x_prev=np.loadtxt('x_prev.dat')
  s_prev=np.loadtxt('s_prev.dat')
  if alf_info['bias']=='bcxstu2026':
    t_prev=np.loadtxt('t_prev.dat')
    u_prev=np.loadtxt('u_prev.dat')

  b=np.zeros((1,nblocks))
  c=np.zeros((nblocks,nblocks))
  x=np.zeros((nblocks,nblocks))
  s=np.zeros((nblocks,nblocks))
  if alf_info['bias']=='bcxstu2026':
    t=np.zeros((nblocks,nblocks))
    u=np.zeros((nblocks,nblocks))

  nparm=0
  for isite in range(0,len(nsubs)):
    n1=nsubs[isite]
    n2=nsubs[isite]*(nsubs[isite]-1)//2;
    for jsite in range(isite,len(nsubs)):
      n3=nsubs[isite]*nsubs[jsite]
      if isite==jsite:
        nparm+=n1+5*n2
        if alf_info['bias']=='bcxstu2026':
          nparm+=4*n2
      elif ms==1:
        nparm+=5*n3
        if alf_info['bias']=='bcxstu2026':
          nparm+=4*n3
      elif ms==2:
        nparm+=n3

  cutlist=np.zeros((nparm,))
  n0=0
  iblock=0
  for isite in range(0,len(nsubs)):
    jblock=iblock
    n1=nsubs[isite]
    n2=nsubs[isite]*(nsubs[isite]-1)//2;
    for jsite in range(isite,len(nsubs)):
      n3=nsubs[isite]*nsubs[jsite]
      if isite==jsite:
        for i in range(0,nsubs[isite]):
          cutlist[n0:n0+1]=cutb
          n0+=1
          for j in range(i+1,nsubs[isite]):
            cutlist[n0:n0+1]=cutc
            n0+=1
            cutlist[n0:n0+2]=cutx
            n0+=2
            cutlist[n0:n0+2]=cuts
            n0+=2
            if alf_info['bias']=='bcxstu2026':
              cutlist[n0:n0+2]=cutt
              n0+=2
              cutlist[n0:n0+2]=cutu
              n0+=2
      elif ms>0:
        for i in range(0,nsubs[isite]):
          for j in range(0,nsubs[jsite]):
            cutlist[n0:n0+1]=cutc2
            n0+=1
            if ms==1:
              cutlist[n0:n0+2]=cutx2
              n0+=2
              cutlist[n0:n0+2]=cuts2
              n0+=2
              if alf_info['bias']=='bcxstu2026':
                cutlist[n0:n0+2]=cutt2
                n0+=2
                cutlist[n0:n0+2]=cutu2
                n0+=2
      jblock+=nsubs[jsite]
    iblock+=nsubs[isite]

  # coeff=C^-1*V;
  coeff=np.loadtxt('OUT.dat')

  coeff_orig=coeff

  scaling=1.5/np.max(np.abs(coeff[0:n0]/cutlist))
  if scaling>1:
    scaling=1
  coeff*=scaling

  print("scaling is:")
  print(scaling)

  ind=0
  iblock=0
  for isite in range(0,len(nsubs)):
    jblock=iblock
    for jsite in range(isite,len(nsubs)):
      if isite==jsite:
        for i in range(0,nsubs[isite]):
          b[0,iblock+i]=coeff[ind]
          ind+=1
          for j in range(i+1,nsubs[isite]):
            c[iblock+i,jblock+j]=coeff[ind]
            ind+=1
            x[iblock+i,jblock+j]=coeff[ind]
            ind+=1
            x[jblock+j,iblock+i]=coeff[ind]
            ind+=1
            s[iblock+i,jblock+j]=coeff[ind]
            ind+=1
            s[jblock+j,iblock+i]=coeff[ind]
            ind+=1
            if alf_info['bias']=='bcxstu2026':
              t[iblock+i,jblock+j]=coeff[ind]
              ind+=1
              t[jblock+j,iblock+i]=coeff[ind]
              ind+=1
              u[iblock+i,jblock+j]=coeff[ind]
              ind+=1
              u[jblock+j,iblock+i]=coeff[ind]
              ind+=1
      elif ms>0:
        for i in range(0,nsubs[isite]):
          for j in range(0,nsubs[jsite]):
            c[iblock+i,jblock+j]=coeff[ind]
            ind+=1
            if ms==1:
              x[iblock+i,jblock+j]=coeff[ind]
              ind+=1
              x[jblock+j,iblock+i]=coeff[ind]
              ind+=1
              s[iblock+i,jblock+j]=coeff[ind]
              ind+=1
              s[jblock+j,iblock+i]=coeff[ind]
              ind+=1
              if alf_info['bias']=='bcxstu2026':
                t[iblock+i,jblock+j]=coeff[ind]
                ind+=1
                t[jblock+j,iblock+i]=coeff[ind]
                ind+=1
                u[iblock+i,jblock+j]=coeff[ind]
                ind+=1
                u[jblock+j,iblock+i]=coeff[ind]
                ind+=1
      jblock+=nsubs[jsite]
    iblock+=nsubs[isite]

  iblock=0
  for isite in range(0,len(nsubs)):
    jblock=iblock
    for jsite in range(isite,len(nsubs)):
      if isite!=jsite:
        for i in range(0,nsubs[isite]):
          b[0,iblock+i]+=c[iblock+i,jblock]
          c[iblock+i,jblock:jblock+nsubs[jsite]]-=c[iblock+i,jblock]
        for j in range(0,nsubs[jsite]):
          b[0,jblock+j]+=c[iblock,jblock+j]
          c[iblock:iblock+nsubs[isite],jblock+j]-=c[iblock,jblock+j]
      jblock+=nsubs[jsite]
    iblock+=nsubs[isite]

  iblock=0
  for isite in range(0,len(nsubs)):
    b[0,iblock:iblock+nsubs[isite]]-=b[0,iblock]
    iblock+=nsubs[isite]


  np.savetxt('b.dat',b,fmt=' %7.2f')
  np.savetxt('c.dat',c,fmt=' %7.2f')
  np.savetxt('x.dat',x,fmt=' %7.2f')
  np.savetxt('s.dat',s,fmt=' %7.2f')
  if alf_info['bias']=='bcxstu2026':
    np.savetxt('t.dat',t,fmt=' %7.2f')
    np.savetxt('u.dat',u,fmt=' %7.2f')
