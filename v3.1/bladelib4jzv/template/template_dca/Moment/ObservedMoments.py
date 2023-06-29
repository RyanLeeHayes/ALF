#! /usr/bin/env python

import numpy as np
import sys, os

nsubs=np.loadtxt('../../nsubs',dtype='int',ndmin=1)
nblocks=np.loadtxt('../../nblocks',dtype='int')
nsites=np.size(nsubs,0)

NF=int(sys.argv[1])
Path=sys.argv[2]

nsubs+=1
nblocks+=nsites

m1=np.zeros((1,nblocks))
m2=np.zeros((nblocks,nblocks))

m1mean=np.zeros((1,nblocks))
m2mean=np.zeros((nblocks,nblocks))

for ifile in range(NF):
  F=np.loadtxt(Path+'/Filter.'+str(ifile)+'.dat',dtype='int')

  Z=np.size(F,0);
  block0i=0
  for i in range(0,nsites):
    for Ai in range(0,nsubs[i]):
      m1[0,block0i+Ai]=np.sum(F[:,i]==Ai)*(1.0/Z)
      block0j=0
      for j in range(0,nsites):
        for Aj in range(0,nsubs[j]):
          print(ifile,block0i+Ai,block0j+Aj)
          m2[block0i+Ai,block0j+Aj]=np.sum(np.logical_and(F[:,i]==Ai,F[:,j]==Aj))*(1.0/Z)
        block0j+=nsubs[j]
    block0i+=nsubs[i]

  np.savetxt(Path+'/m1.'+str(ifile)+'.obs.dat',m1)
  np.savetxt(Path+'/m2.'+str(ifile)+'.obs.dat',m2)

  m1mean+=m1
  m2mean+=m2

np.savetxt(Path+'/m1.obs.dat',m1mean/NF)
np.savetxt(Path+'/m2.obs.dat',m2mean/NF)
