#! /usr/bin/env python

import numpy as np
import sys, os

NF=int(sys.argv[1])
Path=sys.argv[2]
PathOut=sys.argv[3]

nblocks=np.loadtxt('../../nblocks',dtype='int')
nsubs=np.loadtxt('../../nsubs',dtype='int',ndmin=1)
nreps=np.loadtxt('../../nreps',dtype='int')
ncentral=np.loadtxt('../../ncentral',dtype='int')

for ifile in range(NF):
  L=np.loadtxt(Path+'/Lambda.'+str(ifile)+'.'+str(ncentral)+'.dat')

  ind=np.zeros((np.size(L,0),len(nsubs)),dtype='int')
  block0=0
  for i in range(0,len(nsubs)):
    print 'i is',i
    for j in range(0,nsubs[i]):
      ind[L[:,block0+j]>0.99,i]=j+1
    block0+=nsubs[i]

  np.savetxt(PathOut+'/Filter.'+str(ifile)+'.dat',ind,'%1d')

h=np.zeros((nblocks+len(nsubs),))
J=np.zeros((nblocks+len(nsubs),nblocks+len(nsubs)))
np.savetxt(PathOut+'/h.null.dat',h,'%1d')
np.savetxt(PathOut+'/J.null.dat',J,'%1d')
