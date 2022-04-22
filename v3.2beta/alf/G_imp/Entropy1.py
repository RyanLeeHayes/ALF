#! /usr/bin/env python

import sys, os
import numpy as np

Nmc=5000000
Nl=400
c=5.5
dl=1.0/Nl

Ndim=int(sys.argv[1])

histogram=np.zeros((Nl,1))

for i in range(0,Nmc):
  ti=np.random.random_sample((Ndim,))
  unli=np.exp(c*np.sin(np.pi*ti-np.pi/2))
  li=unli/np.sum(unli)

  indi=np.floor(li*Nl).astype('int')
  for j in range(0,Ndim):
    histogram[indi[j],0]+=1

S_k=np.log(histogram)
G_k=-S_k+np.log(np.mean(histogram))

# plot(dl/2:dl:1,histogram1)
# plot((dl/2):dl:1,G12_k')

np.savetxt('G1_'+str(Ndim)+'.dat',G_k)
