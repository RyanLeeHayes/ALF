#! /usr/bin/env python

import sys, os
import numpy as np

Nmc=5000000
Nl=20
c=5.5
dl=1.0/Nl

Ndim=int(sys.argv[1])

histogram_a=np.zeros((Nl,Nl))

for i in range(0,Nmc):
  ti=np.random.random_sample((Ndim,))
  unli=np.exp(c*np.sin(np.pi*ti-np.pi/2))
  li=unli/np.sum(unli)

  indi=np.floor(li*Nl).astype('int')
  for j in range(0,Ndim):
    for k in range(j+1,Ndim):
      histogram_a[indi[j],indi[k]]+=1
histogram=histogram_a+histogram_a.T

S_k=np.log(histogram)
G_k=-S_k+np.log(np.mean(histogram))

# plot(dl/2:dl:1,histogram1)
# plot((dl/2):dl:1,G12_k')

np.savetxt('G2_'+str(Ndim)+'.dat',G_k)
