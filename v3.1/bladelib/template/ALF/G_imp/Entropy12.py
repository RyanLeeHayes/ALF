#! /usr/bin/env python

import sys, os
import numpy as np

Nmc=5000000
Nl=400
c=5.5
cutlsum=0.8
dl=1.0/Nl

Ndim=int(sys.argv[1])

histogram_a=np.zeros((Nl,1))

for i in range(0,Nmc):
  ti=np.random.random_sample((Ndim,))
  unli=np.exp(c*np.sin(np.pi*ti-np.pi/2))
  li=unli/np.sum(unli)

  for j in range(0,Ndim):
    for k in range(j+1,Ndim):
      lisum=li[j]+li[k]
      if lisum > cutlsum:
        ind=np.floor((li[j]/lisum)*Nl).astype('int')
        histogram_a[ind,0]+=1
histogram=histogram_a+histogram_a[::-1]

S_k=np.log(histogram)
G_k=-S_k+np.log(np.mean(histogram))

# plot(dl/2:dl:1,histogram1)
# plot((dl/2):dl:1,G12_k')

np.savetxt('G12_'+str(Ndim)+'.dat',G_k)
