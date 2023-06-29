#! /usr/bin/env python

import numpy as np

i1=np.zeros((15,),dtype=int)
nsites=15

fp=open("Corrected.txt","r")

Gsize=tuple(2*np.ones((nsites,),dtype=int))
G=np.zeros(Gsize)
for line in fp:
  linesplit=line.split()
  for i in range(0,nsites):
    i1[i]=int(linesplit[i])
  G[tuple(i1)]=float(linesplit[nsites])
  # E=float(linesplit[nsites+2])
G=G-G[tuple(np.zeros((nsites,),dtype=int))]

fp.close()

h=np.zeros((nsites,))
for i in range(0,nsites):
  i1=np.zeros((nsites,),dtype=int)
  i1[i]=1
  h[i]=G[tuple(i1)]

J=np.zeros((nsites,nsites))
for i in range(0,nsites):
  for j in range(i+1,nsites):
    i1=np.zeros((nsites,),dtype=int)
    i1[i]=1
    i1[j]=1
    J[i,j]=G[tuple(i1)]-h[i]-h[j]

mh=np.zeros((nsites,))
for i in range(0,nsites):
  G0=np.mean(np.take(G,0,axis=i))
  G1=np.mean(np.take(G,1,axis=i))
  mh[i]=G1-G0
  
np.savetxt("Fields.txt",h,fmt="%8.3f")
np.savetxt("Couplings.txt",J,fmt="%8.3f")
np.savetxt("MeanFields.txt",mh,fmt="%8.3f")
