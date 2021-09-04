#! /usr/bin/env python

import numpy as np

# ddg=0.688*np.array([0, 0, 1,-1,-1, 0, 0, 1, 0, 0])
ddg=0.688*np.array([0, 0, 0, 1,-1,-1,-1, 0, 0, 0, 1, 0, 0, 0, 0])
i1=np.zeros((15,))
nsites=15

fp=open("Result.txt","r")
fpout=open("Corrected.txt","w")

for line in fp:
  linesplit=line.split()
  for i in range(0,nsites):
    i1[i]=int(linesplit[i])
  V=float(linesplit[nsites])+np.dot(i1,ddg)
  E=float(linesplit[nsites+2])
  print(np.dot(i1,ddg))
  
  for i in range(0,nsites):
    fpout.write("%2d " % (i1[i],))
  fpout.write("%8.3f +/- %5.3f\n" % (V,E))

fp.close()
fpout.close()
