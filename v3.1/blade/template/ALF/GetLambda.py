#! /usr/bin/env python

import sys
import numpy as np
from xdrlib import Unpacker
from xdrlib import Packer

if len(sys.argv) < 3:
  print("Error: need an input and output filename")
  quit()

nblocks=np.loadtxt('../nblocks',dtype='int')
Lambdas=np.zeros((0,nblocks))

p=Packer()
p.pack_int(0)
for j in range(0,nblocks):
  p.pack_float(0)
linewidth=len(p.get_buffer())

for ifp in range(2,len(sys.argv)):
  # Lambda=np.loadtxt(sys.argv[ifp])
  fp=open(sys.argv[ifp],"rb")
  fpdata=fp.read()
  lines=len(fpdata)//linewidth
  fp.close()
  Lambda=np.zeros((lines,nblocks))
  p=Unpacker(fpdata)
  for i in range(0,lines):
    p.unpack_int()
    for j in range(0,nblocks):
      Lambda[i,j]=p.unpack_float()
  Lambdas=np.concatenate((Lambdas,Lambda),axis=0)

np.savetxt(sys.argv[1],Lambdas,fmt="%10.6f")
