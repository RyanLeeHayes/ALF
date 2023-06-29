#! /usr/bin/env python

import sys
import numpy as np
from xdrlib import Unpacker
from xdrlib import Packer

if len(sys.argv) < 2:
  print("Error: need an input filename")
  quit()

nblocks=np.loadtxt('../nblocks',dtype='int')

p=Packer()
p.pack_int(0)
for j in range(0,nblocks):
  p.pack_float(0)
linewidth=len(p.get_buffer())

try:
  fp=open(sys.argv[1],"rb")
  fpdata=fp.read()
  lines=len(fpdata)//linewidth
  fp.close()
except:
  try:
    fp=open(sys.argv[1]+'_0',"rb")
    fpdata=fp.read()
    lines=len(fpdata)//linewidth
    fp.close()
  except:
    print(0)
    quit()

print(lines)
