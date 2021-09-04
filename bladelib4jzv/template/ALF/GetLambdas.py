#! /usr/bin/env python

import sys, os
import numpy as np
from subprocess import call

if len(sys.argv)==2:
  production=False
  istep=int(sys.argv[1])
  ndupl=1
elif len(sys.argv)==5:
  production=True
  istep=int(sys.argv[1])
  ndupl=int(sys.argv[2])
  begres=int(sys.argv[3])
  endres=int(sys.argv[4])
else:
  print("Error: Need 1 argument for flattening or 4 arguments for production")
  quit()

nblocks=np.loadtxt('../nblocks',dtype='int')
nsubs=np.loadtxt('../nsubs',dtype='int',ndmin=1)
nreps=np.loadtxt('../nreps',dtype='int')

fp=open('../name','r')
name=fp.readline().strip()
fp.close()

if not os.path.isdir('data'):
  os.mkdir('data')
DIR="data"

# ------------------------------------------------------------------------------

alphabet='abcdefghijklmnopqrstuvwxyz'

for idupl in range(0,ndupl):

  DDIR='../run'+str(istep)
  if production:
    DDIR=DDIR+alphabet[idupl]

  for i in range(0,nreps):
    fnmsin=[]
    if nreps>1:
      reptag="_"+str(i)
    else:
      reptag=""
    if production:
      for j in range(begres,endres):
        fnmsin.append(DDIR+'/res/'+name+'_prod'+str(j+1)+'.lmd'+reptag)
    else:
      fnmsin.append(DDIR+'/res/'+name+'_flat.lmd'+reptag)
    fnmout=DIR+("/Lambda.%d.%d.dat" % (idupl,i))
    call(['../ALF/GetLambda.py',fnmout]+fnmsin)
