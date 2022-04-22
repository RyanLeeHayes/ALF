#! /usr/bin/env python

def GetLambdas(alf_info,istep,ndupl=None,begres=None,endres=None):
  import sys, os
  import numpy as np
  # from subprocess import call
  from alf.GetLambda import GetLambda

  if ndupl==None:
    production=False
    # istep=int(sys.argv[1])
    ndupl=1
  else:
    production=True
    # istep=int(sys.argv[1])
    # ndupl=int(sys.argv[2])
    # begres=int(sys.argv[3])
    # endres=int(sys.argv[4])
  # else:
    # print("Error: Need 1 argument for flattening or 4 arguments for production")
    # quit()

  nblocks=alf_info['nblocks']
  nsubs=alf_info['nsubs']
  nreps=alf_info['nreps']
  name=alf_info['name']

  if not os.path.isdir('data'):
    os.mkdir('data')
  DIR="data"

  # ----------------------------------------------------------------------------

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
      GetLambda(alf_info,fnmout,fnmsin)
