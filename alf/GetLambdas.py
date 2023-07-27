#! /usr/bin/env python

def GetLambdas(alf_info,istep,ndupl=None,begres=None,endres=None):
  """
  Reads alchemical trajectories from binary format

  This routine reads binary alchemical flattening trajectories from
  run[i]/res/[name]_flat.lmd or binary alchemical production trajectories
  from run[i][a]/res/[name]_prod[itt].lmd where [i] is the cycle number,
  [a] is the duplicate letter, [name] is the system name, and [itt] is the
  production chunk, and copies them into human readable trajectories in
  analysis[i]/data/Lambda.[ia].[ir].dat where [ia] is the duplicate index
  and [ir] is the replica index. This routine should be called from the
  analysis[i] directory.

  This routine can be called during flattening or production. Flattening
  versus production is detected by the absence or presence, respectively
  of the three optional parameters.

  Parameters
  ----------
  alf_info : dict
      Dictionary of variables alf needs to run
  istep : int
      The current cycle of alf being analyzed
  ndupl : int, optional
      The number of independent trials run in production. Leave empty to
      signal this is flattening. (defaul is None)
WORKING
  Ff : int
      The final cycle of alf to include in analysis (inclusive)
  skipE : int, optional
      In longer production runs the number of lambda samples may require
      significant amounts of memory to store and analyze. Only alchemical
      frames with index modulus skipE equal to skipE-1 are analyzed.
      (default is 1 to analyze all frames) 
  """

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
