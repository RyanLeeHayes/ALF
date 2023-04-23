#! /usr/bin/env python

def GetVolumes(alf_info,istep,ndupl=None,begres=None,endres=None):
  import sys, os, os.path
  import numpy as np
  # from subprocess import call
  from alf.GetVolume import GetVolume

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
  ncentral=alf_info['ncentral']
  name=alf_info['name']

  if not os.path.isfile('../prep/q'):
    print("No charge file prep/q")
    return
  else:
    q=np.loadtxt('../prep/q')
    ibuff=0
    chargeChange=False
    for i in range(len(nsubs)):
      if not np.all(q[ibuff:(ibuff+nsubs[i])]==q[ibuff]):
        chargeChange=True
      ibuff=ibuff+nsubs[i]
    if not chargeChange:
      print("No charge change")
      return

  if not os.path.isdir('data'):
    os.mkdir('data')
  DIR="data"

  # PSF="../prep/"+name+".psf"
  PSF="../prep/minimized.psf"

  if alf_info['engine'] in ['charmm','bladelib']:
    fmt="dcd"
  elif alf_info['engine'] in ['blade']:
    fmt="xtc"
  else:
    print("Error: unsupported engine type %s" % alf_info['engine'])
    quit()

  # ----------------------------------------------------------------------------

  alphabet='abcdefghijklmnopqrstuvwxyz'

  for idupl in range(0,ndupl):

    DDIR='../run'+str(istep)
    if production:
      DDIR=DDIR+alphabet[idupl]

    datasplit=[]

    for i in range(0,nreps):
      fnmsin=[]
      if nreps>1:
        reptag="_"+str(i)
      else:
        reptag=""
      if production:
        for j in range(begres,endres):
          fnmsin.append(DDIR+'/dcd/'+name+'_prod'+str(j+1)+'.'+fmt+reptag)
      else:
        fnmsin.append(DDIR+'/dcd/'+name+'_flat.'+fmt+reptag)
      fnmout=DIR+("/Volume.%d.%d.dat" % (idupl,i))
      GetVolume(alf_info,fnmout,fnmsin,PSF)
      datasplit.append(np.loadtxt(fnmout))

    nsteps=len(datasplit[-1])

  # ----------------------------------------------------------------------------

    for rep in range(0,nreps):
      fnm=DIR+("/Volume.%d.%d.dat" % (idupl,rep))
      V=np.mean(np.loadtxt(fnm))
      fnm=DIR+("/Volume.%d.%d.dat_NH2O" % (idupl,rep))
      N=np.loadtxt(fnm)
      fnm=DIR+("/Density.%d.%d.dat" % (idupl,rep))
      np.savetxt(fnm,np.reshape(N/V,(1,)),fmt="%12.8f")

  D=np.zeros((ndupl,))
  for idupl in range(0,ndupl):
    D[idupl]=np.loadtxt(DIR+("/Density.%d.%d.dat" % (idupl,ncentral)))
  Davg=1.0/np.mean(1.0/D)
  np.savetxt(DIR+"/Density.dat",np.reshape(Davg,(1,)),fmt="%12.8f")

  gamma_S=0.76414
  kelectric=332.0716
  b_corr=-(2.0*np.pi/3.0)*kelectric*gamma_S*q*Davg
  np.savetxt("b_corr.dat",[b_corr]) # Square brackets save as row vector
