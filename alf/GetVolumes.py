#! /usr/bin/env python

def GetVolumes(alf_info,istep,ndupl=None,begres=None,endres=None):
  """
  Calculate the discrete solvent charge change correction from box volume

  The discrete solvent correction is calculated by reading the box volume
  from trajectory files in run[istep][a]/dcd/[name]_prod[itt].dcd where
  [istep] is the cycle of alf, [a] is the letter index for the [ndupl]
  independent trials, [name] is the name of the system and [itt] is the
  production chunk ranging from [begres+1] to [endres]. The number of
  water molecules from prep/minimized.psf. This routine should be run from
  analysis[istep]. Intermediate files are written to analysis[istep]/data
  and the final correction for each substituent is written to
  analysis[iste]/b_coff.dat based on the charges in alf_info['q']. If
  alf_info['q'] is all zero or is missing, no correction is calculated.
  This correction is automatically added to free energy estimates by
  GetVariance and GetVarianceDCA. This routine passes the conventional
  names of files to the routine alf.GetVolume, which actually reads the
  box size from the trajectory.

  Parameters
  ----------
  alf_info : dict
      Dictionary of variables alf needs to run
  istep : int
      The current cycle of alf being analyzed
  ndupl : int, optional
      The number of independent trials run in production. Leave empty to
      signal this is flattening. (defaul is None)
  begres : int, optional
      The number of chunks of production to discard for equilibration.
      Leave empty for flattening. (default is None)
  endres : int, optional
      The final chunks of production to use for analysis. Leave empty for
      flattening. (default is None)
  """

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

  if not 'q' in alf_info:
    print("No charge list 'q' in alf_info - no charge changing correction will be applied")
    return
  else:
    q=alf_info['q']
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

  if alf_info['engine'] in ['charmm','bladelib','pycharmm']:
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
