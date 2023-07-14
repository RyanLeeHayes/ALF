#! /usr/bin/env python

def GetDCD(istep,ndupl=None,begres=None,endres=None,firstE=0,skipE=1,engine='charmm'):
  """See example at /dfs8/rhayes1_lab/rhayes1/75_MSLD/06_alfdebug/alf-20230707/testing/T4L149U_bladelib when ready to make documentation"""
  import os, os.path
  import alf
  import MDAnalysis as mda

  alf_info=alf.initialize_alf_info(engine)

  if ndupl==None:
    production=False
    ndupl=1
  else:
    production=True

  nreps=alf_info['nreps']
  name=alf_info['name']

  fnmpsf="prep/minimized.psf"

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

    DDIR='run'+str(istep)
    if production:
      DDIR=DDIR+alphabet[idupl]
    if not os.path.exists(DDIR+'/dcdcat'):
      os.mkdir(DDIR+'/dcdcat')

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
        fnmout=(DDIR+'/dcdcat/'+name+'_prod%d-%d_rep%d.%s' % (begres,endres,i,fmt))
      else:
        fnmsin.append(DDIR+'/dcd/'+name+'_flat.'+fmt+reptag)
        fnmout=(DDIR+'/dcdcat/'+name+'_flat_rep%d.%s' % (i,fmt))

      trajin=mda.Universe(fnmpsf,fnmsin,topology_format='psf',format=fmt)
      # trajin.write(fnmout,frames='all')
      system=trajin.select_atoms('all')
      with mda.Writer(fnmout, system.n_atoms) as w:
        for ts in trajin.trajectory[firstE::skipE]:
          print("Saving frame at time %f" % (ts.time,))
          w.write(system)
