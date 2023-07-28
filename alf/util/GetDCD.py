#! /usr/bin/env python

def GetDCD(istep,ndupl=None,begres=None,endres=None,firstE=0,skipE=1,engine='charmm'):
  """
  Concatenate production trajectory files

  Run from main directory. Concatenates dcd (or xtc) trajectories from
  run[istep][a]/dcd/[name]_prod[res].dcd into a single file
  run[istep][a]/dcdcat/[name]_prod[begres]-[endres]_rep[irep].dcd, where
  [istep] is the cycle of alf, [a] is the letter index from the [ndupl]
  indices, [name] is the name of the system in alf_info['name'], [res]
  ranges from [begres+1] to [endres], and [irep] is the replica index.

  Parameters
  ----------
  istep : int
      The cycle of alf for which to concatenate trajectory files
  ndupl : int, optional
      The number of independent trials of production. If missing, this
      will be treated as a flattening run, but there's only one trajectory
      in flattening, so there's no need to concatenate it. Each trial will
      be concatenated into a separate trajectory. (default is None)
  begres : int, optional
      The number of chunks of production to discard for equilibration
      before beginning to write trajectory file. (default is None)
  endres : int, optional
      The number of chunks of production run by alf.runprod. (default is
      None)
  firstE : int, optional
      The first frame of each chunk to include. (default is 0)
  skipE : int, optional
      The frequency to include frames in trajectory, only frames with
      modulus skipE equal to firstE will be included. (default is 1)
  engine : str, optional
      The engine used to run molecular dynamics. Used for determining the
      correct file extension and format. (defaul is 'charmm')
  """

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
