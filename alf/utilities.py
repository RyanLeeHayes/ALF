
def initialize_alf_info(engine='charmm'):
  """
  Wrapper function to read prep/alf_info.py

  Reads prep/alf_info.py to initialize alf_info dictionary. Ensures
  necessary fields are present and formatted correctly. See alf README.md
  for details on a correctly formatted prep/alf_info.py file.

  Parameters
  ----------
  engine : str, optional
      An optional string for the molecular dynamics engine to be used. See
      alf README.md for supported engine strings. (default is 'charmm')

  Returns
  -------
  alf_info : dict
      The alf_info dictionary with parameters most alf routines need to
      run
  """

  import os, sys, shutil, traceback, time, subprocess, random
  import numpy as np
  import alf

  if not os.path.exists('prep/alf_info.py'):
    print("Error, prep/alf_info.py is not defined")
    quit()
  try:
    from prep.alf_info import alf_info
  except Exception:
    print("Error in prep/alf_info.py")
    traceback.print_exc()
    quit()

  # Error checking:
  if not 'name' in alf_info:
    print("Error: need name key in prep/alf_info.py. name = a string that points to the system setup file prep/name.inp")
    quit()
  if not 'nsubs' in alf_info:
    print("Error: need nsubs key in prep/alf_info.py. nsubs = a vector of integers for the number of substituents at each site")
    quit()
  if ( not 'nblocks' in alf_info ) or ( not np.sum(alf_info['nsubs']) == alf_info['nblocks'] ):
    print("Error: need nblocks key in prep/alf_info.py. nblocks = the sum of nsubs")
    quit()
  if not 'ncentral' in alf_info:
    print("Error: need ncentral key in prep/alf_info.py. ncentral is for replica exchange. Use 0 if you are not using replica exchange")
    quit()
  if not 'nreps' in alf_info:
    print("Error: need nreps key in prep/alf_info.py. nreps is for replica exchange. Use 1 if you are not using replica exchange")
    quit()
  if not 'nnodes' in alf_info:
    print("Error: need nnodes key in prep/alf_info.py. nnodes = number of nodes for parallelization, 1 is recommended")
    quit()
  if not 'enginepath' in alf_info:
    print("Error: need enginepath key in prep/alf_info.py. enginepath = string of path to molecular dynamics executable")
    quit()
  if not 'temp' in alf_info:
    print("Error: need temp key in prep/alf_info.py. temp = system temperature in Kelvin")
    quit()

  # Convert to numpy arrays
  alf_info['nsubs']=np.array(alf_info['nsubs'])
  if 'q' in alf_info:
    alf_info['q']=np.array(alf_info['q'])

  engines=['charmm','bladelib','blade','pycharmm']
  if not engine in engines:
    print('Error: unsupported engine')
    print('Supported engine types are:')
    print(engines)
    quit()
  alf_info['engine']=engine

  # Get file extension to copy default flattening and production scripts
  pyengines=['pycharmm']
  csengines=['charmm','bladelib']
  if engine in pyengines: # python script
    fex='py'
  elif engine in csengines: # charmm script
    fex='inp'
  else: # blade script
    fex='inp'

  if not os.path.exists('msld_flat.'+fex):
    shutil.copy(os.path.dirname(__file__)+'/default_scripts/%s_flat.%s' % (engine,fex),'msld_flat.'+fex)
    print("Note: copied default script for msld_flat.%s\nIf you modify or replace this file, it will not be overwritten" % (fex,))
  if not os.path.exists('msld_prod.'+fex):
    shutil.copy(os.path.dirname(__file__)+'/default_scripts/%s_prod.%s' % (engine,fex),'msld_prod.'+fex)
    print("Note: copied default script for msld_prod.%s\nIf you modify or replace this file, it will not be overwritten" % (fex,))

  return alf_info
