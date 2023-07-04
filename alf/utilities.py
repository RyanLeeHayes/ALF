
def initialize_alf_info(engine='charmm'):
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

  engines=['charmm','bladelib','blade']
  if not engine in engines:
    print('Error: unsupported engine')
    print('Supported engine types are:')
    print(engines)
    quit()
  alf_info['engine']=engine

  if not os.path.exists('msld_flat.inp'):
    shutil.copy(os.path.dirname(__file__)+'/default_scripts/%s_flat.inp' % engine,'msld_flat.inp')
    print("Note: copied default script for msld_flat.inp\nIf you modify or replace this file, it will not be overwritten")
  if not os.path.exists('msld_prod.inp'):
    shutil.copy(os.path.dirname(__file__)+'/default_scripts/%s_prod.inp' % engine,'msld_prod.inp')
    print("Note: copied default script for msld_prod.inp\nIf you modify or replace this file, it will not be overwritten")

  return alf_info
