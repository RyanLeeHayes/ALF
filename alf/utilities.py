
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
    quit(1)
  try:
    from prep.alf_info import alf_info
  except Exception:
    print("Error in prep/alf_info.py")
    traceback.print_exc()
    quit(1)

  # Error checking:
  if not 'name' in alf_info:
    print("Error: need name key in prep/alf_info.py. name = a string that points to the system setup file prep/name.inp")
    quit(1)
  if not 'nsubs' in alf_info:
    print("Error: need nsubs key in prep/alf_info.py. nsubs = a vector of integers for the number of substituents at each site")
    quit(1)
  if ( not 'nblocks' in alf_info ) or ( not np.sum(alf_info['nsubs']) == alf_info['nblocks'] ):
    print("Error: need nblocks key in prep/alf_info.py. nblocks = the sum of nsubs")
    quit(1)
  if not 'ncentral' in alf_info:
    print("Error: need ncentral key in prep/alf_info.py. ncentral is for replica exchange. Use 0 if you are not using replica exchange")
    quit(1)
  if not 'nreps' in alf_info:
    print("Error: need nreps key in prep/alf_info.py. nreps is for replica exchange. Use 1 if you are not using replica exchange")
    quit(1)
  if not 'nnodes' in alf_info:
    print("Error: need nnodes key in prep/alf_info.py. nnodes = number of nodes for parallelization, 1 is recommended")
    quit(1)
  if not 'enginepath' in alf_info:
    print("Error: need enginepath key in prep/alf_info.py. enginepath = string of path to molecular dynamics executable")
    quit(1)
  if not 'temp' in alf_info:
    print("Error: need temp key in prep/alf_info.py. temp = system temperature in Kelvin")
    quit(1)

  # Convert to numpy arrays
  alf_info['nsubs']=np.array(alf_info['nsubs'])
  if 'q' in alf_info:
    alf_info['q']=np.array(alf_info['q'])

  # Check loss function
  losses=['linear2018','nonlinear2024']
  if not 'loss' in alf_info:
    # print("Error: need loss key in prep/alf_info.py. Options are")
    # print(losses)
    # quit(1)
    print("Default loss linear2018 added to alf_info")
    alf_info['loss']='linear2018' # default
  if not alf_info['loss'] in losses:
    print("Error: unsupported loss function %s. Options are" % (alf_info['loss'],))
    print(losses)
    quit(1)

  # Check biases
  biases=['bcxs2018','bcxstu2026']
  if not 'bias' in alf_info:
    # print("Error: need bias key in prep/alf_info.py. Options are")
    # print(biases)
    # quit(1)
    print("Default bias bcxs2018 added to alf_info")
    alf_info['bias']='bcxs2018' # default
  if not alf_info['bias'] in biases:
    print("Error: unsupported bias function %s. Options are" % (alf_info['bias'],))
    print(biases)
    quit(1)

  # Check implicit constraints
  impconses=['fnex2011','fnexdozen2024','fnpwise2026','custom']
  if not 'impcons' in alf_info:
    # print("Error: need impcons key in prep/alf_info.py. Options are")
    # print(impconses)
    # quit(1)
    print("Default implicit constraint fnex2011 added to alf_info")
    alf_info['impcons']='fnex2011' # default
  if not alf_info['impcons'] in impconses:
    print("Error: unsupported implicit constraint %s. Options are" % (alf_info['impcons'],))
    print(impconses)
    quit(1)

  # Set implicit constraint options
  if alf_info['impcons'] in ['fnex2011','fnexdozen2024']:
    if not 'impconsopt' in alf_info:
      print("Default implicit constraint options placed in alf_info")
      print("alf_info['impconsopt']={'c':5.5}")
      alf_info['impconsopt']={'c':5.5}
    elif not 'c' in alf_info['impconsopt']:
      print("Default implicit constraint options placed in alf_info")
      print("alf_info['impconsopt']['c']=5.5")
      alf_info['impconsopt']['c']=5.5
  elif alf_info['impcons']=='fnpwise2026':
    if not 'impconsopt' in alf_info:
      print("Default implicit constraint options placed in alf_info")
      print("alf_info['impconsopt']={'width':0.35,'kbeta':100}")
      alf_info['impconsopt']={'width':0.35,'kbeta':100}
    else:
      if not 'width' in alf_info['impconsopt']:
        print("Default implicit constraint width option placed in alf_info")
        print("alf_info['impconsopt']['width']=0.35")
        alf_info['impconsopt']['width']=0.35
      if not 'kbeta' in alf_info['impconsopt']:
        print("Default implicit constraint kbeta option placed in alf_info")
        print("alf_info['impconsopt']['kbeta']=100")
        alf_info['impconsopt']['kbeta']=100

  engines=['charmm','bladelib','blade','pycharmm']
  if not engine in engines:
    print('Error: unsupported engine')
    print('Supported engine types are:')
    print(engines)
    quit(1)
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
