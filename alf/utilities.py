
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

  engines=['charmm','bladelib','blade']
  if not engine in engines:
    print('Error: unsupported engine')
    print('Supported engine types are:')
    print(engines)
    quit()
  alf_info['engine']=engine

  return alf_info
