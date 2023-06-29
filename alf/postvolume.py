
def postvolume(i,eqS,S,N,skipE=1,engine='charmm'):
  import os, sys, shutil, traceback, time, subprocess, random
  import numpy as np
  import alf

  alf_info=alf.initialize_alf_info(engine)

  try:
    os.chdir('analysis%d' % i)
    alf.GetVolumes(alf_info,i,N,eqS,S)
  except Exception:
    sys.stdout.flush()
    sys.stderr.flush()
    traceback.print_exc()
    sys.stdout.flush()
    sys.stderr.flush()
