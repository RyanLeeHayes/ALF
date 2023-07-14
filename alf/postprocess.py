
def postprocess(i,eqS,S,N,skipE=1,boolflat=True,engine='charmm',G_imp=None,ntersite=[0,0]):
  import os, sys, shutil, traceback, time, subprocess, random
  import numpy as np
  import alf

  alf_info=alf.initialize_alf_info(engine)

  home_dir=os.getcwd()

  name=alf_info['name']
  nreps=alf_info['nreps']

  try:
    # Prep the analysis directories
    print('analysis%d started' % i)
    if not os.path.exists('analysis%d' % i):
      os.mkdir('analysis%d' % i)
    shutil.copy('analysis%d/b_sum.dat' % (i-1),'analysis%d/b_prev.dat' % i)
    shutil.copy('analysis%d/c_sum.dat' % (i-1),'analysis%d/c_prev.dat' % i)
    shutil.copy('analysis%d/x_sum.dat' % (i-1),'analysis%d/x_prev.dat' % i)
    shutil.copy('analysis%d/s_sum.dat' % (i-1),'analysis%d/s_prev.dat' % i)
    if not os.path.exists('analysis%d/G_imp' % i):
      if not G_imp:
        G_imp_dir=os.path.dirname(os.path.abspath(__file__))+'/G_imp'
      else:
        G_imp_dir=G_imp
      os.symlink(G_imp_dir,'analysis%d/G_imp' % i)
    os.chdir('analysis%d' % i)

    # Run the analysis
    alf.GetLambdas(alf_info,i,N,eqS,S)
    alf.GetEnergy(alf_info,i,i,skipE)
    fpout=open('output','w')
    fperr=open('error','w')
    subprocess.call([shutil.which('python'),'-c','import alf; alf.RunWham(%d,%f,%d,%d)' % (N*alf_info['nreps'],alf_info['temp'],ntersite[0],ntersite[1])],stdout=fpout,stderr=fperr)
    alf.GetFreeEnergy5(alf_info,ntersite[0],ntersite[1])

    alf.SetVars(alf_info,i+1)
    alf.GetVolumes(alf_info,i,N,eqS,S)
    alf.GetVariance(alf_info,N)
  except Exception:
    sys.stdout.flush()
    sys.stderr.flush()
    traceback.print_exc()
    sys.stdout.flush()
    sys.stderr.flush()
