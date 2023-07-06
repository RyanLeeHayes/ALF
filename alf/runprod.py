
def runprod(step,a,itt0,itt,nsteps=500000,engine='charmm'):
  import os, sys, shutil, traceback, time, subprocess, random
  import numpy as np
  import alf

  alf_info=alf.initialize_alf_info(engine)

  home_dir=os.getcwd()

  run_dir=('run'+str(step)+a)
  nlambdasteps=nsteps//10 # assume nsavl=10

  for i in range(itt0,-1,-1):
    fnm=(run_dir+'/res/%s_prod%d.lmd' % (alf_info['name'],i))
    if i==0 or alf.GetSteps(alf_info,fnm)==nlambdasteps:
      ibeg=i+1
      break
    else:
      print("Error: Run %d incomplete. Going back one step" % i)

  if ibeg==1:
    # Can't just remove run_dir and start over, because otherwise a scripting error might erase 100's of ns of sampling
    if not os.path.exists(run_dir):
      os.mkdir(run_dir)
    if not os.path.exists(run_dir+'/dcd'):
      os.mkdir(run_dir+'/dcd')
    if not os.path.exists(run_dir+'/res'):
      os.mkdir(run_dir+'/res')
    if not os.path.exists(run_dir+'/failed'):
      os.mkdir(run_dir+'/failed')
    if os.path.exists(run_dir+'/variablesprod.inp'):
      os.remove(run_dir+'/variablesprod.inp')
    shutil.copy('variables%d.inp' % step,run_dir+'/variablesprod.inp')
    if not os.path.exists(run_dir+'/prep'):
      os.symlink('../prep',run_dir+'/prep') # ../prep is relative to final path, not current directory
  os.chdir(run_dir)

  for i in range(ibeg,itt+1):
    fnm=('res/%s_prod%d.lmd' % (alf_info['name'],i))
    while alf.GetSteps(alf_info,fnm)!=nlambdasteps:
      if os.path.exists('output_%d' % i):
        os.rename('output_%d' % i,'failed/output_%d' % i)
      if os.path.exists('error_%d' % i):
        os.rename('error_%d' % i,'failed/error_%d' % i)

      try:
        fpout=open('output_%d' % i,'w')
        fperr=open('error_%d' % i,'w')
        if not os.path.exists('../msld_prod.inp'):
          print("Error: msld_prod.inp does not exist.")
        if engine in ['charmm']:
          subprocess.call(['mpirun','-np',str(alf_info['nreps']),'-x','OMP_NUM_THREADS=4','--bind-to','none','--bynode',alf_info['enginepath'],'nsteps=%d' % nsteps,'nsavc=10000','seed=%d' % random.getrandbits(16),'itt=%d' % i,'-i','../msld_prod.inp'],stdout=fpout,stderr=fperr)
        elif engine in ['bladelib']:
          subprocess.call(['mpirun','-np',str(alf_info['nreps']),'-x','OMP_NUM_THREADS=1','--bind-to','none','--bynode',alf_info['enginepath'],'nsteps=%d' % nsteps,'nsavc=10000','seed=%d' % random.getrandbits(16),'itt=%d' % i,'-i','../msld_prod.inp'],stdout=fpout,stderr=fperr)
        elif engine in ['blade']:
          fpin=open('arguments.inp','w')
          fpin.write("variables set nsteps %d\nvariables set itt %d" % (nsteps,i))
          fpin.close()
          subprocess.call(['mpirun','-np',str(alf_info['nreps']),'-x','OMP_NUM_THREADS=1','--bind-to','none','--bynode',alf_info['enginepath'],'../msld_prod.inp'],stdout=fpout,stderr=fperr)
        else:
          print("Error: unsupported engine type %s" % alf_info['engine'])
          quit()
      except Exception:
        sys.stdout.flush()
        sys.stderr.flush()
        traceback.print_exc()
        sys.stdout.flush()
        sys.stderr.flush()
