
def runflat(ni,nf,esteps,nsteps,engine='charmm',G_imp=None,ntersite=[0,0]):
  import os, sys, shutil, traceback, time, subprocess, random
  import numpy as np
  import alf

  alf_info=alf.initialize_alf_info(engine)

  home_dir=os.getcwd()

  iri=max(ni-5,1)

  for i in range(ni,nf+1):
    im5=max(i-4,1)
    N=i-im5+1
    ir=random.randrange(iri,i+1)

    while not os.path.exists('analysis%d/b_sum.dat' % i):
      os.chdir(home_dir)

      if os.path.exists('run%d' % i):
        print('run%d failed' % i)
        if os.path.exists('run%d_failed' % i):
          shutil.rmtree('run%d_failed' % i)
        os.rename('run%d' % i,'run%d_failed' % i)
        time.sleep(15)
 
      try:
        # Prep the run directories
        os.mkdir('run%d' % i)
        os.mkdir('run%d/dcd' % i)
        os.mkdir('run%d/res' % i)
        shutil.copy('variables%d.inp' % i,'run%d/variablesflat.inp' % i)
        os.symlink('../prep','run%d/prep' % i) # ../prep is relative to final path, not current directory
        os.chdir('run%d' % i)

        # Run the simulation
        print('run%d started' % i)
        fpout=open('output','w')
        fperr=open('error','w')
        if not os.path.exists('../msld_flat.inp'):
          print("Error: msld_flat.inp does not exist.")
        if engine in ['charmm']:
          subprocess.call(['mpirun','-np',str(alf_info['nreps']),'-x','OMP_NUM_THREADS=4','--bind-to','none','--bynode',alf_info['enginepath'],'esteps=%d' % esteps,'nsteps=%d' % nsteps,'seed=%d' % random.getrandbits(16),'-i','../msld_flat.inp'],stdout=fpout,stderr=fperr)
        elif engine in ['bladelib']:
          subprocess.call(['mpirun','-np',str(alf_info['nreps']),'-x','OMP_NUM_THREADS=1','--bind-to','none','--bynode',alf_info['enginepath'],'esteps=%d' % esteps,'nsteps=%d' % nsteps,'seed=%d' % random.getrandbits(16),'-i','../msld_flat.inp'],stdout=fpout,stderr=fperr)
        elif engine in ['blade']:
          fpin=open('arguments.inp','w')
          fpin.write("variables set esteps %d\nvariables set nsteps %d" % (esteps,nsteps))
          fpin.close()
          subprocess.call(['mpirun','-np',str(alf_info['nreps']),'-x','OMP_NUM_THREADS=1','--bind-to','none','--bynode',alf_info['enginepath'],'../msld_flat.inp'],stdout=fpout,stderr=fperr)
        else:
          print("Error: unsupported engine type %s" % alf_info['engine'])
          quit()
        os.chdir(home_dir)

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
        alf.GetLambdas(alf_info,i)
        alf.GetEnergy(alf_info,im5,i)
        fpout=open('output','w')
        fperr=open('error','w')
        subprocess.call([shutil.which('python'),'-c','import alf; alf.RunWham(%d,%f,%d,%d)' % (N*alf_info['nreps'],alf_info['temp'],ntersite[0],ntersite[1])],stdout=fpout,stderr=fperr)
        alf.GetFreeEnergy5(alf_info,ntersite[0],ntersite[1])

        alf.SetVars(alf_info,i+1)
        fp=open('../variables%d.inp' % (i+1),'a')
        if engine in ['charmm','bladelib']:
          fp.write('set restartfile = "../run%d/res/%s_flat.res\ntrim restartfile from 2\n' % (ir,alf_info['name']))
        elif engine in ['blade']:
          fp.write('variables set restartfile ../run%d/res/%s_flat.res' % (ir,alf_info['name']))
        fp.close()
      except Exception:
        sys.stdout.flush()
        sys.stderr.flush()
        traceback.print_exc()
        sys.stdout.flush()
        sys.stderr.flush()
      os.chdir(home_dir)
