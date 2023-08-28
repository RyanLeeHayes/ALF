
def runflat(ni,nf,esteps,nsteps,engine='charmm',G_imp=None,ntersite=[0,0]):
  """
  Run several cycles of short simulations followed by bias optimization

  runflat runs many short cycles of simulation followed by flattening. It
  assumes the existence of a directory prep, described in the README.md
  file, and that initialize has been run first to create the
  variables1.inp and analysis0 directories (along with nbhift) equired to
  start the first run. Each cycle [i] a run directory is created called
  run[i], the MD engine is launched with the msld_flat.inp script (you
  may write your own, otherwise a default script will be copied from
  alf/default_scripts), and the alchemical trajectory is saved to
  run[i]/res/[name]_flat.lmd, where [name] is the name of the system in
  alf_info (if you write your own script, make sure the trajectory is
  saved here, similar naming conventions apply for the spatial trajectory
  if using charge change corrections). Next the simulation is analyzed
  in analysis[i], and new biases are chosen using the GetLambdas,
  GetEnergy, RunWham, GetFreeEnergy5, and SetVars routines. These biases
  are saved in analysis[i]/b_sum.dat, analysis[i]/c_sum.dat,
  analysis[i]/x_sum.dat, and analysis[i]/s_sum.dat, and are exported to a
  format the MD engine can read in variables[i+1].inp. Then the [i+1]
  cycle of dynamics is run, followed by its analysis, and so on.

  runflat is designed to be fault-proof against MD crashes, slurm
  cancellation, and other sources of crashes. Consequently, each cycle is
  run repeatedly until it succeeds, and if a cycle has previously
  succeeded, it is not run again. Whether cycle [i] succeeded is
  determined by whether the file analysis[i]/b_sum.dat exists or not. If a
  run is determined to have failed, the previous copy of run[i] is moved
  to run[i]_failed, and the previous copy of run[i]_failed (if any) is
  deleted. Examining run[i]_failed/output and run[i]_failed/error can give
  insight into the causes of the failure. Thus if runflat is run from
  cycle 1 to 100, but is killed in the middle of cycle 43, running runflat
  again with the same arguments will not run cycles 1-42 again, but will
  move the partially completed run43 to run43_failed, and start again at
  cycle 43. If you wish to rerun cycles 1-42, you will need to delete
  b_sum.dat from all the analysis directories, or better yet remove
  analysis[1-43] and run[1-43].

  Parameters
  ----------
  ni : int
      The starting cycle (1 if no previous runs have been performed)
  nf : int
      The final cycle
  esteps : int
      The number of time steps of molecular dynamics simulation to discard
      for equilibration each cycle of flattening
  nsteps : int
      The number of time steps of molecular dynamics simulation to use for
      sampling chemical space each cycle of flattening
  engine : str, optional
      The molecular dynamics engine string, see help(alf) for allowed
      values. (default is 'charmm')
  G_imp : str, optional
      To use a G_imp directory other than the default one in alf/G_imp,
      provide a path. (default is None)
  ntersite : list of two ints, optional
      Flags for whether to use intersite coupling (first element) and
      intersite profiles (second element) in flattening. (default is [0,0]
      for no coupling. If multiple sites are simulated and coupling is
      expected, [0,1] is recommended)
  """

  import os, sys, shutil, traceback, time, subprocess, random
  import numpy as np
  import alf

  alf_info=alf.initialize_alf_info(engine)

  home_dir=os.getcwd()

  iri=max(ni-5,1)

  fex='inp'
  if engine in ['pycharmm']:
    fex='py'

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
        shutil.copy('variables%d.%s' % (i,fex),'run%d/variablesflat.%s' % (i,fex))
        os.symlink('../prep','run%d/prep' % i) # ../prep is relative to final path, not current directory
        os.chdir('run%d' % i)

        # Run the simulation
        print('run%d started' % i)
        fpout=open('output','w')
        fperr=open('error','w')
        if not os.path.exists('../msld_flat.'+fex):
          print("Error: msld_flat.%s does not exist." % fex)
        if engine in ['charmm']:
          subprocess.call(['mpirun','-np',str(alf_info['nreps']),'-x','OMP_NUM_THREADS=4','--bind-to','none','--bynode',alf_info['enginepath'],'esteps=%d' % esteps,'nsteps=%d' % nsteps,'seed=%d' % random.getrandbits(16),'-i','../msld_flat.inp'],stdout=fpout,stderr=fperr)
        elif engine in ['bladelib']:
          subprocess.call(['mpirun','-np',str(alf_info['nreps']),'-x','OMP_NUM_THREADS=1','--bind-to','none','--bynode',alf_info['enginepath'],'esteps=%d' % esteps,'nsteps=%d' % nsteps,'seed=%d' % random.getrandbits(16),'-i','../msld_flat.inp'],stdout=fpout,stderr=fperr)
        elif engine in ['blade']:
          fpin=open('arguments.inp','w')
          fpin.write("variables set esteps %d\nvariables set nsteps %d" % (esteps,nsteps))
          fpin.close()
          subprocess.call(['mpirun','-np',str(alf_info['nreps']),'-x','OMP_NUM_THREADS=1','--bind-to','none','--bynode',alf_info['enginepath'],'../msld_flat.inp'],stdout=fpout,stderr=fperr)
        elif engine in ['pycharmm']:
          fpin=open('arguments.py','w')
          fpin.write("esteps=%d\nnsteps=%d" % (esteps,nsteps))
          fpin.close()
          subprocess.call(['python','../msld_flat.py'],stdout=fpout,stderr=fperr)
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
        np.savetxt('analysis%d/nsubs' % i,np.array(alf_info['nsubs']).reshape((1,-1)),fmt=' %d')

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
        if 'lmalf' in alf_info:
          print("Warning, LMALF is not polished yet")
          # subprocess.call([shutil.which('python'),'-c','import alf; alf.RunWham(%d,0,0)' % (N*alf_info['nreps'])],stdout=fpout,stderr=fperr) # `cat ../ntersiteflat`
          subprocess.call([os.path.dirname(__file__)+'/lmalf/RunWham.sh',str(N*alf_info['nreps']),'0','0'],stdout=fpout,stderr=fperr)
          lambdafiles=[]
          for fileindex in range(N*alf_info['nreps']):
            lambdafiles.append('Lambda/Lambda%d.dat' % (fileindex+1))
          fpcat=open('Lambda/Lambda.dat','w')
          subprocess.call(['cat']+lambdafiles,stdout=fpcat)
          fpcat.close()
          subprocess.call([os.path.dirname(__file__)+'/lmalf/RunLMALF.sh','0','0','Lambda/Lambda.dat','weight.dat','OUT.dat'],stdout=fpout,stderr=fperr)
          print('Warning: coupling flags ignored')
          # alf.GetFreeEnergy5(alf_info,0,0) # `cat ../ntersiteflat`
          alf.GetFreeEnergyLM(alf_info,ntersite[0],ntersite[1])
        else:
          subprocess.call([shutil.which('python'),'-c','import alf; alf.RunWham(%d,%f,%d,%d)' % (N*alf_info['nreps'],alf_info['temp'],ntersite[0],ntersite[1])],stdout=fpout,stderr=fperr)
          alf.GetFreeEnergy5(alf_info,ntersite[0],ntersite[1])

        alf.SetVars(alf_info,i+1)
        fp=open('../variables%d.%s' % (i+1,fex),'a')
        if engine in ['charmm','bladelib']:
          fp.write('set restartfile = "../run%d/res/%s_flat.res\ntrim restartfile from 2\n' % (ir,alf_info['name']))
        elif engine in ['blade']:
          fp.write('variables set restartfile ../run%d/res/%s_flat.res' % (ir,alf_info['name']))
        elif engine in ['pycharmm']:
          fp.write('restartfile=\'../run%d/res/%s_flat.res\'\n' % (ir,alf_info['name']))
        fp.close()
      except Exception:
        sys.stdout.flush()
        sys.stderr.flush()
        traceback.print_exc()
        sys.stdout.flush()
        sys.stderr.flush()
      os.chdir(home_dir)
