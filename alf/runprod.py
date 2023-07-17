
def runprod(step,a,itt0,itt,nsteps=500000,engine='charmm'):
  """
  run a longer production run of lambda dynamics

  runprod runs production lambda dynamics simulations. The routine assumes
  the existence of a prep directory whose format is specified in
  README.md. It is run as if it were another cycle [i] of flattening,
  where [i] is given by the input parameter [step], so analysis[i-1] and
  variables[i].inp are also assumed to exist. Dynamics is run using the
  MD engine input script msld_prod.inp which can be written by the user,
  or will be copied from the default scripts in alf/default_scripts if
  none exists. This input script should take in arguments nsteps for the
  number of MD steps to run per chunk, itt for the index of the current
  chunk (starting at 1) and write lambda trajectories to
  run[i][a]/res/[name]_prod[itt].lmd, where [i] is the cycle index, [a] is
  the duplicate independent trial index, [name] is the name of the system
  alf_info, and [itt] is the chunk index. In order to make estimates of
  statistical reproducibility, multiple independent trials (often 5) are
  run, and each independent trial is given its own letter index [a], with
  the first trial getting 'a', the second trial getting 'b', and so on.
  runprod breaks the simulation up into chunks, by default 1 ns long, and
  each call of runprod is responsible for completing some number of these
  chunks, controlled by the [itt0] and [itt] parameters.

  See README.md for best practices on how these production simulations
  should be run. Specifically, a production simulation should not be more
  than a factor of 4-5 times longer than the previous flattening
  simulation, otherwise new regions of conformational and chemical space
  might open up that shift the optimal values for the biases. The
  postprocess routine may be used to analyze the output of runprod both to
  determine new optimal values of the biases for subsequent cycles of
  production, and to create a Result.txt file containing estimates of dG
  for the chemical processes in the sampled ensemble.

  runprod is intended to run production lambda dynamics simulations in a
  fault-proof manner that is robust against MD crashes and slurm job
  cancellation. Accordingly, each chunk it attempted, starting with
  [itt0]+1, until it succeeds, annd is not rerun if it has already
  succeeded. Success is determined by whether the lambda trajectory
  run[i][a]/res/[name]_prod[itt].lmd exists and contains the correct
  number of steps [nsteps]/10, because msld_prod.inp is assumed to print
  out a lambda frame every 10 time steps. If itt0 or previous chunks also
  failed, a warning will be printed, and those chunks will be rerun. If a
  chunk fails, its output run[i][a]/output_[itt] and error
  run[i][a]/error_[itt] will be moved to the run[i][a]/failed/ directory,
  and can be useful for diagnosing crashes. Each chunk reads a restart
  file from the previous chunk. runprod is intended to be submitted as
  a slurm job array, with a maximum number of simultaneous array elements
  of 1, so that each element of the array completes some number of chunks,
  and a potentially very long simulation is broken up into shorter jobs.
  If something goes wrong, the entire job array may be resubmitted, an
  successfully completed chunks will be skipped and not rerun. If you wish
  to rerun starting from the first chunk, all the previous
  run[i][a]/res/[name]_prod[itt].lmd files for all values of [itt] must be
  deleted, or better yet, delete the entire run[i][a] directory.

  Parameters
  ----------
  step : int
      The current cycle of ALF (1 more than the previous cycle of
      flattening)
  a : char
      The duplicate independent trial character index, starting from 'a'
      and continuing without gaps
  itt0 : int
      The index of the first chunk before the chunks this instance of
      runprod is responsible for. This should be 0 for the first call of
      runprod. Incomplete chunks before and including itt0 will be rerun
      if incomplete, but will give an error message.
  itt : int
      The index of the final chunk this instance of runprod is responsible
      for
  nsteps : int, optional
      The number of time steps of molecular dynamics simulation in each
      chunk of the simulation. By default the simulation uses 1 ns chunks.
      (default is 500000)
  engine : str, optional
      The molecular dynamics engine string, see help(alf) for allowed
      values. (default is 'charmm')
  """

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
