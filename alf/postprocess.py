
def postprocess(i,eqS,S,N,skipE=1,boolflat=True,engine='charmm',G_imp=None,ntersite=[0,0]):
  """
  analyze a longer production run from runflat for improved biases and dG

  postprocess analyzes productions simulations run by runprod. As
  described in the documentation for runprod, runprod is run on a cycle of
  flattening [i], on several independent trials with letter indices [a],
  in several chunks [itt] typically of 1 ns each, with alchemical
  trajectories saved in run[i][a]/res/[name]_prod[itt].lmd, where [name]
  is the name of the system in alf_info. postprocess analyzes these
  alchemical trajectories, discarding some (typically the first quarter)
  for equilibration. The routine estimates new biases in a similar fashion
  to runflat, using the GetLambdas, GetEnergy, RunWham, GetFreeEnergy5,
  and SetVars routines, but only analyzing the independent trials from
  this cycle of flattening/production. The routine also calls GetVariance
  to estimate the free energy change upon chemical perturbation in this
  ensemble using the histogram based estimator. The Potts model estimator
  can be used by finding the appropriate routines in the alf README.md and
  python module documentation. Uncertainties in free energy estimates are
  determined by bootstrapping samples from the independent trials. Free
  energies and uncertainties are saved in a file called
  analysis[i]/Result.txt. The discrete solvent correction for charge
  changing mutations is computed automatically by the GetVolumes routine
  if a key named 'q' exists in alf_info that contains the charge of each
  alchemical group. One should take the difference of these free energies
  from the free energies in another ensemble to obtain the relative free
  energy of interest.

  Parameters
  ----------
  i : int
      The current cycle of ALF (the 'step' with which runflat was run)
  eqS : int
      The number of chunks from runprod to skip for equilibration
  S : int
      The total number of chunks from runprod (chunks are typically 1 ns)
  N : int
      The number of independent trials run by runprod
  skipE : int, optional
      In longer production runs the number of lambda samples may require
      significant amounts of memory to store and analyze. Only alchemical
      frames with index modulus skipE equal to skipE-1 are analyzed.
      (default is 1 to analyze all frames)
  boolflat : bool, optional
      Set to false to skip bias optimization and only estimate free
      energies with the histogram based estimator. (default is True)
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
    np.savetxt('analysis%d/nsubs' % i,np.array(alf_info['nsubs']).reshape((1,-1)),fmt=' %d')
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
      subprocess.call([os.path.dirname(__file__)+'/lmalf/RunLMALF.sh','0','0','Lambda/Lambda.dat','weight.dat','OUT.dat'])
      print('Warning: coupling flags ignored')
      # alf.GetFreeEnergy5(alf_info,0,0) # `cat ../ntersiteflat`
      alf.GetFreeEnergyLM(alf_info,ntersite[0],ntersite[1])
    else:
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
