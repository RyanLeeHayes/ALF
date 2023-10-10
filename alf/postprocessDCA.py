
def SetupDCA(i,NF,FREQ,engine='charmm'):
  """
  Sets up directories for Pott model estimator analysis

  Assumes postprocess has been run first, Uses same eqS and S chunks as
  were used by postprocess.

  This is the first routine in the Potts model estimator, call FilterDCA
  next.

  Parameters
  ----------
  i : int
      The current cycle of ALF (the 'step' with which runflat was run)
  NF : int
      The number of independent trials run by runprod
  FREQ : int
      In longer production runs the number of lambda samples may require
      significant amounts of memory to store and analyze. Only alchemical
      frames with index modulus FREQ equal to FREQ-1 are analyzed. 1
      analyzes all frames.
  engine : str, optional
      The molecular dynamics engine string, see help(alf) for allowed
      values. (default is 'charmm')
  """

  import os, sys, shutil, traceback, time, subprocess, random
  import numpy as np
  import alf

  alf_info=alf.initialize_alf_info(engine)

  home_dir=os.getcwd()

  # Prep the analysis directories
  if not os.path.exists('dca%d' % i):
    os.mkdir('dca%d' % i)
  shutil.copy('analysis%d/b_sum.dat' % (i-1),'dca%d/b_prev.dat' % i)
  shutil.copy('analysis%d/c_sum.dat' % (i-1),'dca%d/c_prev.dat' % i)
  shutil.copy('analysis%d/x_sum.dat' % (i-1),'dca%d/x_prev.dat' % i)
  shutil.copy('analysis%d/s_sum.dat' % (i-1),'dca%d/s_prev.dat' % i)
  np.savetxt('dca%d/nsubs' % i,np.array(alf_info['nsubs']).reshape((1,-1)),fmt=' %d')
  if not os.path.exists('analysis%d' % i):
    print("Error: analysis%d does not exist, run alf.postprocess first" % i)
    quit()
  if os.path.isfile('analysis%d/b_corr.dat' % i):
    shutil.copy('analysis%d/b_corr.dat' % i,'dca%d/b_corr.dat' % i)
  if not os.path.exists('dca%d/data' % i):
    os.mkdir('dca%d/data' % i)
  time.sleep(15)

def FilterDCA(i,iNF,NF,FREQ,engine='charmm'):
  """
  Filters continuous lambda trajectory into discrete lambda states

  At each time step (row of output) and each site (column of output) the
  index block with a lambda value greater than 0.99 is printed, starting
  from an index of 1 for the first site. If no lambda value is above this
  threshhold, a 0 is printed for the alchemical intermediate state. This
  function should be called NF times, once with each of the iNF values
  from 0 to NF-1.

  This is the second routine in the Potts model estimator. SetupDCA should
  be called before it, and MomentDCA should be called after it.

  Parameters
  ----------
  i : int
      The current cycle of ALF (the 'step' with which runflat was run)
  iNF : int
      The independent trial this call of FilterDCA is responsible for
      analyzing
  NF : int
      The number of independent trials run by runprod
  FREQ : int
      In longer production runs the number of lambda samples may require
      significant amounts of memory to store and analyze. Only alchemical
      frames with index modulus FREQ equal to FREQ-1 are analyzed. 1
      analyzes all frames.
  engine : str, optional
      The molecular dynamics engine string, see help(alf) for allowed
      values. (default is 'charmm')
  """

  import os, sys, shutil, traceback, time, subprocess, random
  import numpy as np
  import alf

  alf_info=alf.initialize_alf_info(engine)

  home_dir=os.getcwd()

  os.chdir('dca%d' % i)

  exe=os.path.dirname(os.path.abspath(__file__))+'/dca/Filter'
  fnmin=('../analysis%d/data/Lambda.%d.%d.dat' % (i,iNF,alf_info['ncentral']))
  fnmout=('data/Filter.%d.dat' % iNF)
  subprocess.call(['mpirun','-n','1','--map-by','node','--bind-to','none','-x','OMP_NUM_THREADS=1',exe,str(FREQ),fnmin,fnmout])
  print("WARNING: direct standard error somewhere")
  time.sleep(15)

def MomentDCA(i,iNF,NF,FREQ,engine='charmm'):
  """
  Computes moments of the discrete lambda states classified by FilterDCA

  Computes the 1D moments (the probability of observing each discrete
  lambda state) and the 2D moments (the probability of observing each pair
  of discrete lambda states, and saves them to dca[i]/m1.[iNF].obs.dat and
  dca[i]/m2.[iNF].obs.dat respectively. These moments are not strictly
  required for PLM Potts model analysis, only for LM Potts model analysis,
  but they are read by the subsequent function, and should be called to
  prevent errors. This function should be called NF times, once with each
  of the iNF values from 0 to NF-1.

  This is the third routine in the Potts model estimator. FilterDCA should
  be called before it and BSMomentDCA should be called after it.

  Parameters
  ----------
  i : int
      The current cycle of ALF (the 'step' with which runflat was run)
  iNF : int
      The independent trial this call of MomentDCA is responsible for
      analyzing
  NF : int
      The number of independent trials run by runprod
  FREQ : int
      In longer production runs the number of lambda samples may require
      significant amounts of memory to store and analyze. Only alchemical
      frames with index modulus FREQ equal to FREQ-1 are analyzed. 1
      analyzes all frames.
  engine : str, optional
      The molecular dynamics engine string, see help(alf) for allowed
      values. (default is 'charmm')
  """

  import os, sys, shutil, traceback, time, subprocess, random
  import numpy as np
  import alf

  alf_info=alf.initialize_alf_info(engine)

  home_dir=os.getcwd()

  os.chdir('dca%d' % i)

  exe=os.path.dirname(os.path.abspath(__file__))+'/dca/Moment'
  fnmin=('data/Filter.%d.dat' % iNF)
  fnmout1=('data/m1.%d.obs.dat' % iNF)
  fnmout2=('data/m2.%d.obs.dat' % iNF)
  subprocess.call(['mpirun','-n','1','--map-by','node','--bind-to','none','-x','OMP_NUM_THREADS=1',exe,fnmin,fnmout1,fnmout2])
  print("WARNING: direct standard error somewhere")
  time.sleep(15)

def BSMomentDCA(i,NF,FREQ,NBS=50,engine='charmm'):
  """
  Sets up bootstrapping for calculation of uncertainty in Potts model

  Bootstraps the 1D and 2D moments calculated previously, whether (LM) or
  not (PLM) they are needed. Also creates a randomized bootstrapping list
  of which independent trials to combine that is used by PLM. 50 is
  recommended for the number of bootstrapping samples. This function
  should be called once.

  This is the fourth routine in the Potts model estimator. MomentDCA
  should be called before it and either LMDCA (few sites) or PLMDCA
  (many sites) should be called after it.

  Parameters
  ----------
  i : int
      The current cycle of ALF (the 'step' with which runflat was run)
  NF : int
      The number of independent trials run by runprod
  FREQ : int
      In longer production runs the number of lambda samples may require
      significant amounts of memory to store and analyze. Only alchemical
      frames with index modulus FREQ equal to FREQ-1 are analyzed. 1
      analyzes all frames.
  NBS : int, optional
      The number of bootstrap samples to take (default is 50)
  engine : str, optional
      The molecular dynamics engine string, see help(alf) for allowed
      values. (default is 'charmm')
  """

  import os, sys, shutil, traceback, time, subprocess, random
  import numpy as np
  import alf

  alf_info=alf.initialize_alf_info(engine)

  home_dir=os.getcwd()

  os.chdir('dca%d' % i)

  alf.BootstrapMomentsDCA(alf_info,NF,'data',NBS=NBS)
  time.sleep(15)

def LMDCA(i,iBS,NF,FREQ,NBS=50,engine='charmm'):
  """
  Performs likelihood maximization for Potts model estimator

  Performs likelihood maximization for either the original data or one of
  the bootstrap samples. NBS=50 is recommended for the number of
  bootstrapping samples. This function should be called [NBS+1] times.

  This is the fifth routine in the Potts model estimator. BSMomentDCA
  should be called before it and FinishDCA should be called after it.

  Parameters
  ----------
  i : int
      The current cycle of ALF (the 'step' with which runflat was run)
  iBS : int
      The bootstrap index this call of LMDCA is responsible for. 0 to BS-1
      are responsible for the bootstrap samples, BS is responsible for the
      the original dataset
  NF : int
      The number of independent trials run by runprod
  FREQ : int
      In longer production runs the number of lambda samples may require
      significant amounts of memory to store and analyze. Only alchemical
      frames with index modulus FREQ equal to FREQ-1 are analyzed. 1
      analyzes all frames.
  NBS : int, optional
      The number of bootstrap samples to take (default is 50)
  engine : str, optional
      The molecular dynamics engine string, see help(alf) for allowed
      values. (default is 'charmm')
  """

  import os, sys, shutil, traceback, time, subprocess, random
  import numpy as np
  import alf

  alf_info=alf.initialize_alf_info(engine)

  home_dir=os.getcwd()

  os.chdir('dca%d' % i)

  exe=os.path.dirname(os.path.abspath(__file__))+'/dca/LM'
  if iBS==NBS:
    fnmout1='data/h.LM.dat'
    fnmout2='data/J.LM.dat'
    fnmin1='data/m1.obs.dat' 
    fnmin2='data/m2.obs.dat'
  else:
    fnmout1=('data/h.bs%d.LM.dat' % iBS)
    fnmout2=('data/J.bs%d.LM.dat' % iBS)
    fnmin1=('data/m1.bs%d.obs.dat' % iBS)
    fnmin2=('data/m2.bs%d.obs.dat' % iBS)
  subprocess.call(['mpirun','-n','1','--map-by','node','--bind-to','none','-x','OMP_NUM_THREADS=8',exe,fnmout1,fnmout2,fnmin1,fnmin2])
  time.sleep(15)

def PLMDCA(i,iBS,NF,FREQ,NBS=50,engine='charmm'):
  """
  Performs psuedolikelihood maximization for Potts model estimator

  Performs pseudolikelihood maximization for either the original data or
  one of the bootstrap samples. NBS=50 is recommended for the number of
  bootstrapping samples. This function should be called [BS+1] times.

  This is the fifth routine in the Potts model estimator. BSMomentDCA
  should be called before it and FinishDCA should be called after it.

  Parameters
  ----------
  i : int
      The current cycle of ALF (the 'step' with which runflat was run)
  iBS : int
      The bootstrap index this call of PLMDCA is responsible for. 0 to
      BS-1 are responsible for the bootstrap samples, BS is responsible
      for the the original dataset
  NF : int
      The number of independent trials run by runprod
  FREQ : int
      In longer production runs the number of lambda samples may require
      significant amounts of memory to store and analyze. Only alchemical
      frames with index modulus FREQ equal to FREQ-1 are analyzed. 1
      analyzes all frames.
  NBS : int, optional
      The number of bootstrap samples to take (default is 50)
  engine : str, optional
      The molecular dynamics engine string, see help(alf) for allowed
      values. (default is 'charmm')
  """

  import os, sys, shutil, traceback, time, subprocess, random
  import numpy as np
  import alf

  alf_info=alf.initialize_alf_info(engine)

  home_dir=os.getcwd()

  os.chdir('dca%d' % i)

  exe=os.path.dirname(os.path.abspath(__file__))+'/dca/PLM'
  if iBS==NBS:
    fnmout1='data/h.PLM.dat'
    fnmout2='data/J.PLM.dat'
    bsindices=np.arrange(NF,dtype='int')
  else:
    fnmout1=('data/h.bs%d.PLM.dat' % iBS)
    fnmout2=('data/J.bs%d.PLM.dat' % iBS)
    bsindices=np.loadtxt('data/bs%d.dat' % iBS,dtype='int')
  fnmsin=[]
  for bsindex in bsindices:
    fnmsin.append('data/Filter.%d.dat' % bsindex)
  print(fnmsin)
  subprocess.call(['mpirun','-n','1','--map-by','node','--bind-to','none','-x','OMP_NUM_THREADS=1',exe,fnmout1,fnmout2]+fnmsin)
  time.sleep(15)

def FinishDCA(i,NF,FREQ,NBS=50,engine='charmm'):
  """
  Gives the results of the Potts model estimator

  This routine calls GetVarianceDCA to provide a Results.txt file if there
  are sufficiently few sites. This routine also calls GetModelDCA to
  create h and J matrices from which dG for arbitrary sequences can be
  determined. These h and J matrices are stored in
  dca[i]/data/h.bias.[LM].dat and dca[i]/data/J.bias.[LM]dat, and various
  bootstrapped values are stored in dca[i]/data/h.bs[iBS].[LM].dat and
  dca[i]/data/h.bs[iBS].[LM].dat , where [iBS] is the bootstrap index, and
  [LM] is either 'LM' or 'PLM' indicating how the parameters were
  determined. Uncertainties for arbitrary sequences may be computed with
  these bootstrapped h and J samples. The free energy of a sequence S is
  the sum of the h values for each site in it, plus 0.5*J for each pair of
  sites in the sequence (because each i,j site pair is counted twice, as
  both i,j and j,i). Includes the discrete solvent charge change
  correction if it was computed correctly by postprocess. This function
  should be called once.

  This is the final routine in the Potts model estimator. LMDCA or PLMDCA
  should be called before it.

  Parameters
  ----------
  i : int
      The current cycle of ALF (the 'step' with which runflat was run)
  NF : int
      The number of independent trials run by runprod
  FREQ : int
      In longer production runs the number of lambda samples may require
      significant amounts of memory to store and analyze. Only alchemical
      frames with index modulus FREQ equal to FREQ-1 are analyzed. 1
      analyzes all frames.
  NBS : int, optional
      The number of bootstrap samples to take (default is 50)
  engine : str, optional
      The molecular dynamics engine string, see help(alf) for allowed
      values. (default is 'charmm')
  """

  import os, sys, shutil, traceback, time, subprocess, random
  import numpy as np
  import alf

  alf_info=alf.initialize_alf_info(engine)

  home_dir=os.getcwd()

  os.chdir('dca%d' % i)

  alf.GetVarianceDCA(alf_info,NF,'data',NBS=NBS)
  alf.GetModelDCA(alf_info,NF,'data',NBS=NBS)
