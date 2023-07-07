
def SetupDCA(i,NF,FREQ,engine='charmm'):
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
  print("WARNING: Check whether analysis%d exists")
  if os.path.isfile('analysis%d/b_corr.dat' % i):
    shutil.copy('analysis%d/b_corr.dat' % i,'dca%d/b_corr.dat' % i)
  if not os.path.exists('dca%d/data' % i):
    os.mkdir('dca%d/data' % i)
  time.sleep(15)

def FilterDCA(i,iNF,NF,FREQ,engine='charmm'):
  import os, sys, shutil, traceback, time, subprocess, random
  import numpy as np
  import alf

  alf_info=alf.initialize_alf_info(engine)

  home_dir=os.getcwd()

  os.chdir('dca%d' % i)

  exe=os.path.dirname(os.path.abspath(__file__))+'/dca/Filter'
  fnmin=('../analysis%d/data/Lambda.%d.%d.dat' % (i,iNF,alf_info['ncentral']))
  fnmout=('data/Filter.%d.dat' % iNF)
  subprocess.call(['mpirun','-n','1','-bynode','--bind-to','none','-x','OMP_NUM_THREADS=1',exe,str(FREQ),fnmin,fnmout])
  print("WARNING: direct standard error somewhere")
  time.sleep(15)

def MomentDCA(i,iNF,NF,FREQ,engine='charmm'):
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
  subprocess.call(['mpirun','-n','1','-bynode','--bind-to','none','-x','OMP_NUM_THREADS=1',exe,fnmin,fnmout1,fnmout2])
  print("WARNING: direct standard error somewhere")
  time.sleep(15)

def BSMomentDCA(i,NF,FREQ,engine='charmm'):
  import os, sys, shutil, traceback, time, subprocess, random
  import numpy as np
  import alf

  alf_info=alf.initialize_alf_info(engine)

  home_dir=os.getcwd()

  os.chdir('dca%d' % i)

  exe=os.path.dirname(os.path.abspath(__file__))+'/dca/BootstrapMoments.py'
  subprocess.call([exe,str(NF),'data'])
  time.sleep(15)

def LMDCA(i,iBS,NF,FREQ,engine='charmm'):
  import os, sys, shutil, traceback, time, subprocess, random
  import numpy as np
  import alf

  alf_info=alf.initialize_alf_info(engine)

  home_dir=os.getcwd()

  os.chdir('dca%d' % i)

  print("Warning, BS is hardcoded")
  BS=50
  exe=os.path.dirname(os.path.abspath(__file__))+'/dca/LM'
  if iBS==BS:
    fnmout1='data/h.LM.dat'
    fnmout2='data/J.LM.dat'
    fnmin1='data/m1.obs.dat' 
    fnmin2='data/m2.obs.dat'
  else:
    fnmout1=('data/h.bs%d.LM.dat' % iBS)
    fnmout2=('data/J.bs%d.LM.dat' % iBS)
    fnmin1=('data/m1.bs%d.obs.dat' % iBS)
    fnmin2=('data/m2.bs%d.obs.dat' % iBS)
  subprocess.call(['mpirun','-n','1','-bynode','--bind-to','none','-x','OMP_NUM_THREADS=8',exe,fnmout1,fnmout2,fnmin1,fnmin2])
  time.sleep(15)

def PLMDCA(i,iBS,NF,FREQ,engine='charmm'):
  import os, sys, shutil, traceback, time, subprocess, random
  import numpy as np
  import alf

  alf_info=alf.initialize_alf_info(engine)

  home_dir=os.getcwd()

  os.chdir('dca%d' % i)

  print("Warning, BS is hardcoded")
  BS=50
  exe=os.path.dirname(os.path.abspath(__file__))+'/dca/PLM'
  if iBS==BS:
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
  subprocess.call(['mpirun','-n','1','-bynode','--bind-to','none','-x','OMP_NUM_THREADS=1',exe,fnmout1,fnmout2]+fnmsin)
  time.sleep(15)

def FinishDCA(i,NF,FREQ,engine='charmm'):
  import os, sys, shutil, traceback, time, subprocess, random
  import numpy as np
  import alf

  alf_info=alf.initialize_alf_info(engine)

  home_dir=os.getcwd()

  os.chdir('dca%d' % i)

  alf.GetVarianceDCA(alf_info,NF,'data')
  alf.GetModelDCA(alf_info,NF,'data')
