#! /usr/bin/env python

def SetGimp(alf_info):
  """
  Creates the G_imp directory

  The G_imp directory defines the free energy profiles of a flat landscape
  so this free energy is subtracted from all free energy profiles before
  flattening them. The free energy has units of kT. Run from the analysis0
  directory.

  Parameters
  ----------
  alf_info : dict
      Dictionary of variables alf needs to run
  """

  import os, shutil, os.path
  import subprocess

  # Don't generate G_imp if using custom implicit constraints, the user will supply it
  if alf_info['impcons']=='custom':
    return

  if not os.path.exists('../G_imp'):
    os.mkdir('../G_imp')
  shutil.copy('../prep/nsubs','nsubs')
  fpout=open('../G_imp/output','w')
  fperr=open('../G_imp/error','w')
  if alf_info['loss']=='linear2018':
    nbins2=20
  elif alf_info['loss']=='nonlinear2024':
    nbins2=16
  if alf_info['impcons'] in ['fnex2011','fnexdozen2024']:
    options=[str(alf_info['impconsopt']['c'])]
  elif alf_info['impcons']=='fnpwise2026':
    options=[str(alf_info['impconsopt']['kbeta']),str(alf_info['impconsopt']['width'])]
  subprocess.call(['time',os.path.dirname(__file__)+'/bin/impcons',str(10000000),'gimp',str(nbins2),alf_info['impcons']]+options,stdout=fpout,stderr=fperr)
  fpout.close()
  fperr.close()

def linear2018(alf_info,N,ntersite,fpout,fperr):
  """
  Wrapper for runlinear2018

  A wrapper is needed because when runlinear2018 is called directly, it
  doesn't release the GPU upon completion, which then precludes the
  molecular dynamics engine from using the GPU, so it must be called with
  subprocess.

  Parameters
  ----------
  alf_info : dict
      Dictionary of variables alf needs to run
  N : int
      The number of flattening cycles or independent trials to combine
      with wham
  ntersite : list
      The first element of the ntersite list controls whether intersite
      coupling parameters are optimized. 0 is no, 1 optimizes c, x, and s
      coupling biases, 2 only optimizes c coupling biases. The second
      element of the ntersite list, controlling whether intersite free
      energy profiles are flattened. 0 is no, 1 is yes.
  fpout : file
      File pointer to write standard output to
  fperr : file
      File pointer to write standard error to
  """

  import subprocess, shutil

  subprocess.call(['time',shutil.which('python'),'-c','import alf; alf.runlinear2018("%s",%d,%f,%d,%d)' % (alf_info['bias'],N*alf_info['nreps'],alf_info['temp'],ntersite[0],ntersite[1])],stdout=fpout,stderr=fperr)

def runlinear2018(bias,nf,temp,nts0,nts1):
  """
  Compute profiles and how changes in biases will affect penalty function

  This is a python wrapper for the GPU code in alf/bin/linear2018. This
  code uses the output of alf.GetEnergy to combine samples from
  simulations with different biases into a coherent model of the free
  energy landscape using wham or mbar equations. Free energy profiles are
  then computed for various projections of this free energy landscape as a
  function of various lambda variables. The square of the deviation of
  these profiles from their average values is the penalty function
  optimized by ALF using a linear approximation. The linear effect of
  changes in bias parameters on the free energy profiles is computed, and
  the squaring by the penalty function gives a penalty function that is a
  quadratic function of changes in the bias parameters. The second
  derivatives of this penalty function are saved to
  analysis[i]/multisite/C.dat and the first derivatives are save to
  analysis[i]/multisite/V.dat.

  This routine should be run from the analysis[i] directory, and should be
  followed by a call to GetFreeEnergy5 to invert the C.dat and V.dat
  matrices. This routine is called by alf.runflat and alf.postprocess.

  alf.PlotFreeEnergy5 may be used to plot the free energy profiles
  computed by wham to assess convergence.

  Parameters
  ----------
  bias : string
      The bias style (bcxs2018 or bcxstu2026)
  nf : int
      The number of flattening cycles or independent trials to combine
      with wham
  temp : float
      The temperature at which the simulation was run
  nts0 : int
      The first element of the ntersite list, controlling whether
      intersite coupling parameters are optimized. 0 is no, 1 optimizes
      c, x, and s coupling biases, 2 only optimizes c coupling biases
  nts1 : int
      The second element of the ntersite list, controlling whether
      intersite free energy profiles are flattened. 0 is no, 1 is yes.
  """

  import os
  import ctypes

  if not os.path.exists('multisite'):
    os.mkdir('multisite')

  # ../dWHAMdV_mso/wham $1 Energy Lambda $2 $3
  whamlib = ctypes.CDLL(os.path.dirname(__file__)+'/bin/liblinear2018.so')
  pywham = whamlib.wham # declare as a c function to prevent c++ name mangling
  pywham.argtypes=[ctypes.c_char_p,ctypes.c_int,ctypes.c_double,ctypes.c_int,ctypes.c_int]
  pywham(bytes(bias,encoding='ascii'),nf,temp,nts0,nts1)

def nonlinear2024(alf_info,N,ntersite,fpout,fperr):
  """
  Compute profiles and how changes in biases affect 2024 loss function

  Parameters
  ----------
  alf_info : dict
      Dictionary of variables alf needs to run
  N : int
      The number of flattening cycles or independent trials to combine
      with wham
  ntersite : list
      The first element of the ntersite list controls whether intersite
      coupling parameters are optimized. 0 is no, 1 optimizes c, x, and s
      coupling biases, 2 only optimizes c coupling biases. The second
      element of the ntersite list, controlling whether intersite free
      energy profiles are flattened. 0 is no, 1 is yes.
  fpout : file
      File pointer to write standard output to
  fperr : file
      File pointer to write standard error to
  """

  import subprocess
  import os
  import numpy as np

  if not os.path.exists('multisite'):
    os.mkdir('multisite')

  # Run wham and get the weights
  subprocess.call(['time',os.path.dirname(__file__)+'/bin/whamweight',str(N*alf_info['nreps']),str(alf_info['temp']),str(ntersite[0]),str(ntersite[1])],stdout=fpout,stderr=fperr)

  # Generate monte carlo samples
  samples=len(np.loadtxt('Lambda/weight.dat'))
  if alf_info['impcons'] in ['fnex2011','fnexdozen2024','fnpwise2026']:
    if alf_info['impcons'] in ['fnex2011','fnexdozen2024']:
      options=[str(alf_info['impconsopt']['c'])]
    elif alf_info['impcons']=='fnpwise2026':
      options=[str(alf_info['impconsopt']['kbeta']),str(alf_info['impconsopt']['width'])]
    subprocess.call(['time',os.path.dirname(__file__)+'/bin/impcons',str(samples),'sample',alf_info['impcons']]+options,stdout=fpout,stderr=fperr)
  elif alf_info['impcons']=='custom':
    # still have to generate mc_Lambda.dat by tiling mc_Lambda.dat from G_imp
    mc_Lambda=np.loadtxt('../G_imp/mc_Lambda.dat')
    mc_Lambdas=np.zeros((0,alf_info['nblocks']))
    while mc_Lambdas.shape[0]<samples:
      mc_Lambdas=np.vstack((mc_Lambdas,mc_Lambda))
    np.savetxt('mc_Lambda.dat',mc_Lambdas[0:samples,:])

  # Run loss function
  subprocess.call(['time',os.path.dirname(__file__)+'/bin/nonlinear2024',alf_info['bias'],str(alf_info['temp']),str(ntersite[0]),str(ntersite[1])],stdout=fpout,stderr=fperr)
