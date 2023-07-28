#! /usr/bin/env python

def RunWham(nf,temp,nts0,nts1):
  """
  Compute profiles and how changes in biases will affect penalty function

  This is a python wrapper for the GPU code in alf/wham/wham. This code
  uses the output of alf.GetEnergy to combine samples from simulations
  with different biases into a coherent model of the free energy
  landscape using wham or mbar equations. Free energy profiles are then
  computed for various projections of this free energy landscape as a
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
      The seconnd element of the ntersite list, controlling whether
      intersite free energy profiles are flattened. 0 is no, 1 is yes.
  """

  import os
  import ctypes

  if not os.path.exists('multisite'):
    os.mkdir('multisite')

  # ../dWHAMdV_mso/wham $1 Energy Lambda $2 $3
  whamlib = ctypes.CDLL(os.path.dirname(__file__)+'/libwham.so')
  pywham = whamlib.wham # declare as a c function to prevent c++ name mangling
  pywham.argtypes=[ctypes.c_int,ctypes.c_double,ctypes.c_int,ctypes.c_int]
  pywham(nf,temp,nts0,nts1)
