
def initialize(engine='charmm',minimize=True):
  """
  set up the initial data structures required by runflat

  This function sets up the files required for the first run of runflat.
  It expects a prep directory formatted according to the specifications
  in README.md. The data These files required by runflat are analysis0,
  which contains starting values of 0 for all the bias parameters,
  variables1.inp which formats those biases to be read by the MD engine,
  and nbshift, which described how the bias parameters change for
  neighboring replicas when using replica exchange. Note that while
  the contents of nbshift are irrelevant if not using replica exchange,
  GetEnergy still requires that they be present. This routine calls the
  InitVars routine.

  Parameters
  ----------
  engine : str, optional
      The molecular dynamics engine string, see help(alf) for allowed
      values. (default is 'charmm')
  minimize : bool, optional
      A flag for whether to perform a minimization on the first cycle of
      ALF. True performs the minimization. If minimization is not
      performed, the default scripts expect that the prep directory
      contain a starting structure named minimized.pdb or minimized.crd.
      (default is True)
  """

  import os, sys, shutil, traceback, time, subprocess, random
  import numpy as np
  import alf

  alf_info=alf.initialize_alf_info(engine)

  home_dir=os.getcwd()

  alf.InitVars(alf_info,minimize=minimize)
  os.chdir(home_dir)
