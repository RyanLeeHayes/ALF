
def initialize(engine='charmm',minimize=True):
  import os, sys, shutil, traceback, time, subprocess, random
  import numpy as np
  import alf

  alf_info=alf.initialize_alf_info(engine)

  home_dir=os.getcwd()

  alf.InitVars(alf_info,minimize=minimize)
  os.chdir(home_dir)
