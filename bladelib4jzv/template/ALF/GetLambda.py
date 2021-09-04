#! /usr/bin/env python

import sys
import numpy as np
from scipy.io import FortranFile

if len(sys.argv) < 3:
  print("Error: need an input and output filename")
  quit()

nblocks=np.loadtxt('../nblocks',dtype='int')+1
Lambdas=np.zeros((0,nblocks-1))

for ifp in range(2,len(sys.argv)):
  fp=FortranFile(sys.argv[ifp],'r')

  # The header and icntrl array are read in as a single record
  # Read the icntrl array (length 20) and extract key variables

  header = (fp.read_record([('hdr',np.string_,4),('icntrl',np.int32,20)]))
  hdr = header['hdr'][0]
  icntrl = header['icntrl'][0][:]
  nfile = icntrl[0]     # Total number of dynamcis steps in lambda file
  npriv = icntrl[1]     # Number of steps preceding this run
  nsavl = icntrl[2]     # Save frequency for lambda in file
  nblocks = icntrl[6]   # Total number of blocks = env + subsite blocks
  nsitemld = icntrl[10] # Total number of substitution sites (R-groups) in MSLD

  # Time step for dynamics in AKMA units
  delta4 = (fp.read_record(dtype=np.float32))

  # Title in trajectoory file 
  title = (fp.read_record([('h',np.int32,1),('title',np.string_,80)]))[0][1]

  # Unused in current processing
  nbiasv = (fp.read_record(dtype=np.int32))
  junk = (fp.read_record(dtype=np.float32))

  # Array (length nblocks) indicating which subsites below
  # to which R-substitiution site
  isitemld = (fp.read_record(dtype=np.int32))

  # Temeprature used in lambda dynamics thermostat
  temp = (fp.read_record(dtype=np.float32))

  # Unsed data for this processing
  junk3 = (fp.read_record(dtype=np.float32))

  Lambda=np.zeros((nfile,nblocks-1))

  for i in range(nfile):
    # Read a line of lambda values
    lambdav = (fp.read_record(dtype=np.float32))
    theta = (fp.read_record(dtype=np.float32))
    Lambda[i,:]=lambdav[1:]

  fp.close()

  Lambdas=np.concatenate((Lambdas,Lambda),axis=0)

np.savetxt(sys.argv[1],Lambdas,fmt="%10.6f")


