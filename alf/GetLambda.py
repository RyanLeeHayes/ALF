#! /usr/bin/env python

def GetLambdaCharmm(alf_info,fnmout,fnmsin):
  import sys, os
  import numpy as np
  from scipy.io import FortranFile

  nblocks=alf_info['nblocks']
  Lambdas=np.zeros((0,nblocks))

  for fnmin in fnmsin:
    if not os.path.exists(fnmin):
      print('Error, %s does not exist, molecular dynamics probably failed, check run output and run error for clues' % (fnmin,))
    fp=FortranFile(fnmin,'r')

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
    # title = (fp.read_record([('h',np.int32,1),('title',np.string_,80)]))[0][1]
    title = (fp.read_record([('h',np.int32),('title',np.string_,80)]))[0][1]

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

  np.savetxt(fnmout,Lambdas,fmt="%10.6f")



def GetLambdaBlade(alf_info,fnmout,fnmsin):
  import sys, os
  import numpy as np
  from xdrlib import Unpacker
  from xdrlib import Packer

  nblocks=alf_info['nblocks']
  Lambdas=np.zeros((0,nblocks))

  p=Packer()
  p.pack_int(0)
  for j in range(0,nblocks):
    p.pack_float(0)
  linewidth=len(p.get_buffer())

  for fnmin in fnmsin:
    if not os.path.exists(fnmin):
      print('Error, %s does not exist, molecular dynamics probably failed, check run output and run error for clues' % (fnmin,))
    # Lambda=np.loadtxt(sys.argv[ifp])
    fp=open(fnmin,"rb")
    fpdata=fp.read()
    lines=len(fpdata)//linewidth
    fp.close()
    Lambda=np.zeros((lines,nblocks))
    p=Unpacker(fpdata)
    for i in range(0,lines):
      p.unpack_int()
      for j in range(0,nblocks):
        Lambda[i,j]=p.unpack_float()
    Lambdas=np.concatenate((Lambdas,Lambda),axis=0)

  np.savetxt(fnmout,Lambdas,fmt="%10.6f")



def GetLambda(alf_info,fnmout,fnmsin):
  if alf_info['engine'] in ['charmm','bladelib']:
    GetLambdaCharmm(alf_info,fnmout,fnmsin)
  elif alf_info['engine'] in ['blade']:
    GetLambdaBlade(alf_info,fnmout,fnmsin)
  else:
    print("Error: unsupported engine type %s" % alf_info['engine'])
    quit()
