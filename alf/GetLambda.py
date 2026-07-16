#! /usr/bin/env python

def GetLambdaCharmm(alf_info,fnmout,fnmsin):
  """
  Reads alchemical trajectories from CHARMM binary format

  This routine is called by the routine GetLambda to read CHARMM binary
  alchemical trajectories and write them to human readable format.

  Parameters
  ----------
  alf_info : dict
      Dictionary of variables alf needs to run
  fnmout : str
      The filename for the human readable output
  fnmsin : list of str
      A list of the filenames of the binary input files
  """

  import os, subprocess

  for fnmin in fnmsin:
    if not os.path.exists(fnmin):
      print('Error, %s does not exist, molecular dynamics probably failed, check run output and run error for clues' % (fnmin,))

  exe=os.path.dirname(os.path.abspath(__file__))+'/bin/GetLambda'
  subprocess.check_output([exe,fnmout,*fnmsin])

def GetLambdaCharmmSlow(alf_info,fnmout,fnmsin):
  """
  Reads alchemical trajectories from CHARMM binary format

  This routine is called by the routine GetLambda to read CHARMM binary
  alchemical trajectories and write them to human readable format.

  Parameters
  ----------
  alf_info : dict
      Dictionary of variables alf needs to run
  fnmout : str
      The filename for the human readable output
  fnmsin : list of str
      A list of the filenames of the binary input files
  """

  import sys, os
  import numpy as np
  from scipy.io import FortranFile

  nblocks=alf_info['nblocks']
  # Lambdas=np.zeros((0,nblocks))

  fpout=open(fnmout,'w')

  for fnmin in fnmsin:
    if not os.path.exists(fnmin):
      print('Error, %s does not exist, molecular dynamics probably failed, check run output and run error for clues' % (fnmin,))
    fp=FortranFile(fnmin,'r')

    # The header and icntrl array are read in as a single record
    # Read the icntrl array (length 20) and extract key variables

    header = (fp.read_record([('hdr',np.bytes_,4),('icntrl',np.int32,20)]))
    hdr = header['hdr'][0]
    icntrl = header['icntrl'][0][:]
    nfile = icntrl[0]     # Total number of dynamcis steps in lambda file
    npriv = icntrl[1]     # Number of steps preceding this run
    nsavl = icntrl[2]     # Save frequency for lambda in file
    nblocks = icntrl[6]   # Total number of blocks = env + subsite blocks
    nsitemld = icntrl[10] # Total number of substitution sites (R-groups) in MSLD

    # Time step for dynamics in AKMA units
    delta4 = (fp.read_record(dtype=np.float32))

    # Title in trajectory file 
    # title = (fp.read_record([('h',np.int32,1),('title',np.bytes_,80)]))[0][1]
    title = (fp.read_record([('h',np.int32),('title',np.bytes_,80)]))[0][1]

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

    np.savetxt(fpout,Lambda,fmt="%10.6f")
    # Lambdas=np.concatenate((Lambdas,Lambda),axis=0)

  fpout.close()
  # np.savetxt(fnmout,Lambdas,fmt="%10.6f")



def GetLambdaBlade(alf_info,fnmout,fnmsin):
  """
  Reads alchemical trajectories from BLaDE binary format

  This routine is called by the routine GetLambda to read standalone BLaDE
  binary alchemical trajectories and write them to human readable format.

  Parameters
  ----------
  alf_info : dict
      Dictionary of variables alf needs to run
  fnmout : str
      The filename for the human readable output
  fnmsin : list of str
      A list of the filenames of the binary input files
  """

  import os
  import numpy as np

  nblocks=alf_info['nblocks']
  Lambdas=np.zeros((0,nblocks))

  # BLaDE writes each record as XDR: one big-endian int32 index followed by
  # nblocks big-endian float32 lambda values, with no padding between fields.
  dtype=np.dtype([('index','>i4'),('lambdas','>f4',(nblocks,))])

  for fnmin in fnmsin:
    if not os.path.exists(fnmin):
      print('Error, %s does not exist, molecular dynamics probably failed, check run output and run error for clues' % (fnmin,))
    fp=open(fnmin,"rb")
    fpdata=fp.read()
    fp.close()
    lines=len(fpdata)//dtype.itemsize
    records=np.frombuffer(fpdata,dtype=dtype,count=lines)
    Lambdas=np.concatenate((Lambdas,records['lambdas'].astype(np.float64)),axis=0)

  np.savetxt(fnmout,Lambdas,fmt="%10.6f")



def GetLambda(alf_info,fnmout,fnmsin):
  """
  Reads alchemical trajectories from binary format

  This routine is called by the routine GetLambdas to read binary
  alchemical trajectories and write them to human readable format. Based
  on the contents of alf_info['engine'], this routine either wraps the
  GetLambdaCharmm routine for reading CHARMM binary alchemical
  trajectories or the GetLambdaBlade routine for reading standalone BLaDE
  binary alchemical trajectories.

  Parameters
  ----------
  alf_info : dict
      Dictionary of variables alf needs to run
  fnmout : str
      The filename for the human readable output
  fnmsin : list of str
      A list of the filenames of the binary input files
  """

  if alf_info['engine'] in ['charmm','bladelib','pycharmm']:
    GetLambdaCharmm(alf_info,fnmout,fnmsin)
  elif alf_info['engine'] in ['blade']:
    GetLambdaBlade(alf_info,fnmout,fnmsin)
  else:
    print("Error: unsupported engine type %s" % alf_info['engine'])
    quit(1)
