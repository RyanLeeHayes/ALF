#! /usr/bin/env python

def GetStepsCharmm(alf_info,fnm):
  import numpy as np
  from scipy.io import FortranFile

  try:
    fp=FortranFile(fnm,'r')
  except IOError:
    try:
      fp=FortranFile(fnm+'_0','r')
    except IOError:
      return(0)

  try:
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

  except:
    return(0)

  # Can't trust it. Make sure the frames are actually there.
  # print(nfile)

  nread=0
  for i in range(nfile):
    try:
      # Read a line of lambda values
      lambdav = (fp.read_record(dtype=np.float32))
      theta = (fp.read_record(dtype=np.float32))
      # Lambda[i,:]=lambdav[1:]
    except:
      break
    nread+=1

  return(nread)

  fp.close()



def GetStepsBlade(alf_info,fnm):
  import numpy as np
  from xdrlib import Unpacker
  from xdrlib import Packer

  nblocks=alf_info['nblocks']

  p=Packer()
  p.pack_int(0)
  for j in range(0,nblocks):
    p.pack_float(0)
  linewidth=len(p.get_buffer())

  try:
    fp=open(fnm,"rb")
    fpdata=fp.read()
    lines=len(fpdata)//linewidth
    fp.close()
  except:
    try:
      fp=open(fnm+'_0',"rb")
      fpdata=fp.read()
      lines=len(fpdata)//linewidth
      fp.close()
    except:
      return(0)

  return(lines)



def GetSteps(alf_info,fnm):
  if alf_info['engine'] in ['charmm','bladelib']:
    return GetStepsCharmm(alf_info,fnm)
  elif alf_info['engine'] in ['blade']:
    return GetStepsBlade(alf_info,fnm)
  else:
    print("Error: unsupported engine type %s" % alf_info['engine'])
    quit()
