#! /usr/bin/env python

def RunWham(nf,nts0,nts1):
  import os
  import ctypes

  if not os.path.exists('multisite'):
    os.mkdir('multisite')

  # ../dWHAMdV_mso/wham $1 Energy Lambda $2 $3
  whamlib = ctypes.CDLL('/home/rhaye/41_MSLD/51_ALF-3.1/v3.1.3/alf/alf/wham/libwham.so')
  pywham = whamlib.main
  pywham.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.c_int]
  pywham(nf,nts0,nts1)
