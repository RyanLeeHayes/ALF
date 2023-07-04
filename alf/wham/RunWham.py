#! /usr/bin/env python

def RunWham(nf,temp,nts0,nts1):
  import os
  import ctypes

  if not os.path.exists('multisite'):
    os.mkdir('multisite')

  # ../dWHAMdV_mso/wham $1 Energy Lambda $2 $3
  whamlib = ctypes.CDLL(os.path.dirname(__file__)+'/libwham.so')
  pywham = whamlib.wham # declare as a c function to prevent c++ name mangling
  pywham.argtypes=[ctypes.c_int,ctypes.c_double,ctypes.c_int,ctypes.c_int]
  pywham(nf,temp,nts0,nts1)
