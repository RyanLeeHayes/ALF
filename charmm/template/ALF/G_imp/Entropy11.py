#! /usr/bin/env python

import sys, os
import numpy as np

Nmc=5000000
Nl=20
c=5.5
dl=1.0/Nl

Ndim=int(sys.argv[1])

for i in range(2,Ndim+1):
  G_ki=np.loadtxt('G1_'+str(i)+'.dat')
  h_i=np.sum(np.reshape(np.exp(-G_ki),(Nl,Nl),order='C'),axis=1)
  S_ki=np.log(h_i)
  G_ki=-S_ki+np.log(np.mean(h_i))
  G_ki=np.reshape(G_ki,(Nl,1))
  G_ki=np.tile(G_ki,(1,Nl))
  for j in range(2,Ndim+1):
    G_kj=np.loadtxt('G1_'+str(j)+'.dat')
    h_j=np.sum(np.reshape(np.exp(-G_kj),(Nl,Nl),order='C'),axis=1)
    S_kj=np.log(h_j)
    G_kj=-S_kj+np.log(np.mean(h_j))
    G_kj=np.reshape(G_kj,(1,Nl))
    G_kj=np.tile(G_kj,(Nl,1))
    np.savetxt('G1_'+str(i)+'_'+str(j)+'.dat',G_ki+G_kj)
