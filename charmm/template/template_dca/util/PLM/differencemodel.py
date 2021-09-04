#! /usr/bin/env python

import sys, os
import numpy as np
import math

directories=['../../rbd/dca214/data','../../complex/dca214/data']
bs=50

if not os.path.isdir('model'):
  os.mkdir('model')

nsubs=np.loadtxt(directories[1]+'/../../nsubs',dtype='int',ndmin=1)

h=[[],[]]
J=[[],[]]
density=[[],[]]

for i in range(2):
  h[i]=np.loadtxt(directories[i]+'/h.bias.dat')
  J[i]=np.loadtxt(directories[i]+'/J.bias.dat')
  density[i]=np.loadtxt(directories[i]+'/Density.dat')
h_diff=h[1]-h[0]
J_diff=J[1]-J[0]
density_diff=density[1]-density[0]
dQ=np.loadtxt('q.dat')
gamma_S=0.76414
kelectric=332.0716
correction=-(2.0*np.pi/3.0)*kelectric*gamma_S*dQ*density_diff

np.savetxt('model/h.bias.dat',h_diff,fmt=' %10.6f')
np.savetxt('model/J.bias.dat',J_diff,fmt=' %10.6f')
np.savetxt('model/h.charge.dat',correction,fmt=' %10.6f')

for i in range(2):
  h[i]=np.loadtxt(directories[i]+'/h.model.dat')
  J[i]=np.loadtxt(directories[i]+'/J.model.dat')
h_diff=h[1]-h[0]+correction
J_diff=J[1]-J[0]

np.savetxt('model/h.model.dat',h_diff,fmt=' %10.6f')
np.savetxt('model/J.model.dat',J_diff,fmt=' %10.6f')

block0i=0
for i in range(len(nsubs)):
  block0j=0
  for j in range(len(nsubs)):
    for isub in range(nsubs[i]):
      J0=J_diff[block0i+isub,block0j]
      J_diff[block0i+isub,block0j:(block0j+nsubs[j])]-=J0
      h_diff[block0i+isub]+=0.5*J0
    for jsub in range(nsubs[j]):
      J0=J_diff[block0i,block0j+jsub]
      J_diff[block0i:(block0i+nsubs[i]),block0j+jsub]-=J0
      h_diff[block0j+jsub]+=0.5*J0
    block0j+=nsubs[j]
  block0i+=nsubs[i]
block0i=0
for i in range(len(nsubs)):
  h0=h_diff[block0i]
  h_diff[block0i:(block0i+nsubs[i])]-=h0
  block0i+=nsubs[i]
np.savetxt('model/Fields.model.dat',h_diff[1::2],fmt=' %10.6f')
np.savetxt('model/Couplings.model.dat',J_diff[1::2,1::2],fmt=' %10.6f')

for iB in range(bs):
  for i in range(2):
    h[i]=np.loadtxt(directories[i]+'/h.bs'+str(iB)+'.model.dat')
    J[i]=np.loadtxt(directories[i]+'/J.bs'+str(iB)+'.model.dat')
  h_diff=h[1]-h[0]+correction
  J_diff=J[1]-J[0]

  np.savetxt('model/h.bs'+str(iB)+'.model.dat',h_diff,fmt=' %10.6f')
  np.savetxt('model/J.bs'+str(iB)+'.model.dat',J_diff,fmt=' %10.6f')
