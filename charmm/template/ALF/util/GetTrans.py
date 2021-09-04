#! /usr/bin/env python

import numpy as np

nblocks=np.loadtxt('../nblocks',dtype='int')
nsubs=np.loadtxt('../nsubs',dtype='int',ndmin=1)


for tag in range(0,5):

  data=np.loadtxt("data/Lambda."+str(tag)+".0.dat")

  ibuff=0
  for i in range(0,len(nsubs)):

    trans=np.zeros((nsubs[i],nsubs[i]))

    i_curr=-1
    for t in range(0,data.shape[0]): # row in data:
      i_prev=i_curr
      for j in range(0,nsubs[i]):
        if data[t,j+ibuff]>0.8:
          i_curr=j;
      if i_prev>=0 and i_prev!=i_curr:
        trans[i_prev,i_curr]+=1

    for j in range(0,nsubs[i]):
      total=0
      line=""
      for k in range(0,nsubs[i]):
        total+=trans[j,k]
        line+=(" %3d" % (trans[j,k],))
      line+=("   (%4d)" % (total,))
      print line
    print "\n"

    ibuff+=nsubs[i]
