#! /usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

# import alf
# from prep.alf_info import alf_info
# alf.initialize('blade')

Emid=np.arange(1.0/800,1,1.0/400)
Emid2=np.arange(1.0/40,1,1.0/20)
EmidX,EmidY=np.meshgrid(Emid2,Emid2)

# nsubs=alf_info['nsubs']
# nblocks=alf_info['nblocks']
nsubs=np.loadtxt('../prep/nsubs',dtype='int',ndmin=1)
nblocks=np.sum(nsubs)

ntersite=[0,0]
msprof=ntersite[1]

iG=1
iF=1

G1={}
G12={}
G2={}
G11={}
iblock=0
for isite in range(len(nsubs)):
  jblock=iblock
  for jsite in range(isite,len(nsubs)):
    if isite==jsite:

      # LEG={};
      plt.figure(iF)
      iF=iF+1
      # hold off
      LEG=[]
      for i in range(nsubs[isite]):
        # LEG{i}=num2str(i);
        LEG.append(str(i+1))
        G1[i+iblock]=np.loadtxt('multisite/G'+str(iG)+'.dat')
        iG=iG+1
        # subplot(1,1,1)
        # plot(Emid,G1{i+iblock})
        plt.plot(Emid,G1[i+iblock])
        plt.xlabel('lambda')
        plt.ylabel('Free energy [kcal/mol]')
        # hold on
      # hold off
      # legend('location','best',LEG)
      plt.legend(LEG)

      plt.figure(iF)
      iF=iF+1
      for i in range(nsubs[isite]):
        G12[i+iblock]={}
        for j in range((i+1),nsubs[isite]):
          G12[i+iblock][j+iblock]=np.loadtxt('multisite/G'+str(iG)+'.dat')
          iG=iG+1
          # subplot(nsubs(isite)-1,nsubs(isite)-1,(i-1)*(nsubs(isite)-1)+(j-1))
          # hold off
          # plot(Emid,G12{i+iblock,j+iblock})
          plt.subplot(nsubs[isite]-1,nsubs[isite]-1,i*(nsubs[isite]-1)+j)
          plt.plot(Emid,G12[i+iblock][j+iblock])

      if nsubs[isite]>2:
        plt.figure(iF)
        iF=iF+1
        for i in range(nsubs[isite]):
          G2[i+iblock]={}
          for j in range((i+1),nsubs[isite]):
            G2[i+iblock][j+iblock]=np.reshape(np.loadtxt('multisite/G'+str(iG)+'.dat'),(20,20))
            iG=iG+1
            # subplot(nsubs(isite)-1,nsubs(isite)-1,(i-1)*(nsubs(isite)-1)+(j-1))
            # hold off
            # surf(Emid2,Emid2,G2{i+iblock,j+iblock})
            pltax=plt.subplot(nsubs[isite]-1,nsubs[isite]-1,i*(nsubs[isite]-1)+j,projection='3d')
            pltax.plot_surface(EmidX,EmidY,G2[i+iblock][j+iblock])

    elif msprof:

      plt.figure(iF)
      iF=iF+1
      for i in range(nsubs[isite]):
        G11[i+iblock]={}
        for j in range(nsubs[jsite]):
          G11[i+iblock][j+jblock]=np.reshape(np.loadtxt('multisite/G'+str(iG)+'.dat'),(20,20))
          iG=iG+1
          # subplot(nsubs(isite),nsubs(jsite),(i-1)*(nsubs(jsite))+(j))
          # hold off
          # surf(Emid2,Emid2,G11{i+iblock,j+jblock})
          plt.subplot(nsubs[isite],nsubs[jsite],i*nsubs[jsite]+j)
          plt.plot_surface(EmidX,EmidY,G11[i+iblock][j+jblock])

    jblock=jblock+nsubs[jsite]
  iblock=iblock+nsubs[isite]

plt.show()

