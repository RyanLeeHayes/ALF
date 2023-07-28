#! /usr/bin/env python

def PlotFreeEnergy5(directory=None,ntersite=[0,0]):
  """
  Plots free energy profiles for a cycle of alf

  This should be run after alf.postprocess has completed because it
  displays the free energy profiles computed by alf.RunWham

  Parameters
  ----------
  directory : str, optional
      A string for the analysis directory of the cycle of interest. If
      blank, analysis will be performed in this directory. (default is
      None)
  ntersite : list of int, optional
      The ntersite list used during postprocessing on this cycle of alf.
      If the second element of the list is incorrect, multisite systems
      will not display correctly. (default is [0,0])
  """

  import sys, os
  import numpy as np
  import matplotlib.pyplot as plt

  # import alf
  # from prep.alf_info import alf_info
  # alf.initialize('blade')

  Emid=np.arange(1.0/800,1,1.0/400)
  Emid2=np.arange(1.0/40,1,1.0/20)
  EmidX,EmidY=np.meshgrid(Emid2,Emid2)

  DIR=os.getcwd()
  if directory:
    os.chdir(directory)

  # nsubs=alf_info['nsubs']
  # nblocks=alf_info['nblocks']
  nsubs=np.loadtxt('nsubs',dtype='int',ndmin=1)
  nblocks=np.sum(nsubs)

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

        plt.figure(iF)
        iF=iF+1
        LEG=[]
        for i in range(nsubs[isite]):
          LEG.append(str(i+1))
          G1[i+iblock]=np.loadtxt('multisite/G'+str(iG)+'.dat')
          iG=iG+1
          # plot(Emid,G1{i+iblock})
          plt.plot(Emid,G1[i+iblock])
          plt.xlabel('lambda')
          plt.ylabel('Free energy [kcal/mol]')
        plt.title('Site %d:   1D profiles' % (isite+1,))
        plt.legend(LEG)

        pltfg=plt.figure(iF)
        iF=iF+1
        for i in range(nsubs[isite]):
          G12[i+iblock]={}
          for j in range((i+1),nsubs[isite]):
            G12[i+iblock][j+iblock]=np.loadtxt('multisite/G'+str(iG)+'.dat')
            iG=iG+1
            # plot(Emid,G12{i+iblock,j+iblock})
            plt.subplot(nsubs[isite]-1,nsubs[isite]-1,i*(nsubs[isite]-1)+j)
            plt.plot(Emid,G12[i+iblock][j+iblock])
        plt.suptitle('Site %d:   Transition profiles' % (isite+1,))
        pltfg.supxlabel('Column+1 on at x=0')
        pltfg.supylabel('Row substituent on at x=1')

        if nsubs[isite]>2:
          plt.figure(iF)
          iF=iF+1
          for i in range(nsubs[isite]):
            G2[i+iblock]={}
            for j in range((i+1),nsubs[isite]):
              G2[i+iblock][j+iblock]=np.reshape(np.loadtxt('multisite/G'+str(iG)+'.dat'),(20,20))
              iG=iG+1
              # surf(Emid2,Emid2,G2{i+iblock,j+iblock})
              pltax=plt.subplot(nsubs[isite]-1,nsubs[isite]-1,i*(nsubs[isite]-1)+j,projection='3d')
              pltax.plot_surface(EmidX,EmidY,G2[i+iblock][j+iblock])
          plt.suptitle('Site %d:   2D profiles' % (isite+1,))

      elif msprof:

        plt.figure(iF)
        iF=iF+1
        for i in range(nsubs[isite]):
          G11[i+iblock]={}
          for j in range(nsubs[jsite]):
            G11[i+iblock][j+jblock]=np.reshape(np.loadtxt('multisite/G'+str(iG)+'.dat'),(20,20))
            iG=iG+1
            # surf(Emid2,Emid2,G11{i+iblock,j+jblock})
            pltax=plt.subplot(nsubs[isite],nsubs[jsite],i*nsubs[jsite]+j+1,projection='3d')
            pltax.plot_surface(EmidX,EmidY,G11[i+iblock][j+jblock])
        plt.suptitle('Sites %d and %d:   2D coupling profiles' % (isite+1,jsite+1))

      jblock=jblock+nsubs[jsite]
    iblock=iblock+nsubs[isite]

  plt.show()
  os.chdir(DIR)

if __name__ == '__main__':
  import sys
  if len(sys.argv)==1:
    PlotFreeEnergy5()
  elif len(sys.argv)==2:
    PlotFreeEnergy5(sys.argv[1])
  elif len(sys.argv)==3:
    PlotFreeEnergy5(ntersite=[int(sys.argv[1]),int(sys.argv[2])])
  elif len(sys.argv)==4:
    PlotFreeEnergy5(sys.argv[1],ntersite=[int(sys.argv[2]),int(sys.argv[3])])
