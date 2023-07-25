#! /usr/bin/env python

def GetEnergy(alf_info,Fi,Ff,skipE=1):
  """
  Compute energies for wham/mbar reweighting

  This routine selects several simulations for reweighting with wham/mbar
  and computes the relative bias energy from each simulation on the
  alchemical trajectory for each simulation. Different replicas in
  replica exchange are considered different simulations. This routine
  should be run inside analysis[Ff]. Biases are taken from the apropriate
  analysis[i] directory and alchemical trajectories are taken from the
  output of the routine alf.GetLambdas in analysis[i]/data. Copies
  alchemical trajectories from data directories analysis[Fi]/data to
  analysis[Ff]/data into analysis[Ff]/Lambda/Lambda[*].dat and copies
  bias energies into analysis[Ff]/Energy/ESim[*].dat . Energies are
  computed relative to the central replica bias in the [Ff] cycle.

  This routine can be called during flattening or production. During
  flattening, Ff-Fi=4 is recommended (include 5 cycles in analysis),
  during production Ff-Fi=0 is recommended (include 1 cycle in analysis).
  Flattening versus production is detected by the presence of the run[Ff]
  directory. During production this directory does not exist, instead,
  run[Ff]a, and the other independent trial directories exist.

  Parameters
  ----------
  alf_info : dict
      Dictionary of variables alf needs to run
  Fi : int
      The first cycle of alf to include in analysis (inclusive)
  Ff : int
      The final cycle of alf to include in analysis (inclusive)
  skipE : int, optional
      In longer production runs the number of lambda samples may require
      significant amounts of memory to store and analyze. Only alchemical
      frames with index modulus skipE equal to skipE-1 are analyzed.
      (default is 1 to analyze all frames) 
  """

  import sys, os
  import numpy as np

  # Fi=int(sys.argv[1])
  # Ff=int(sys.argv[2])
  NF=Ff-Fi+1;

  # if len(sys.argv)>3:
  #   skipE=int(sys.argv[3])
  # else:
  #   skipE=1

  alphabet='abcdefghijklmnopqrstuvwxyz'

  if os.path.isdir('../run'+str(Ff)):
    production=False
    ndupl=1
  else:
    production=True
    ndupl=0
    for i in range(0,len(alphabet)):
      if os.path.isdir('../run'+str(Ff)+alphabet[i]):
        ndupl+=1
    if ndupl==0:
      print("Error, not flattening or production")
      quit()

  nblocks=alf_info['nblocks']
  nsubs=alf_info['nsubs']
  nreps=alf_info['nreps']
  ncentral=alf_info['ncentral']

  b_shift=np.loadtxt('../nbshift/b_shift.dat')
  c_shift=np.loadtxt('../nbshift/c_shift.dat')
  x_shift=np.loadtxt('../nbshift/x_shift.dat')
  s_shift=np.loadtxt('../nbshift/s_shift.dat')

  Lambda=[]
  b=[]
  c=[]
  x=[]
  s=[]
  for i in range(0,NF):
    DIR='../analysis'+str(Fi+i)
    for j in range(0,ndupl):
      for k in range(0,nreps):
        if production:
          Lambda.append(np.loadtxt(DIR+'/data/Lambda.'+str(j)+'.'+str(k)+'.dat')[(skipE-1)::skipE,:])
        else:
          Lambda.append(np.loadtxt(DIR+'/data/Lambda.'+str(j)+'.'+str(k)+'.dat'))
        b.append(np.loadtxt(DIR+'/b_prev.dat')+b_shift*(k-ncentral))
        c.append(np.loadtxt(DIR+'/c_prev.dat')+c_shift*(k-ncentral))
        x.append(np.loadtxt(DIR+'/x_prev.dat')+x_shift*(k-ncentral))
        s.append(np.loadtxt(DIR+'/s_prev.dat')+s_shift*(k-ncentral))

  if not os.path.isdir('Lambda'):
    os.mkdir('Lambda')
  if not os.path.isdir('Energy'):
    os.mkdir('Energy')

  E=[]
  for i in range(0,NF*ndupl*nreps):
    E.append([])
    for j in range(0,NF*ndupl*nreps):
      bi=b[i]
      ci=c[i]
      xi=x[i]
      si=s[i]
      Lj=Lambda[j]
      E[-1].append(np.reshape(np.dot(Lj,-bi),(-1,1)))
      E[-1][-1]+=np.sum(np.dot(Lj,-ci)*Lj,axis=1,keepdims=True)
      E[-1][-1]+=np.sum(np.dot(1-np.exp(-5.56*Lj),-xi)*Lj,axis=1,keepdims=True)
      E[-1][-1]+=np.sum(np.dot(Lj/(Lj+0.017),-si)*Lj,axis=1,keepdims=True)

  for i in range(0,NF*ndupl*nreps):
    Ei=E[nreps*(NF*ndupl-1)+ncentral][i]
    for j in range(0,NF*ndupl*nreps):
      Ei=np.concatenate((Ei,E[j][i]),axis=1)
    np.savetxt('Energy/ESim'+str(i+1)+'.dat',Ei,fmt='%12.5f')

  for i in range(0,NF*ndupl*nreps):
    Li=Lambda[i]
    np.savetxt('Lambda/Lambda'+str(i+1)+'.dat',Li,fmt='%10.6f')
