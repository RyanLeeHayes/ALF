#! /usr/bin/env python

def GetTrans(istep,ndupl=None,engine='charmm',lc=0.8):
  """
  Count transitions in alchemical trajectories

  Run from the main directory. Counts transitions in the
  analysis[istep]/data/Lambda.[id].0.dat human readable alchemical
  trajectory files at the 0.99 lambda cutoff threshold, where [istep] is
  the cycle of alf and [id] the the independent trial index.

  Parameters
  ----------
  istep : int
      The cycle of alf for which to count transitions
  ndupl : int, optional
      The number of independent trials of production. If missing, this
      will be treated as a flattening run. (default is None)
  engine : str, optional
      The engine used to run molecular dynamics. (defaul is 'charmm')
  lc : float, optional
      The lambda cutoff above which lambda must rise for a transition to
      count. 0.8 is a common value, another common value is 0.99, which
      is a more stringent lambda cutoff that results in transition rates a
      factor of two or three lower
  """

  import numpy as np
  import alf

  alf_info=alf.initialize_alf_info(engine)

  nsubs=alf_info['nsubs']
  nblocks=alf_info['nblocks']

  if ndupl==None:
    ndupl=1

  for tag in range(0,ndupl):

    data=np.loadtxt("analysis"+str(istep)+"/data/Lambda."+str(tag)+".0.dat")

    ibuff=0
    for i in range(0,len(nsubs)):

      trans=np.zeros((nsubs[i],nsubs[i]),dtype='int')

      i_curr=-1
      for t in range(0,data.shape[0]): # row in data:
        i_prev=i_curr
        for j in range(0,nsubs[i]):
          if data[t,j+ibuff]>lc:
            i_curr=j;
        if i_prev>=0 and i_prev!=i_curr:
          trans[i_prev,i_curr]+=1

      print("Trial %d, site %d" % (tag,i))
      print("Transition matrix")
      print(trans)
      print("Transitions from")
      print(np.sum(trans,axis=1))
      print("Transition to")
      print(np.sum(trans,axis=0))

      ibuff+=nsubs[i]
