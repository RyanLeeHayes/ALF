#! /usr/bin/env python

import numpy as np

exptags=['ecoli_000000000000000',
'ancccons_111111111111111',
'MostStableChimera_010110100111011',
'LeastStableChimera_100001011001001',
'MostStableSingle_R41L_000010000000000',
'2ndMostStableSingle_H62P_000000100000000',
'LeastStableSingle_I66T_000000010000000']

expseq=np.array([[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
[0,1,0,1,1,0,1,0,0,1,1,1,0,1,1],
[1,0,0,0,0,1,0,1,1,0,0,1,0,0,1],
[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0]])

expresults=np.array([
[ 0.00,0.40],
[-2.13,0.60],
[-1.37,0.32],
[ 0.38,0.53],
[ 0.11,0.13],
[-0.81,0.95],
[ 0.29,0.27]])

predresults=np.zeros(expresults.shape)

fp1=open("Corrected.txt","r")

lines1=fp1.readlines()

nsites=len(lines1[0].split())-3
i1=np.zeros((nsites,),dtype='int')
entry=np.zeros((2,))
for i in range(0,len(lines1)):
  line1=lines1[i].split()

  for j in range(0,nsites):
    i1[j]=int(line1[j])
  entry[0]=float(line1[nsites])
  entry[1]=float(line1[nsites+2])

  for j in range(len(exptags)):
    if np.array_equal(i1,expseq[j]):
      predresults[j,:]=entry

fp1.close()

print(np.concatenate((expresults,predresults),axis=1))
print(np.corrcoef(expresults[:,0],y=predresults[:,0],rowvar=False)[0,1])
d=predresults[:,0]-expresults[:,0]
print(np.sqrt(np.mean(d**2)-np.mean(d)**2))
