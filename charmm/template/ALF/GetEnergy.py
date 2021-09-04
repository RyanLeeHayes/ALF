#! /usr/bin/env python

import sys, os
import numpy as np

Fi=int(sys.argv[1])
Ff=int(sys.argv[2])
NF=Ff-Fi+1;

if len(sys.argv)>3:
  skipE=int(sys.argv[3])
else:
  skipE=1

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

nblocks=np.loadtxt('../nblocks',dtype='int')
nsubs=np.loadtxt('../nsubs',dtype='int',ndmin=1)
nreps=np.loadtxt('../nreps',dtype='int')
ncentral=np.loadtxt('../ncentral',dtype='int')

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
