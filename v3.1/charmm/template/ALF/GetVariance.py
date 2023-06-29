#! /usr/bin/env python

import sys, os
import numpy as np
import copy

NF=int(sys.argv[1])

nblocks=np.loadtxt('../nblocks',dtype='int')
nsubs=np.loadtxt('../nsubs',dtype='int',ndmin=1)
nreps=np.loadtxt('../nreps',dtype='int')
nlig=np.prod(nsubs)

b=np.loadtxt('b_prev.dat')
c=np.loadtxt('c_prev.dat')
x=np.loadtxt('x_prev.dat')
s=np.loadtxt('s_prev.dat')

b_shift=np.loadtxt('../nbshift/b_shift.dat')
c_shift=np.loadtxt('../nbshift/c_shift.dat')
x_shift=np.loadtxt('../nbshift/x_shift.dat')
s_shift=np.loadtxt('../nbshift/s_shift.dat')
ncentral=np.loadtxt('../ncentral',dtype='int')

f=np.loadtxt('f.dat')

kT=0.001987*298

G=np.zeros((NF,nlig))
ind=np.zeros((nlig,len(nsubs)),dtype='int')
for i in range(1,nlig):
  ind[i,:]=ind[i-1,:]
  for j in range(len(nsubs)-1,-1,-1):
    if (ind[i,j]+1)<nsubs[j]:
      ind[i,j]+=1
      break
    else:
      ind[i,j]=0
blk=copy.deepcopy(ind)
for i in range(1,len(nsubs)):
  blk[:,i:]+=nsubs[i-1]

Eall=np.zeros((nreps,NF,nlig))
Eshift=np.zeros((nreps,NF,nlig))
lndenom=np.zeros((nreps,NF,nlig))
nframes=np.zeros((nreps,NF))
for irep in range(0,nreps):
  for i in range(0,NF):
    isim=i*nreps+irep
    L=np.loadtxt('data/Lambda.'+str(i)+'.'+str(irep)+'.dat')
    nframes[irep,i]=L.shape[0]
    for j in range(0,nlig):
      LList=np.zeros((nblocks,))
      LList[blk[j,:]]=1
      # E=np.sum(b[blk[j,:]])
      Eall[irep,i,j]=np.dot(LList,-b)+np.dot(np.dot(LList,-c),LList)+np.dot(np.dot(1-np.exp(-5.56*LList),-x),LList)+np.dot(np.dot(LList/(LList+0.017),-s),LList)
      Eshift[irep,i,j]=(irep-ncentral)*(np.dot(LList,-b_shift)+np.dot(np.dot(LList,-c_shift),LList)+np.dot(np.dot(1-np.exp(-5.56*LList),-x_shift),LList)+np.dot(np.dot(LList/(LList+0.017),-s_shift),LList))
      # lndenom[irep,i,j]=np.log(nframes[irep,i])+f[isim]-(Eall[irep,i,j]+Eshift[irep,i,j])/kT
      lndenom[irep,i,j]=np.log(nframes[irep,i])+f[isim]-(Eshift[irep,i,j])/kT

PkeepA=np.zeros((nreps,NF,nlig))
for irep in range(0,nreps):
  Pkeep=np.zeros((NF,nlig))
  G=np.zeros((NF,nlig))
  for i in range(0,NF):
    print(i)
    L=np.loadtxt('data/Lambda.'+str(i)+'.'+str(irep)+'.dat')
    for j in range(0,nlig):
      P=np.sum(np.all(L[:,blk[j,:]]>0.99,axis=1))
      Pkeep[i,j]=P
      G[i,j]=-Eall[irep,i,j]-Eshift[irep,i,j]-kT*np.log(P)
  PkeepA[irep,:,:]=Pkeep

  np.savetxt('G.'+str(irep)+'.dat',G)

  Gmin=np.min(G,axis=0)
  Value=Gmin-kT*np.log(np.mean(np.exp(-(G-Gmin)/kT),axis=0))
  Value-=Value[0]

  np.random.seed(2401)
  GS=np.zeros((50,nlig))
  for i in range(0,50):
    GS[i,:]=Gmin-kT*np.log(np.mean(np.exp(-(G[np.random.randint(0,NF,(NF,1)),:]-Gmin)/kT),axis=0))
  Error=np.std(GS,axis=0)

  fp=open('Result.'+str(irep)+'.txt','w')

  for i in range(0,nlig):
    for j in range(0,len(nsubs)):
      fp.write('%2d ' % ind[i,j])
    fp.write('%8.3f +/- %5.3f\n' % (Value[i],Error[i]))

  fp.close()

# WORKING
print(PkeepA)
print(Eall)
print(Eshift)
print(lndenom)
# Pkeep=np.sum(PkeepA,axis=0)
for i in range(0,NF):
  for j in range(0,nlig):
    # Actual expressions
    # G[i,j]=Eall[ncentral,i,j]-kT*np.log(np.sum(PkeepA[:,i,j]*np.exp(-Eshift[ncentral,i,j]/kT)/np.sum(np.exp(lndenom[:,i,j]))))
    # More efficient? - exp(-Eshift(central)/kT)=1`
    G[i,j]=-Eall[ncentral,i,j]-kT*np.log(np.sum(PkeepA[:,i,j]))+kT*np.log(np.sum(np.exp(lndenom[:,i,j])))

np.savetxt('G.dat',G)

Gmin=np.min(G,axis=0)
Value=Gmin-kT*np.log(np.mean(np.exp(-(G-Gmin)/kT),axis=0))
Value-=Value[0]

np.random.seed(2401)
GS=np.zeros((50,nlig))
for i in range(0,50):
  GS[i,:]=Gmin-kT*np.log(np.mean(np.exp(-(G[np.random.randint(0,NF,(NF,1)),:]-Gmin)/kT),axis=0))
Error=np.std(GS,axis=0)

fp=open('Result.txt','w')

for i in range(0,nlig):
  for j in range(0,len(nsubs)):
    fp.write('%2d ' % ind[i,j])
  fp.write('%8.3f +/- %5.3f\n' % (Value[i],Error[i]))

fp.close()
