#! /usr/bin/env python

import sys
import numpy as np
from subprocess import call

if len(sys.argv)<2:
  print("Error, missing step argument")
  quit()

Step=int(sys.argv[1])

nblocks=np.loadtxt('../nblocks',dtype='int')
nsubs=np.loadtxt('../nsubs',dtype='int',ndmin=1)
nreps=np.loadtxt('../nreps',dtype='int')
nnodes=np.loadtxt('../nnodes',dtype='int')
ncentral=np.loadtxt('../ncentral',dtype='int')

fp=open('../name','r')
name=fp.readline().strip()
fp.close()


b_prev=np.loadtxt('b_prev.dat')
b=np.loadtxt('b.dat')
b_sum=b_prev+b
b_sum=np.reshape(b_sum,(1,-1))
np.savetxt('b_sum.dat',b_sum,fmt=' %7.2f')

fp=open("fb_est.inp","w")
ibuff=0
for i in range(0,len(nsubs)):
  for j in range(0,nsubs[i]):
    line=("set lams%ds%d = %8.2f\n" % (i+1,j+1,b_sum[0,ibuff+j]))
    fp.write(line)
  ibuff+=nsubs[i]
fp.close()

c_prev=np.loadtxt('c_prev.dat')
c=np.loadtxt('c.dat')
c_sum=c_prev+c
np.savetxt('c_sum.dat',c_sum,fmt=' %7.2f')

fp=open("vb_est.inp","w")
ibuff=0
for si in range(0,len(nsubs)):
  jbuff=ibuff
  for sj in range(si,len(nsubs)):
    for i in range(0,nsubs[si]):
      ii=i+ibuff
      j0=0
      if si==sj:
        j0=i+1
      for j in range(j0,nsubs[sj]):
        jj=j+jbuff
        line=("set cs%ds%ds%ds%d = %8.2f\n" % (si+1,i+1,sj+1,j+1,-c_sum[ii,jj]))
        fp.write(line)
    jbuff+=nsubs[sj]
  ibuff+=nsubs[si]
fp.close()

x_prev=np.loadtxt('x_prev.dat')
x=np.loadtxt('x.dat')
x_sum=x_prev+x
np.savetxt('x_sum.dat',x_sum,fmt=' %7.2f')

fp=open("xb_est.inp","w")
ibuff=0
for si in range(0,len(nsubs)):
  jbuff=0
  for sj in range(0,len(nsubs)):
    for i in range(0,nsubs[si]):
      ii=i+ibuff
      for j in range(0,nsubs[sj]):
        jj=j+jbuff
        if ii!=jj:
          line=("set xs%ds%ds%ds%d = %8.2f\n" % (si+1,i+1,sj+1,j+1,-x_sum[ii,jj]))
          fp.write(line)
    jbuff+=nsubs[sj]
  ibuff+=nsubs[si]
fp.close()

s_prev=np.loadtxt('s_prev.dat')
s=np.loadtxt('s.dat')
s_sum=s_prev+s
np.savetxt('s_sum.dat',s_sum,fmt=' %7.2f')

fp=open("sb_est.inp","w")
ibuff=0
for si in range(0,len(nsubs)):
  jbuff=0
  for sj in range(0,len(nsubs)):
    for i in range(0,nsubs[si]):
      ii=i+ibuff
      for j in range(0,nsubs[sj]):
        jj=j+jbuff
        if ii!=jj:
          line=("set ss%ds%ds%ds%d = %8.2f\n" % (si+1,i+1,sj+1,j+1,-s_sum[ii,jj]))
          fp.write(line)
    jbuff+=nsubs[sj]
  ibuff+=nsubs[si]
fp.close()

fp=open("parm.inp","w")
fp.write("set sysname = \""+name+"\n")
fp.write("trim sysname from 2\n")
fp.write("set nnodes = "+str(nnodes)+"\n")
fp.write("set nreps = "+str(nreps)+"\n")
fp.write("set ncentral = "+str(ncentral)+"\n")
fp.write("set nblocks = "+str(nblocks)+"\n")
fp.write("set nsites = "+str(len(nsubs))+"\n")
for i in range(0,len(nsubs)):
  fp.write("set nsubs"+str(i+1)+" = "+str(nsubs[i])+"\n")
fp.close()

fp=open('../variables'+str(Step)+'.inp','w')
call(['cat','fb_est.inp','vb_est.inp','xb_est.inp','sb_est.inp','parm.inp'],stdout=fp)
fp.close()
