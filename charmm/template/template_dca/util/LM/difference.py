#! /usr/bin/env python

import math

fp1=open("../../run_unfolded/ResultU.txt","r")
fp2=open("../../run/1_400ns/dca37/Result.txt","r")
fp3=open("Result.txt","w")

lines1=fp1.readlines()
lines2=fp2.readlines()

nsites=len(lines1[0].split())-3
for i in range(0,len(lines1)):
  line1=lines1[i].split()
  line2=lines2[i].split()

  i1=[]
  for j in range(0,nsites):
    i1.append(int(line1[j]))
  V=float(line2[nsites])-float(line1[nsites])
  E=math.sqrt(float(line2[nsites+2])**2 + float(line1[nsites+2])**2)

  for j in range(0,nsites):
    fp3.write("%2d " % (i1[j],))
  fp3.write("%8.3f +/- %5.3f\n" % (V,E))

fp1.close()
fp2.close()
fp3.close()
