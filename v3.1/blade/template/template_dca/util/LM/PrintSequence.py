#! /usr/bin/env python

import sys
import math

conv321={'ALA':'A','HSE':'B','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','HSD':'H','ILE':'I','HSP':'J','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}

fp_pdb=open("../../run_ms/1_pH5.5/prep/minimized.pdb","r")

seq0=''
for line in fp_pdb:
  if len(line)==77:
    segid=line[72:76]
    if segid=='PROT':
      resname=line[16:21].strip()
      resn=conv321[resname]
      resid=int(line[21:26])
      if resid>len(seq0)+1:
        print('Error, resid '+str(resid)+' resname '+resname+' is beyond end of sequence '+seq0)
        quit()
      else:
        seq0=seq0[0:resid-1]+resn+seq0[resid:]
      # print(line)
      # print(seq0)

# fp_setup=open("../../run_ms/1_pH5.5/prep/ec2acc.inp","r")
fp_setup=open("mut.dat","r")

mut_resid=[]
mut_resn=[]
for line in fp_setup:
  mut_resid.append(int(line.split()[0]))
  mut_resn.append([])
  for resn in line.split()[1:]:
    if resn=='0':
      mut_resn[-1].append(seq0[mut_resid[-1]-1])
    else:
      mut_resn[-1].append(resn.upper())
print(mut_resid)
print(mut_resn)

if len(sys.argv)<=len(mut_resid):
  print('Error, need '+str(len(mut_resid))+' integer arguments')
  quit()
mut=[]
for i in range(len(mut_resid)):
  try:
    mut.append(int(sys.argv[i+1]))
  except:
    print('Error converting argument '+str(i+1))
    quit()
  if mut[i]<0 or mut[i]>=len(mut_resn[i]):
    print('Error converting argument '+str(i+1)+', value is outside bounds')
    quit()

seq=seq0[0:]
for i in range(len(mut_resid)):
  seq=seq[0:mut_resid[i]-1]+mut_resn[i][mut[i]]+seq[mut_resid[i]:]
  # print(seq)

for i in range(len(seq)):
  if seq[i]=='B':
    seq=seq[0:i]+'H'+seq[i+1:]
    # print(seq)
  elif seq[i]=='J':
    seq=seq[0:i]+'H'+seq[i+1:]
    # print(seq)
  elif seq[i]=='_':
    seq=seq[0:i]+seq[i+1:]
    # print(seq)
print(seq)
