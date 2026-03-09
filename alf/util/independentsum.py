#! /usr/bin/env python

def result_independentsum(inlist,out,inlistname=None,outname=None):
  """
  Concatenate Result.txt files in different directories

  Provides a Result.txt file in the out directory by stacking together all
  possible permutations of Result.txt files in the inlist directories.
  Uses for combining results from different pentapeptides into a single
  model of the unfolded state or combining protein and ligand results to
  determine the free energy of the unbound state upon both protein and
  ligand mutation.

  Parameters
  ----------
  inlist : list of strings
      The directory paths to be concatenated together
  out : string
      The directory path in which to place the output Result.txt
  inlistname : list of strings, optional
      List of file names in directory if other than Result.txt
  outname : string, optional
      File name for writing, if other than Result.txt
  """

  import math

  fpin=[]
  for id in range(len(inlist)):
    if inlistname==None:
      fpin.append(open(inlist[id]+"/Result.txt","r"))
    else:
      fpin.append(open(inlist[id]+"/"+inlistname[id],"r"))
  if outname==None:
    fpout=open(out+"/Result.txt","w")
  else:
    fpout=open(out+"/"+outname,"w")

  i1=[]
  V=[]
  E=[]
  indices=[]
  nsubs=[]
  nsites=[]
  for i in range(0,len(fpin)):
    i1.append([])
    V.append([])
    E.append([])
    indices.append(0)
    lines=fpin[i].readlines()
    nsites.append(len(lines[0].split())-3)
    for j in range(0,len(lines)):
      line=lines[j].split()

      i1[i].append([])
      for k in range(0,nsites[i]):
        i1[i][j].append(int(line[k]))
      V[i].append(float(line[nsites[i]]))
      E[i].append(float(line[nsites[i]+2]))
    for j in range(0,nsites[i]):
      nsubs.append(i1[i][-1][j]+1)
  

  while indices[0]<len(i1[0]):
    Vi=0
    Ei=0
    for i in range(0,len(fpin)):
      for j in range(0,nsites[i]):
        fpout.write("%2d " % (i1[i][indices[i]][j],))
      Vi+=V[i][indices[i]]
      Ei=math.sqrt(Ei**2 + E[i][indices[i]]**2)
    fpout.write("%8.3f +/- %5.3f\n" % (Vi,Ei))
    for i in range(0,len(fpin)):
      indices[len(fpin)-1-i]+=1
      if indices[len(fpin)-1-i]==len(i1[len(fpin)-1-i]) and (len(fpin)-1-i)!=0:
        indices[len(fpin)-1-i]=0
      else:
        break

  for i in range(0,len(fpin)):
    fpin[i].close()
  fpout.close()

