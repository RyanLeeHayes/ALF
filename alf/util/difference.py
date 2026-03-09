#! /usr/bin/env python

def result_difference(reference,target,out,referencename=None,targetname=None,outname=None):
  """
  Obtain ddG by taking difference between target minus reference.

  Parameters
  ----------
  reference : string
      The directory in which the reference Result.txt file is located
  target : string
      The directory in which the target Result.txt file is located
  out : string
      The directory path in which to place the output Result.txt
  referencename : string, optional
      Name of file in reference directory if other than Result.txt
  targetname : string, optional
      Name of file in target directory if other than Result.txt
  outname : string, optional
      File name for writing, if other than Result.txt
  """

  import math

  if referencename==None:
    fp1=open(reference+"/Result.txt","r")
  else:
    fp1=open(reference+"/"+referencename,"r")
  if targetname==None:
    fp2=open(target+"/Result.txt","r")
  else:
    fp2=open(target+"/"+targetname,"r")
  if outname==None:
    fp3=open(out+"/Result.txt","w")
  else:
    fp3=open(out+"/"+outname,"w")

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

