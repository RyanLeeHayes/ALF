#! /usr/bin/env python

import sys
import numpy as np

import MDAnalysis as mda
# See https://docs.python.org/2/tutorial/modules.html for instructions on submodules:
# import MDAnalysis.analysis.distances
# from MDAnalysis.analysis import contacts

# See alse https://www.mdanalysis.org/docs/documentation_pages/coordinates/init.html

if len(sys.argv) < 4:
  print("Error: need an input, psf, and output filename")
  quit()

Volumes=np.zeros((0,1))

for ixtc in range(3,len(sys.argv)):
  trajin=mda.Universe(sys.argv[2],format="PSF")
  trajin.load_new(sys.argv[ixtc],format="XTC")

  lines=trajin.trajectory.n_frames
  Volume=np.zeros((lines,1))

  i=0
  for ts in trajin.trajectory:
    if (ts.dimensions[3]!=90 or ts.dimensions[4]!=90 or ts.dimensions[5]!=90):
      print("GetVolume.py rquires an orthorhombic box")
    Volume[i]=ts.dimensions[0]*ts.dimensions[1]*ts.dimensions[2]
    i+=1
  Volumes=np.concatenate((Volumes,Volume),axis=0)

selection=trajin.select_atoms("resname TIP3 and name OH2")

np.savetxt(sys.argv[1],Volumes,fmt="%10.6f")
np.savetxt(sys.argv[1]+"_NH2O",np.reshape(selection.n_atoms,(1,)),fmt="%d")
