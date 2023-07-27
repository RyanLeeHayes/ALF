##
## pyCHARMM script drafted from many examples
##
##   ** MSLD with BLADE (no OMM support) **
##

import os
import sys
import numpy as np
import pandas

##############################################
# Load pyCHARMM libraries

import pycharmm
import pycharmm.generate as gen
import pycharmm.ic as ic
import pycharmm.coor as coor
import pycharmm.energy as energy
import pycharmm.dynamics as dyn
import pycharmm.nbonds as nbonds
import pycharmm.minimize as minimize
import pycharmm.crystal as crystal
import pycharmm.image as image
import pycharmm.psf as psf
import pycharmm.read as read
import pycharmm.write as write
import pycharmm.settings as settings
import pycharmm.cons_harm as cons_harm
import pycharmm.cons_fix as cons_fix
import pycharmm.select as select
import pycharmm.shake as shake
import pycharmm.scalar as scalar
from pycharmm.lib import charmm as libcharmm


##############################################
# Set up global parameters

# nonbonded conditions
nb_fswitch = False          # normal fswitching functions
nb_pme = True               # normal PME

# dynamics conditions                
blade = True

# dynamics variables
cpt_on   = True             # run with CPT for NPT?
timestep = 0.002            # ps
# ns1      = 500000           # number of MD steps per 1 ns
# total_ns = 1                # total number of production ns sampling
# nequil = int(ns1*(1/10))    # equil for 100 ps
# nprod  = ns1*total_ns       # prod sampling for 5 ns
nsavc  = 1000               # dcd save frequency

##############################################
# Stream variables and system setup files
exec(open('arguments.py').read())
exec(open('variablesflat.py').read())
# sysname=alf_info['name'] # Set by variablesflat.py
exec(open('prep/'+sysname+'.py').read())


##############################################
# Set NonBonded settings & SP energy calc
cutnb = 14.0
cutim = cutnb
ctofnb = 12.0
ctonnb = 10.0

## nbond switching
## use a dictionary so that it becomes easy to switch between w/ vs w/o PME
nbonds_dict = {'cutnb':cutnb,'cutim':cutim,
           'ctonnb':ctonnb,'ctofnb':ctofnb,
           'atom':True,'vatom':True,
           'cdie':True,'eps':1.0,
           'inbfrq':-1, 'imgfrq':-1}

if nb_pme:
    nbonds_dict['switch']=True
    nbonds_dict['vfswitch']=True
    nbonds_dict['ewald']=True
    nbonds_dict['pmewald']=True
    nbonds_dict['kappa']=0.32
    nbonds_dict['fftx']=pmegrid
    nbonds_dict['ffty']=pmegrid
    nbonds_dict['fftz']=pmegrid
    nbonds_dict['order']=6

elif nb_fswitch:
    nbonds_dict['fswitch']=True
    nbonds_dict['vfswitch']=True
    nbonds_dict['ewald']=False
    nbonds_dict['pmewald']=False

else: 
    print("NonBonded Parameter Error - both pme and switch are false")
    pycharmm.lingo.charmm_script('stop')

nbonds=pycharmm.NonBondedScript(**nbonds_dict)
nbonds.run()
energy.show()

##############################################
# Minimize the system

if minimizeflag==True:
  minimize.run_sd(nstep=250,nprint=50,step=0.005,tolenr=1e-3,tolgrd=1e-3)
  energy.show()
  #minimize.run_abnr(nstep=250,nprint=50,tolenr=1e-3,tolgrd=1e-3)
  #energy.show()

  # write out psf, crd, pdb files
  write.psf_card('prep/minimized.psf')
  write.coor_card('prep/minimized.crd')
  write.coor_pdb('prep/minimized.pdb')

## Read in minimized coordinates
read.coor_card('prep/minimized.crd',resid=True)


##############################################
# Set up and run Dynamics

# dynamics conditions
if blade:
    useblade = 'prmc pref 1 iprs 100 prdv 100'
    gscal = 0.1
    ntrfrq=0
    leap = True
    openmm = False
else: 
    print("MSLD can only be run with BLADE - exiting...")
    pycharmm.lingo.charmm_script('stop')

# set shake
shake.on(bonh=True,fast=True,tol=1e-7)
dyn.set_fbetas(np.full((psf.get_natom()),gscal,dtype=float))

# initialize blade
pycharmm.lingo.charmm_script('energy blade')

# set up output directories
if not os.path.isdir('res'): os.system('mkdir res')
if not os.path.isdir('dcd'): os.system('mkdir dcd')

# set up dynamics dictionary of parameters
dynamics_dict = {'cpt':cpt_on,'leap':True,'langevin':True,
    'timestep':timestep,
    'nsavc':nsavc,
    'nsavl':10,  # frequency for saving lambda values in lamda-dynamics
    'nprint': 1000, # Frequency to write to output
    'iprfrq': 1000, # Frequency to calculate averages
    'ntrfrq':ntrfrq,
    'firstt':temp,'finalt':temp,'tstruct':temp,'tbath':temp,
    'iasors': 1,'iasvel':1,'iscvel': 0,'iscale': 0,
    'ihtfrq':0,'ieqfrq':0,'ichecw': 0,
    'inbfrq':-1,'imgfrq':-1,'ihbfrq':0,'ilbfrq':0,
    'echeck': -1}


if cpt_on:
    dynamics_dict['pconstant'] = True
    dynamics_dict['pmass'] = psf.get_natom()*0.12
    dynamics_dict['pref'] = 1.0
    dynamics_dict['pgamma'] = 20.0
    dynamics_dict['hoover'] = True
    dynamics_dict['reft'] = temp
    dynamics_dict['tmass'] = 1000


if blade:
    dynamics_dict['omm'] = False
    dynamics_dict['blade'] = useblade


# MD equilibration
dcd_file = pycharmm.CharmmFile(file_name='dcd/{}_{}.dcd'.format(sysname,'heat'), 
               file_unit=1,formatted=False,read_only=False)
res_file = pycharmm.CharmmFile(file_name='res/{}_{}.res'.format(sysname,'heat'), 
               file_unit=2,formatted=True,read_only=False)
lam_file = pycharmm.CharmmFile(file_name='res/{}_{}.lmd'.format(sysname,'heat'), 
               file_unit=3,formatted=False,read_only=False)

if not 'restartfile' in locals():
    dynamics_dict['start']  = True
    dynamics_dict['restart']= False
    dynamics_dict['iunrea'] = -1
else:
    prv_rest = pycharmm.CharmmFile(file_name=restartfile, 
                   file_unit=4,formatted=True,read_only=False)
    dynamics_dict['start']  = False
    dynamics_dict['restart']= True
    dynamics_dict['iunrea'] = prv_rest.file_unit
dynamics_dict['nstep']  = esteps
dynamics_dict['isvfrq'] = esteps # Frequency to save restart file
dynamics_dict['iunwri'] = res_file.file_unit
dynamics_dict['iuncrd'] = dcd_file.file_unit
dynamics_dict['iunldm'] = lam_file.file_unit

equil_dyn = pycharmm.DynamicsScript(**dynamics_dict)
equil_dyn.run()

if 'restartfile' in locals():
    prv_rest.close()
dcd_file.close()
res_file.close()
lam_file.close()

write.coor_pdb('dcd/{}_fframe.{}.pdb'.format(sysname,'equil')) # write out final frame


# MD production
dcd_file = pycharmm.CharmmFile(file_name='dcd/{}_{}.dcd'.format(sysname,'flat'), 
               file_unit=1,formatted=False,read_only=False)
res_file = pycharmm.CharmmFile(file_name='res/{}_{}.res'.format(sysname,'flat'), 
               file_unit=2,formatted=True,read_only=False)
prv_rest = pycharmm.CharmmFile(file_name='res/{}_{}.res'.format(sysname,'heat'), 
               file_unit=4,formatted=True,read_only=False)
lam_file = pycharmm.CharmmFile(file_name='res/{}_{}.lmd'.format(sysname,'flat'), 
               file_unit=3,formatted=False,read_only=False)

dynamics_dict['start']  = False
dynamics_dict['restart']= True
dynamics_dict['nstep']  = nsteps
dynamics_dict['isvfrq'] = nsteps # Frequency to save restart file
dynamics_dict['iunrea'] = prv_rest.file_unit
dynamics_dict['iunwri'] = res_file.file_unit
dynamics_dict['iuncrd'] = dcd_file.file_unit
dynamics_dict['iunldm'] = lam_file.file_unit

prod_dyn = pycharmm.DynamicsScript(**dynamics_dict)
prod_dyn.run()

dcd_file.close()
res_file.close()
lam_file.close()

# write.coor_pdb('dcd/{}_fframe.{}.pdb'.format(sysname,'flat')) # write out final frame

#if openmm: pycharmm.lingo.charmm_script('omm clear')
#if blade: pycharmm.lingo.charmm_script('blade off')

# # collect lambda statistics
# proc_lam = pycharmm.CharmmFile(file_name='res/{}_{}.lmd'.format(sysname,'flat'), 
#            file_unit=33,formatted=False,read_only=False)
# pycharmm.lingo.charmm_script('traj lamb print ctlo 0.95 cthi 0.99 first {} nunit {}'.format(proc_lam.file_unit,1))


##############################################
# FINISHED

pycharmm.lingo.charmm_script('stop')


