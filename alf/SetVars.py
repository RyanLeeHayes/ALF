#! /usr/bin/env python

def SetVarsCharmm(alf_info,Step,minimize=False):
  """
  Writes out a variables[Step].inp file to read biases into MD engine

  For the charmm (and bladelib) engine. Called by alf.SetVars.

  Parameters
  ----------
  alf_info : dict
      Dictionary of variables alf needs to run. alf_info['engine']
      determines the format of the variables[Step].inp file
  Step : int
      The next cycle of alf for which the new biases are being written
  minimize : bool, optional
      A boolean flag indicating whether or not to run minimization on this
      cycle of alf. (default is False)
  """

  import numpy as np

  nblocks=alf_info['nblocks']
  nsubs=alf_info['nsubs']
  nreps=alf_info['nreps']
  ncentral=alf_info['ncentral']
  name=alf_info['name']
  nnodes=alf_info['nnodes']
  temp=alf_info['temp']

  fp=open('../variables'+str(Step)+'.inp','w')
  fp.write("* Variables from step %d of ALF\n*\n\n" % (Step,))

  b_prev=np.loadtxt('b_prev.dat')
  b=np.loadtxt('b.dat')
  b_sum=b_prev+b
  b_sum=np.reshape(b_sum,(1,-1))
  np.savetxt('b_sum.dat',b_sum,fmt=' %7.2f')

  ibuff=0
  for i in range(0,len(nsubs)):
    for j in range(0,nsubs[i]):
      line=("set lams%ds%d = %8.2f\n" % (i+1,j+1,b_sum[0,ibuff+j]))
      fp.write(line)
    ibuff+=nsubs[i]

  c_prev=np.loadtxt('c_prev.dat')
  c=np.loadtxt('c.dat')
  c_sum=c_prev+c
  np.savetxt('c_sum.dat',c_sum,fmt=' %7.2f')

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

  x_prev=np.loadtxt('x_prev.dat')
  x=np.loadtxt('x.dat')
  x_sum=x_prev+x
  np.savetxt('x_sum.dat',x_sum,fmt=' %7.2f')

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

  s_prev=np.loadtxt('s_prev.dat')
  s=np.loadtxt('s.dat')
  s_sum=s_prev+s
  np.savetxt('s_sum.dat',s_sum,fmt=' %7.2f')

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

  fp.write("set sysname = \""+name+"\n")
  fp.write("trim sysname from 2\n")
  fp.write("set nnodes = "+str(nnodes)+"\n")
  fp.write("set nreps = "+str(nreps)+"\n")
  fp.write("set ncentral = "+str(ncentral)+"\n")
  fp.write("set nblocks = "+str(nblocks)+"\n")
  fp.write("set nsites = "+str(len(nsubs))+"\n")
  for i in range(0,len(nsubs)):
    fp.write("set nsubs"+str(i+1)+" = "+str(nsubs[i])+"\n")
  fp.write("set temp = "+str(temp)+"\n")
  fp.write("set minimizeflag = "+str(int(minimize==True))+"\n")
  fp.write("\n")
  fp.close()



def SetVarsBlade(alf_info,Step,minimize=False):
  """
  Writes out a variables[Step].inp file to read biases into MD engine

  For the blade engine. Called by alf.SetVars.

  Parameters
  ----------
  alf_info : dict
      Dictionary of variables alf needs to run. alf_info['engine']
      determines the format of the variables[Step].inp file
  Step : int
      The next cycle of alf for which the new biases are being written
  minimize : bool, optional
      A boolean flag indicating whether or not to run minimization on this
      cycle of alf. (default is False)
  """

  import numpy as np

  nblocks=alf_info['nblocks']
  nsubs=alf_info['nsubs']
  nreps=alf_info['nreps']
  ncentral=alf_info['ncentral']
  name=alf_info['name']
  nnodes=alf_info['nnodes']
  temp=alf_info['temp']

  fp=open('../variables'+str(Step)+'.inp','w')

  b_prev=np.loadtxt('b_prev.dat')
  b=np.loadtxt('b.dat')
  b_sum=b_prev+b
  b_sum=np.reshape(b_sum,(1,-1))
  np.savetxt('b_sum.dat',b_sum,fmt=' %7.2f')

  ibuff=0
  for i in range(0,len(nsubs)):
    for j in range(0,nsubs[i]):
      line=("variables set lams%ds%d %8.2f\n" % (i+1,j+1,-b_sum[0,ibuff+j]))
      fp.write(line)
    ibuff+=nsubs[i]

  c_prev=np.loadtxt('c_prev.dat')
  c=np.loadtxt('c.dat')
  c_sum=c_prev+c
  np.savetxt('c_sum.dat',c_sum,fmt=' %7.2f')

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
          line=("variables set cs%ds%ds%ds%d %8.2f\n" % (si+1,i+1,sj+1,j+1,-c_sum[ii,jj]))
          fp.write(line)
      jbuff+=nsubs[sj]
    ibuff+=nsubs[si]

  x_prev=np.loadtxt('x_prev.dat')
  x=np.loadtxt('x.dat')
  x_sum=x_prev+x
  np.savetxt('x_sum.dat',x_sum,fmt=' %7.2f')

  ibuff=0
  for si in range(0,len(nsubs)):
    jbuff=0
    for sj in range(0,len(nsubs)):
      for i in range(0,nsubs[si]):
        ii=i+ibuff
        for j in range(0,nsubs[sj]):
          jj=j+jbuff
          if ii!=jj:
            line=("variables set xs%ds%ds%ds%d %8.2f\n" % (si+1,i+1,sj+1,j+1,-x_sum[ii,jj]))
            fp.write(line)
      jbuff+=nsubs[sj]
    ibuff+=nsubs[si]

  s_prev=np.loadtxt('s_prev.dat')
  s=np.loadtxt('s.dat')
  s_sum=s_prev+s
  np.savetxt('s_sum.dat',s_sum,fmt=' %7.2f')

  ibuff=0
  for si in range(0,len(nsubs)):
    jbuff=0
    for sj in range(0,len(nsubs)):
      for i in range(0,nsubs[si]):
        ii=i+ibuff
        for j in range(0,nsubs[sj]):
          jj=j+jbuff
          if ii!=jj:
            line=("variables set ss%ds%ds%ds%d %8.2f\n" % (si+1,i+1,sj+1,j+1,-s_sum[ii,jj]))
            fp.write(line)
      jbuff+=nsubs[sj]
    ibuff+=nsubs[si]

  fp.write("variables set sysname "+name+"\n")
  fp.write("variables set nnodes "+str(nnodes)+"\n")
  fp.write("variables set nreps "+str(nreps)+"\n")
  fp.write("variables set ncentral "+str(ncentral)+"\n")
  fp.write("variables set nblocks "+str(nblocks)+"\n")
  fp.write("variables set nsites "+str(len(nsubs))+"\n")
  for i in range(0,len(nsubs)):
    fp.write("variables set nsubs"+str(i+1)+" "+str(nsubs[i])+"\n")
  fp.write("variables set temp "+str(temp)+"\n")
  if minimize:
    print("Warning: Implicitly or explicitly requested minimization with blade. Blade does not minimize. Set optional argument minimize=False in alf.initialize()\n")
  fp.close()



def SetVarsPycharmm(alf_info,Step,minimize=False):
  """
  Writes out a variables[Step].inp file to read biases into MD engine

  For the pycharmm engine. Called by alf.SetVars.

  Parameters
  ----------
  alf_info : dict
      Dictionary of variables alf needs to run. alf_info['engine']
      determines the format of the variables[Step].inp file
  Step : int
      The next cycle of alf for which the new biases are being written
  minimize : bool, optional
      A boolean flag indicating whether or not to run minimization on this
      cycle of alf. (default is False)
  """

  import numpy as np
  import yaml
  import copy

  nblocks=alf_info['nblocks']
  nsubs=alf_info['nsubs']
  nreps=alf_info['nreps']
  ncentral=alf_info['ncentral']
  name=alf_info['name']
  nnodes=alf_info['nnodes']
  temp=alf_info['temp']

  fp=open('../variables'+str(Step)+'.py','w')

  fp.write("import yaml\n")
  fp.write("import numpy as np\n")

  bias={}

  sub0=np.cumsum(nsubs)-nsubs

  b_prev=np.loadtxt('b_prev.dat')
  b=np.loadtxt('b.dat')
  b_sum=b_prev+b
  b_sum=np.reshape(b_sum,(1,-1))
  np.round(b_sum,decimals=2)
  np.savetxt('b_sum.dat',b_sum,fmt=' %7.2f')

  for i in range(0,len(nsubs)):
    for j in range(0,nsubs[i]):
      key=f'lams{i+1}s{j+1}'
      bias[key]=b_sum[0,sub0[i]+j].tolist()

  c_prev=np.loadtxt('c_prev.dat')
  c=np.loadtxt('c.dat')
  c_sum=c_prev+c
  np.round(c_sum,decimals=2)
  np.savetxt('c_sum.dat',c_sum,fmt=' %7.2f')

  for si in range(0,len(nsubs)):
    for sj in range(si,len(nsubs)):
      for i in range(0,nsubs[si]):
        j0=(i+1 if si==sj else 0)
        for j in range(j0,nsubs[sj]):
          key=f'cs{si+1}s{i+1}s{sj+1}s{j+1}'
          bias[key]=-c_sum[sub0[si]+i,sub0[sj]+j].tolist()

  x_prev=np.loadtxt('x_prev.dat')
  x=np.loadtxt('x.dat')
  x_sum=x_prev+x
  np.round(x_sum,decimals=2)
  np.savetxt('x_sum.dat',x_sum,fmt=' %7.2f')

  for si in range(0,len(nsubs)):
    for sj in range(0,len(nsubs)):
      for i in range(0,nsubs[si]):
        for j in range(0,nsubs[sj]):
          if sub0[si]+i!=sub0[sj]+j:
            key=f'xs{si+1}s{i+1}s{sj+1}s{j+1}'
            bias[key]=-x_sum[sub0[si]+i,sub0[sj]+j].tolist()

  s_prev=np.loadtxt('s_prev.dat')
  s=np.loadtxt('s.dat')
  s_sum=s_prev+s
  np.round(s_sum,decimals=2)
  np.savetxt('s_sum.dat',s_sum,fmt=' %7.2f')

  for si in range(0,len(nsubs)):
    for sj in range(0,len(nsubs)):
      for i in range(0,nsubs[si]):
        for j in range(0,nsubs[sj]):
          if sub0[si]+i!=sub0[sj]+j:
            key=f'ss{si+1}s{i+1}s{sj+1}s{j+1}'
            bias[key]=-s_sum[sub0[si]+i,sub0[sj]+j].tolist()

  bias['b']=b_sum.tolist()
  bias['c']=c_sum.tolist()
  bias['x']=x_sum.tolist()
  bias['s']=s_sum.tolist()

  fp.write("bias_string=\"\"\"\n")
  yaml.dump(bias,fp)
  fp.write("\"\"\"\n")
  fp.write("bias=yaml.load(bias_string,Loader=yaml.Loader)\n")
  fp.write("bias['b']=np.array(bias['b'])\n")
  fp.write("bias['c']=np.array(bias['c'])\n")
  fp.write("bias['x']=np.array(bias['x'])\n")
  fp.write("bias['s']=np.array(bias['s'])\n")

  alf_info_copy=copy.deepcopy(alf_info)
  alf_info_copy['nsubs']=alf_info_copy['nsubs'].tolist()
  alf_info_copy['nblocks']=alf_info_copy['nblocks'].tolist()
  if 'q' in alf_info_copy:
    alf_info_copy['q']=alf_info_copy['q'].tolist()

  fp.write("alf_info_string=\"\"\"\n")
  yaml.dump(alf_info_copy,fp)
  fp.write("\"\"\"\n")
  fp.write("alf_info=yaml.load(alf_info_string,Loader=yaml.Loader)\n")
  fp.write("alf_info['nsubs']=np.array(alf_info['nsubs'])\n")
  fp.write("if 'q' in alf_info:\n")
  fp.write("  alf_info['q']=np.array(alf_info['q'])\n")

  fp.write("sysname='"+name+"'\n")
  fp.write("nnodes="+str(nnodes)+"\n")
  fp.write("nreps="+str(nreps)+"\n")
  fp.write("ncentral="+str(ncentral)+"\n")
  fp.write("nblocks="+str(nblocks)+"\n")
  fp.write("nsites="+str(len(nsubs))+"\n")
  fp.write("nsubs="+str(nsubs.tolist())+"\n")
  fp.write("nsubs=np.array(nsubs)\n")
  for i in range(0,len(nsubs)):
    fp.write("nsubs"+str(i+1)+"="+str(nsubs[i])+"\n")
  fp.write("temp="+str(temp)+"\n")

  if minimize==True:
    fp.write("minimizeflag=True\n")
  else:
    fp.write("minimizeflag=False\n")

  fp.close()



def SetVars(alf_info,Step,minimize=False):
  """
  Writes out a variables[Step].inp file to read biases into MD engine

  Creates a variables[Step].inp file containing bias parameters, for the
  molecular dynamics engine to read in on the [Step] cycle of alf. This
  routine should be called from analysis[Step-1]. This is a wrapper
  routine, and based on the value of alf_info['engine'], either
  SetVarsCharmm, SetVarsBlade, or SetVarsPycharmm will be called.

  Parameters
  ----------
  alf_info : dict
      Dictionary of variables alf needs to run. alf_info['engine']
      determines the format of the variables[Step].inp file
  Step : int
      The next cycle of alf for which the new biases are being written
  minimize : bool, optional
      A boolean flag indicating whether or not to run minimization on this
      cycle of alf. (default is False)
  """

  if alf_info['engine'] in ['charmm','bladelib']:
    SetVarsCharmm(alf_info,Step,minimize=minimize)
  elif alf_info['engine'] in ['blade']:
    SetVarsBlade(alf_info,Step,minimize=minimize)
  elif alf_info['engine'] in ['pycharmm']:
    SetVarsPycharmm(alf_info,Step,minimize=minimize)
  else:
    print("Error: unsupported engine type %s" % alf_info['engine'])
    quit()



def InitVars(alf_info,minimize=True):
  """
  Prepares directory for first cycle of ALF

  Creates the analysis0 directory containing initial bias parameters, as
  well as nbshift for replica exchange. Calls alf.SetVars to create
  variables1.inp to input these bias parameters to the molecular dynamics
  engine. Called by alf.initialize.

  Parameters
  ----------
  alf_info : dict
      Dictionary of variables alf needs to run
  minimize : bool, optional
      A boolean flag indicating whether or not to run minimization on the
      first cycle of alf. (default is True)
  """

  import sys, os
  import numpy as np
  from subprocess import call

  nblocks=alf_info['nblocks']

  b=np.zeros([1,nblocks])
  c=np.zeros([nblocks,nblocks])

  if not os.path.exists('analysis0'):
    os.mkdir('analysis0')
  np.savetxt('analysis0/b_prev.dat',b)
  np.savetxt('analysis0/b.dat',b)
  np.savetxt('analysis0/c_prev.dat',c)
  np.savetxt('analysis0/c.dat',c)
  np.savetxt('analysis0/x_prev.dat',c)
  np.savetxt('analysis0/x.dat',c)
  np.savetxt('analysis0/s_prev.dat',c)
  np.savetxt('analysis0/s.dat',c)

  # Don't change nbshift if it exists
  if not os.path.exists('nbshift'):
    os.mkdir('nbshift')
    np.savetxt('nbshift/b_shift.dat',b)
    np.savetxt('nbshift/c_shift.dat',c)
    np.savetxt('nbshift/x_shift.dat',c)
    np.savetxt('nbshift/s_shift.dat',c)

  os.chdir('analysis0')
  SetVars(alf_info,1,minimize=minimize)
