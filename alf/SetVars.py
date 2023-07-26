#! /usr/bin/env python

def SetVarsCharmm(alf_info,Step,minimize=False):
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
  fp.write("set minimize = "+str(int(minimize==True))+"\n")
  fp.write("\n")
  fp.close()



def SetVarsBlade(alf_info,Step,minimize=False):
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
  import numpy as np
  import yaml
  import copy

  fp=open('../variables'+str(Step)+'.py','w')

  fp.write("import yaml\n")
  fp.write("import numpy as np\n")

  bias={}

  b_prev=np.loadtxt('b_prev.dat')
  b=np.loadtxt('b.dat')
  b_sum=b_prev+b
  b_sum=np.reshape(b_sum,(1,-1))
  np.round(b_sum,decimals=2)
  np.savetxt('b_sum.dat',b_sum,fmt=' %7.2f')

  bias['b']=b_sum.tolist()

  c_prev=np.loadtxt('c_prev.dat')
  c=np.loadtxt('c.dat')
  c_sum=c_prev+c
  np.round(c_sum,decimals=2)
  np.savetxt('c_sum.dat',c_sum,fmt=' %7.2f')

  bias['c']=c_sum.tolist()

  x_prev=np.loadtxt('x_prev.dat')
  x=np.loadtxt('x.dat')
  x_sum=x_prev+x
  np.round(x_sum,decimals=2)
  np.savetxt('x_sum.dat',x_sum,fmt=' %7.2f')

  bias['x']=x_sum.tolist()

  s_prev=np.loadtxt('s_prev.dat')
  s=np.loadtxt('s.dat')
  s_sum=s_prev+s
  np.round(s_sum,decimals=2)
  np.savetxt('s_sum.dat',s_sum,fmt=' %7.2f')

  bias['s']=s_sum.tolist()

  fp.write("bias_string=\"\"\"\n")
  yaml.dump(bias,fp)
  fp.write("\"\"\"\n")
  fp.write("bias=yaml.load(bias_string,Loader=yaml.Loader)\n")
  fp.write("bias['b']=np.array(bias['b'])\n")
  fp.write("bias['c']=np.array(bias['c'])\n")
  fp.write("bias['x']=np.array(bias['x'])\n")
  fp.write("bias['s']=np.array(bias['s'])\n")

  nsubs=alf_info['nsubs']

  ibuff=0
  lamss={}
  for i in range(0,len(nsubs)):
    lamss[i]={}
    for j in range(0,nsubs[i]):
      lamss[i][j]=b_sum[0,ibuff+j].tolist()
    ibuff+=nsubs[i]

  ibuff=0
  cssss={}
  for si in range(0,len(nsubs)):
    jbuff=ibuff
    cssss[si]={}
    for sj in range(si,len(nsubs)):
      cssss[si][sj]={}
      for i in range(0,nsubs[si]):
        ii=i+ibuff
        cssss[si][sj][i]={}
        j0=0
        if si==sj:
          j0=i+1
        for j in range(j0,nsubs[sj]):
          jj=j+jbuff
          cssss[si][sj][i][j]=-c_sum[ii,jj].tolist()
      jbuff+=nsubs[sj]
    ibuff+=nsubs[si]

  ibuff=0
  xssss={}
  for si in range(0,len(nsubs)):
    xssss[si]={}
    jbuff=0
    for sj in range(0,len(nsubs)):
      xssss[si][sj]={}
      for i in range(0,nsubs[si]):
        ii=i+ibuff
        xssss[si][sj][i]={}
        for j in range(0,nsubs[sj]):
          jj=j+jbuff
          if ii!=jj:
            xssss[si][sj][i][j]=-x_sum[ii,jj].tolist()
      jbuff+=nsubs[sj]
    ibuff+=nsubs[si]

  ibuff=0
  sssss={}
  for si in range(0,len(nsubs)):
    sssss[si]={}
    jbuff=0
    for sj in range(0,len(nsubs)):
      sssss[si][sj]={}
      for i in range(0,nsubs[si]):
        ii=i+ibuff
        sssss[si][sj][i]={}
        for j in range(0,nsubs[sj]):
          jj=j+jbuff
          if ii!=jj:
            sssss[si][sj][i][j]=-s_sum[ii,jj].tolist()
      jbuff+=nsubs[sj]
    ibuff+=nsubs[si]

  fp.write("lamss_string=\"\"\"\n")
  yaml.dump(lamss,fp)
  fp.write("\"\"\"\n")
  fp.write("lamss=yaml.load(lamss_string,Loader=yaml.Loader)\n")
  fp.write("cssss_string=\"\"\"\n")
  yaml.dump(cssss,fp)
  fp.write("\"\"\"\n")
  fp.write("cssss=yaml.load(cssss_string,Loader=yaml.Loader)\n")
  fp.write("xssss_string=\"\"\"\n")
  yaml.dump(xssss,fp)
  fp.write("\"\"\"\n")
  fp.write("xssss=yaml.load(xssss_string,Loader=yaml.Loader)\n")
  fp.write("sssss_string=\"\"\"\n")
  yaml.dump(sssss,fp)
  fp.write("\"\"\"\n")
  fp.write("sssss=yaml.load(sssss_string,Loader=yaml.Loader)\n")

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

  if minimize==True:
    fp.write("minimizeflag=True\n")
  else:
    fp.write("minimizeflag=False\n")

  fp.close()



def SetVars(alf_info,Step,minimize=False):
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
