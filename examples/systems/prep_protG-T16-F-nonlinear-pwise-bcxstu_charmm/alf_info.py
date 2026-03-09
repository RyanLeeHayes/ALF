import numpy as np
alf_info={}
alf_info['name']='protG'
alf_info['nsubs']=np.loadtxt('prep/nsubs',dtype='int',ndmin=1)
alf_info['nblocks']=np.sum(alf_info['nsubs'])
alf_info['ncentral']=0
alf_info['nreps']=1
alf_info['nnodes']=1
alf_info['enginepath']='/data/homezvol0/rhayes1/CHARMM_EXE/gnu/charmm'
alf_info['temp']=298.15
alf_info['loss']='nonlinear2024'
alf_info['bias']='bcxstu2026'
alf_info['impcons']='fnpwise2026'
alf_info['impconsopt']={'width':0.35,'kbeta':100}
alf_info['q']=np.loadtxt('prep/q',ndmin=1)
