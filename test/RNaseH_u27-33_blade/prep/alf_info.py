import numpy as np
alf_info={}
alf_info['name']='u27-33'
alf_info['nsubs']=np.loadtxt('prep/nsubs',dtype='int',ndmin=1)
alf_info['nblocks']=np.sum(alf_info['nsubs'])
alf_info['ncentral']=0
alf_info['nreps']=1
alf_info['nnodes']=1
alf_info['enginepath']='/data/homezvol0/rhayes1/blade/build/blade'
alf_info['temp']=298.15
