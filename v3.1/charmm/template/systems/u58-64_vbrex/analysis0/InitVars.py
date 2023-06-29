#! /usr/bin/env python

import sys
import numpy as np
from subprocess import call

nblocks=np.loadtxt('../nblocks',dtype='int')

b=np.zeros([1,nblocks])
c=np.zeros([nblocks,nblocks])

np.savetxt('b_prev.dat',b)
np.savetxt('b.dat',b)
np.savetxt('c_prev.dat',c)
np.savetxt('c.dat',c)
np.savetxt('x_prev.dat',c)
np.savetxt('x.dat',c)
np.savetxt('s_prev.dat',c)
np.savetxt('s.dat',c)

call(['../ALF/SetVars.py','1'])
