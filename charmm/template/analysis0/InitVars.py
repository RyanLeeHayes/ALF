#! /usr/bin/env python

import sys, os
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

if not os.path.exists('../nbshift'):
  os.mkdir('../nbshift')
np.savetxt('../nbshift/b_shift.dat',b)
np.savetxt('../nbshift/c_shift.dat',c)
np.savetxt('../nbshift/x_shift.dat',c)
np.savetxt('../nbshift/s_shift.dat',c)

call(['../ALF/SetVars.py','1'])
