#! /usr/bin/env python

from MouseTrap import *

from obspy.core import AttribDict
paz = AttribDict({'sensitivity': 600000000.0, 'poles': [(-0.1103+0.111j), (-0.1103-0.111j), (-86.3+0j), (-241+178j), (-241-178j), (-535+719j), (-535-719j)], 'gain': 110400.0, 'zeros': [0j, 0j, (-68.8+0j), (-323+0j), (-2530+0j)]})

m1 = mouse()
m1.create(paz, 4000, 0.1, 100)

import matplotlib.pyplot as plt
t = np.arange(0, len(m1.mouse)/10., 0.1)
plt.plot(t,m1.mouse)
plt.show()
