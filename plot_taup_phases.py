#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tool to plot phases

Kyle Smith
January 29,2018
"""

import obspy
from obspy.io.sac import SACTrace
import obspy.signal.rotate as rotate
import os
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
from obspy.taup import TauPyModel
import matplotlib.pyplot as plt
import math

phases=["P"]
phases2=["p"]
sourcedepth = 100
taup_model = "ak135"
dists = list(np.arange(0,180,1))
model = TauPyModel(model=taup_model)
Phase1arrivals = []
Phase2arrivals = []
for dist_deg in dists[:]:
    temparr = model.get_travel_times(source_depth_in_km=sourcedepth,distance_in_degree=dist_deg,phase_list=[phases[0]])
    temparr2 = model.get_travel_times(source_depth_in_km=sourcedepth,distance_in_degree=dist_deg,phase_list=[phases2[0]])
    #print(len(temparr))
    if len(temparr)==0:
        Phase1arrivals.append(math.nan)
    else:
        Phase1arrivals.append(temparr[0].time)
    if len(temparr2)==0:
        Phase2arrivals.append(math.nan)
    else:
        Phase2arrivals.append(temparr2[0].time)
dkm = dists
#dkm = [a*111 for a in dists]
L2, = plt.plot(dkm,Phase2arrivals,'ro',label=phases2[0])
L1, = plt.plot(dkm,Phase1arrivals,'b*',label=phases[0])
plt.xlim(0, 180)
plt.legend(handles=[L1,L2])
plt.ylabel("Time after EQ origin time (s)")
plt.xlabel("Distance (deg)")
plt.show()
#print(Phase1arrivals)
