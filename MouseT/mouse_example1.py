#! /usr/bin/env python

import numpy as np
from obspy.core import read
from obspy.arclink import Client #ArcLink
from obspy.sac.sacio import attach_paz # load poles and zeros from SAC PZ file
from MouseTrap import *

# constants
t_start_origin = 120 # records start-time is `t_start_origin` seconds before the event origin time

# read waveform file
st = read('VANNI.mseed')

# read poles and zeros
paz_file = 'SAC_PZs_CH_VANNI_HHZ__2009.303.00.00.00.0000_99999.9999.24.60.60.99999'
attach_paz(st[0], paz_file, tovel=True)
paz = st[0].stats.paz

# demean, integrate, check signal-to-noise ratio
error = PrepareRecord(st, t_start_origin)
if error:
	print ('    %s' % error)
	exit()

# create synthetic mouse
length = max(len(st[0]), len(st[1]), len(st[2]))
dt = st[0].stats.delta
mouse = mouse()
mouse.create(paz, length*2, dt, length*dt, 1)

# fit waveform by synthetic mouse
mouse.fit_3D(st, t_min=105, t_max=210)
mouse.plot(st, outfile='VANNI.png', xmax=300, ylabel='raw displacement')

if mouse.exist(t_start_origin):
	print('=== MOUSE DETECTED ===')
	print('time of onset:   %6.1f s' % mouse.onset)
	print('amplitude:   %10.2e m s^-2' % mouse.amplitude)
	print('phi:             %6.1f deg' % (mouse.phi*180./np.pi))
	print('theta:           %6.1f deg' % (mouse.theta*180./np.pi))
	print('fit:            %7.2f' % mouse.fit)
else:
	print('=== no mouse detected ===')
