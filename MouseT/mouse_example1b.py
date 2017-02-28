#! /usr/bin/env python

"""
Example of possible usage of MouseTrap module: detection of "mouse" disturbances in a long record.

It internally cuts record into overlapping time windows and try to detect one disturbance in each time window.

Limitations:
 - If two disturbances are present in one time window, only one of them is found.
 - If a disturbance is closer than `t_mean` to the starttime of the record, it might not be found and the fit value is misleading
 

"""

import numpy as np
from obspy.core import read
from obspy.clients.arclink import Client #ArcLink
from obspy.io.sac.sacpz import attach_paz # load poles and zeros from SAC PZ file
from MouseTrap import *

# constants
t_fit = 10*60 # length of time window, where the mouse onset is searched (second)
t_mean = 2*60 # time for calculating mean value in the bebinning of each time window (second)
mouse_len = 5*60 # syntetic mouse length in second
mouse_onset = 2*60 # the onset of the syntetic mouse is `mouse_onset` seconds after synthetic mouse starttime
bitrate = 10 # Hz, downsample the record to this sampling
t_before = max(mouse_onset, t_mean)

# read waveform file
st = read('VANNI_1hr.mseed')
t_start = max(st[0].stats.starttime, st[1].stats.starttime, st[2].stats.starttime)
t_len = min(st[0].stats.endtime, st[1].stats.endtime, st[2].stats.endtime) - t_start

# read poles and zeros
paz_file = 'SAC_PZs_CH_VANNI_HHZ__2009.303.00.00.00.0000_99999.9999.24.60.60.99999'
attach_paz(st[0], paz_file, tovel=True)
paz = st[0].stats.paz

factor = int(st[0].stats.sampling_rate / bitrate) # will be downsampled by this factor
dt = st[0].stats.delta * factor # then the time delta will be this
m1 = mouse(fit_time_before = mouse_onset, fit_time_after = mouse_len-mouse_onset) # new object of MouseTrap
m1.create(paz, int(mouse_len/dt), dt, mouse_onset) # create synthetic mouse bite of type 1 = step in acceleration

# cut records into time windows and try to find the mouse in each of them
#for t in range (0, int(t_len-t_before-mouse_len+mouse_onset), int(t_fit-t_before-mouse_len+mouse_onset)):
t = 0
while t < t_len-t_before-mouse_len+mouse_onset:
	if t > t_len-t_fit: # last time window can overlap with the previous one, otherwise it could be shorter and make troubles
		t = t_len-t_fit
	st2 = st.slice(t_start+t, t_start+t_fit+t)
	demean(st2, t_mean) # remove mean calculated from the beginning of the record
	ToDisplacement(st2, 10) # integrate, subsample to 10 Hz

	# fit waveform by synthetic m1
	m1.fit_3D(st2, t_min=mouse_onset, t_max=t_fit-mouse_len+mouse_onset)
	m1.plot(st2, outfile='VANNI_{0:0>4.0f}'.format(t), xmin=mouse_onset-20, xmax=t_fit, ylabel='raw displacement [count]')
	onset, amp, phi, theta, fit = m1.params(degrees=True)
	
	T = t # for printing output
	if amp > 1e-7 and fit > 0.5: # mouse existence criterion, THIS IS ONLY EXAMPLE
		mouse_present = True
		t += onset+mouse_len-mouse_onset
	else:
		mouse_present = False
		t += t_fit-t_before-mouse_len+mouse_onset
	print('=== block t = %6.1f ===' % T)
	print('%s' %{True:'MOUSE', False:'no mouse'}[mouse_present])
	print('time of onset:   %6.1f s'% (onset+T))
	print('amplitude:   %10.2e m s^-2' % amp)
	print('phi:             %6.1f deg' % phi)
	print('theta:           %6.1f deg' % theta)
	print('fit:            %7.2f'% fit)
