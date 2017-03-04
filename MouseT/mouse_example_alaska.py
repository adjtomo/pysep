#! /usr/bin/env python

import numpy as np
from obspy.core import read
from obspy.clients.arclink import Client #ArcLink
from obspy.io.sac.sacpz import attach_paz # load poles and zeros from SAC PZ file
from MouseTrap import *
import glob

# constants
t_start_origin = 100 # records start-time is `t_start_origin` seconds before the event origin time

# Event info
ddir = '/home/vipul/REPOSITORIES/GEOTOOLS/python_util/util_data_syn/20161208101813000/'
net = '*'
station = 'BPAW'

# read waveform file
waveform_dir = ddir + 'RAW/'
st = read(waveform_dir + '*' + net + '*' + station + '*')

# read poles and zeros
resp_dir = ddir + 'resp/'
paz_file = resp_dir + '*' + station + '*' + 'resp'
paz_list = glob.glob(paz_file)
attach_paz(st[0], paz_list[0], tovel=True)
paz = st[0].stats.paz

# demean, integrate, check signal-to-noise ratio
error = PrepareRecord(st, t_start_origin)
if error:
	print ('    %s' % error)
	exit()

# create synthetic mouse
length = max(len(st[0]), len(st[1]), len(st[2]))
dt = st[0].stats.delta
mouse = mouse(fit_time_before = 100,fit_time_after = 400)
#mouse = mouse()
mouse_length = length*2
mouse_onset = 100
print(mouse_length,mouse_onset)
mouse.create(paz, mouse_length, dt, mouse_onset, 2)

# fit waveform by synthetic mouse
mouse.fit_3D(st, t_min=80, t_max=120)
mouse.plot(st, distance=5e5, outfile=station+'.eps', xmin = 0, xmax=600, ylabel='raw displacement')
mouse.remove(st)

if mouse.exist(t_start_origin):
	print('=== MOUSE DETECTED ===')
	print('time of onset:   %6.1f s' % mouse.onset_found)
	print('amplitude:   %10.2e m s^-2' % mouse.amplitude)
	print('phi:             %6.1f deg' % (mouse.phi*180./np.pi))
	print('theta:           %6.1f deg' % (mouse.theta*180./np.pi))
	print('fit:            %7.2f' % mouse.fit)
else:
	print('=== no mouse detected ===')
