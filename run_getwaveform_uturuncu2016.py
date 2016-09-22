#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to fetch waveform data from IRIS for all 63 events in

    C. Alvizuri and C. Tape. Full moment tensors for small events (Mw < 3) at
    Uturuncu volcano, Bolivia. Geophys. J. Int., 206(3):1761{1783, 2016. doi:
    10.1093/gji/ggw247.

    Usage:
        python run_getwaveform_alvizuri2016.py

20160920 cralvizuri <cralvizuri@alaska.edu>
"""

import obspy
import getwaveform_iris
from obspy.clients.fdsn import Client
from obspy.core.event import Event, Origin, Magnitude
import util_helpers as uh

# flag to request data or not
send_request = True

# read catalog file to build event object
fmt_cat_uturuncu = "/home/carltape/REPOSITORIES/manuscripts/alvizuri/papers/2014fmt/data/uturuncu_mech.txt"
skip_headers = 22   # starts at index = 0

# settings for preparing data for CAP
rotate = True
output_cap_weight_file = True
remove_response = True
detrend = True
demean = True
output_event_info = True
# pre-filter for deconvolution
fminc = 1/20   # see cutoff = [0.05 20] in rs_bolivia.m
fmaxc = 20
pre_filt=(0.5*fminc, fminc, fmaxc, 2.0*fmaxc)    
resample_freq = 0       # For Uturuncu we use the original sample rate (100)
scale_factor = 10**2    # for CAP use 10**2 (to convert m/s to cm/s)

# parameters for extracting waveforms
min_dist = 0 
max_dist = 500
tbefore_sec = 50    # see (tstart, tend) in rs_bolivia.m
tafter_sec = 150
network = 'XP'
station = '*'
channel = 'HH*'

# read event info from the catalogs, create event objects, request data to IRIS
f = open(fmt_cat_uturuncu, "r")
lines = f.readlines()[skip_headers:]
f.close()

# NOTE this loop sends a request for each line in the catalog (63 events)
for line in lines:
    line_elements = line.split()
    eid = line_elements[-1]
    otime = uh.eid2otime(eid)
    elon = line_elements[6]
    elat = line_elements[7]
    edep = float(line_elements[8]) * 1000.0 # meters
    emag = line_elements[16]

    # create event object
    orig = Origin()
    orig.longitude = elon
    orig.latitude = elat
    orig.depth = edep
    orig.time = otime
    mag = Magnitude()
    mag.mag = emag
    mag.magnitude_type = "Mw"
    ev = Event()
    ev.origins.append(orig)
    ev.magnitudes.append(mag)

    if send_request:
        # get waveforms
        client = Client("IRIS")
        getwaveform_iris.run_get_waveform(c = client, event = ev, 
                                      min_dist = min_dist, max_dist = max_dist, 
                                      before = tbefore_sec, after = tafter_sec, 
                                      network = network, station = station, channel = channel, 
                                      resample_freq = resample_freq, ifrotate = rotate,
                                      ifCapInp = output_cap_weight_file, 
                                      ifRemoveResponse = remove_response,
                                      ifDetrend = detrend, ifDemean = demean, 
                                      ifEvInfo = output_event_info, 
                                      scale_factor = scale_factor, pre_filt = pre_filt)
    else:
        print(orig)
