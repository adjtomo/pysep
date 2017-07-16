#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Template script for fetching waveforms from IRIS using the utilities in
getwaveform*py 

This example is setup for the 2008 Mount Carmel, Illinois, event (Mw5.2)

    Usage:
        python run_getwaveform_illinois.py

Contributors: C. Alvizuri, V. Silwal, C. Tape
"""

import obspy
import getwaveform
from obspy.clients.fdsn import Client
from obspy.core.event import Event, Origin, Magnitude
from obspy.core import UTCDateTime
import util_helpers as uh
import sys

# parameters for extracting waveforms
send_request = True     # comment
rotateRTZ = True        # comment
rotateUVW = False       # works only if 'rotateRTZ = True'
output_cap_weight_file = True   # comments
remove_response = True          # here
iplot_response = False          # would
detrend = True                  # be
demean = True                   # nice
output_event_info = True        # here too
fminc = 1/20                    # see cutoff = [0.05 20] in rs_bolivia.m
fmaxc = 20                      # and here
pre_filt = (0.5*fminc, fminc, fmaxc, 2.0*fmaxc)    
resample_TF = False     # resample_freq is not used if this is False
resample_freq = 0       # For Uturuncu we use the original sample rate (100)
scale_factor = 10**2    # for CAP use 10**2 (to convert m/s to cm/s)
idb = 1                 # default: =1-IRIS; =2-AEC; =3-LLNL

min_dist = 0        # more comments here as well
max_dist = 500      # and so on..
tbefore_sec = 50    # (hope this catches on)
tafter_sec = 150    # 
network = 'IU,NM'
station = '*'
channel = 'BH?,LH?'

# Event info. This is built by the user. (What are other alternatives?)
origins_quakes = {
        "Illinois": ["2008-04-18T09:36:59.110", "-87.886", "38.452", "14300.0", "5.2 "], 
        }

# NOTE this loop sends a request for each line in the catalog (63 events)
for line in origins_quakes.items():

    # Build the Magnitude(), Event(), Origin() objects to submit the request
    # to the data center (IRIS, etc)

    # Get event info from the user (this could change)
    otime = UTCDateTime(line[1][0])
    # Another option is to use a simpler evid and convert to UTC format, eg
    # eid = 200804180936 
    # otime = uh.eid2otime(eid)
    elon = line[1][1]
    elat = line[1][2]
    edep = line[1][3]
    emag = line[1][4]

    # Create the event objects
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

    # Request the data
    if send_request:
        # get waveforms
        client = Client("IRIS")
        getwaveform.run_get_waveform(c = client, # consider add some comments
                event = ev,                      #
                ref_time_place=ev,
                min_dist = min_dist,
                max_dist = max_dist, 
                before = tbefore_sec, after = tafter_sec, 
                network = network, station = station, channel = channel, 
                resample_freq = resample_freq, 
                pre_filt = pre_filt,
                ifrotateRTZ = rotateRTZ,
                ifrotateUVW = rotateUVW,
                ifCapInp = output_cap_weight_file, 
                ifRemoveResponse = remove_response,
                ifDetrend = detrend, ifDemean = demean, 
                ifEvInfo = output_event_info, 
                scale_factor = scale_factor,
                iplot_response = iplot_response,
                idb = idb
                )
    else:
        print(orig)

