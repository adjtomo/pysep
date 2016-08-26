#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tools for interfacing IRIS data, ObsPy, and SAC input/output.
"""
from __future__ import print_function

import os

import numpy as np
import obspy
from obspy.clients.fdsn import Client
from scipy import signal

from util_write_cap import *

def resample(st, freq):
    """
    Custom resampling with a very sharp zerophase filter.
    """
    new_nyquist = 0.5 * freq
    for tr in st:
        current_nyquist = 0.5 * tr.stats.sampling_rate
        # Filter if necessary.
        if new_nyquist < current_nyquist:
            zerophase_chebychev_lowpass_filter(trace=tr, freqmax=new_nyquist)
        tr.detrend("linear")
        tr.taper(max_percentage=0.02, type="hann")
        tr.data = np.require(tr.data, requirements=["C"])
        # "Perfect" sinc resampling.
        tr.interpolate(sampling_rate=freq, method="lanczos", a=8,
                       window="blackman")

def run_get_waveform(event,
                     min_dist=20, max_dist=300, before=100, after=300,
                     network='*', channel='BH*', resample_freq=20.0, 
                     ifrotate=True, ifCapInp=True, ifRemoveResponse=True,
                     ifDetrend=True, ifDemean=True, ifEvInfo=True,
                     scale_factor=10.0**2,
                     pre_filt=(0.005, 0.006, 10.0, 15.0)):
    """
    Get SAC waveforms for an event

    basic usage:
        run_get_waveform(event)

    event -       obspy Event object

    min_dist - minimum station distance (default = 20)
    max_dist - maximum station distance (default =300)
    before -   time window length before the event time (default= 100)
    after  -   time window length after the event time (default = 300)
    network -  network codes of potential stations (default=*)
    channel -  component(s) to get, accepts comma separated (default='BH*')
    resample_freq- sampling frequency to resample files (default 20.0, 0 for
                                                     no resampling)
    ifrotate - Boolean, if true will output sac files rotated to baz
               unrotated sac files will also be written
    ifCapInp - Boolean, make weight files for CAP
    ifEvInfo - Boolean, output 'ev_info.dat' containg event info (True)
    ifRemoveResponse - Boolean, will remove response (True)
    ifDetrend - Boolean, will remove linear trend from data (True)
    ifDemean  - Boolean, will insult the data (True)
    scale_factor - scale all data by one value (10.0**2)
                    This usually puts the data in the units required by CAP
                    From m/s to cm/s
    pre_filt  - list, corner frequencies of filter to apply before deconv
                a good idea when deconvolving (ifRemoveResponse=True)
    """
    # BK network doesn't return data when using the IRIS client.
    # this option switches to NCEDC if BK is 
    if "BK" in network:
        print("\nWARNING. BK network. Using NCEDC client")
        c = Client("NCEDC")
    else:
        print("\nUsing IRIS client")
        c = Client("IRIS")

    evtime = event.preferred_origin().time

    print("Download stations...")
    stations = c.get_stations(network=network, starttime=evtime - before,
                              endtime=evtime + after, channel=channel,
                              level="response")

    print(stations)

    sta_limit_distance(event, stations, min_dist=min_dist, max_dist=max_dist)

    print("Downloading waveforms...")
    bulk_list = make_bulk_list_from_stalist(
        stations, evtime - before, evtime + after, cmp=channel)
    stream = c.get_waveforms_bulk(bulk_list)
    print(stream)

    if ifDemean:
        stream.detrend('demean')

    if ifDetrend:
        stream.detrend('linear')

    if ifRemoveResponse:
        stream.remove_response(inventory=stations, pre_filt=pre_filt,
                               output="VEL")

    if scale_factor > 0:
        for tr in stream.traces:
            tr.data = tr.data * scale_factor

    stream.detrend('demean')

    st2 = add_sac_metadata(stream, ev=event, stalist=stations)

    time_shift_sac(st2, -1 * before)

    if resample_freq != 0:
        print("\n--> !! WARNING !! Resampling...\n")
        resample(st2, freq=resample_freq)

    # Now do some QA: throw out traces with missing data
    # keep a log with status of trace extraction
    # this is mirrored in llnl_tool.py and iris_tools.py
    outlog = "get_data_status_IRIS.log"
    fid = open(outlog, "a")
    fid.write("\n--------------------\n%s\n" % event.short_str())

    for tr in st2:
        fid.write("\n%s %s %s %s %s %s %6s %.2f sec" % (evtime, \
                tr.stats.network, tr.stats.station, tr.stats.channel, \
                tr.stats.starttime, tr.stats.endtime, tr.stats.npts, \
                float(tr.stats.npts / tr.stats.sampling_rate)))
        if tr.stats.npts < tr.stats.sampling_rate * (before + after):
            print("WARNING. missing data for station %s" % tr.stats.station)
            print("WARNING. removing this station")
            fid.write(" -- data missing. removing.")
            st2.remove(tr)

    write_stream_sac(st2, evtime)

    if ifrotate:
        rotate_and_write_stream(st2, evtime)

    if ifCapInp:
        write_cap_weights(st2, evtime)

    if ifEvInfo:
        write_ev_info(event, evtime)

    #Fix b and e sac headers
    correct_sac_tshift(evtime.strftime('%Y%m%d%H%M%S%f')[:-3]+'/', \
            before, after)

