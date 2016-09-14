#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tools for interfacing IRIS data, ObsPy, and SAC input/output.
"""
from __future__ import print_function

import os

import obspy
from obspy.clients.fdsn import Client
from scipy import signal

from util_write_cap import *

def run_get_waveform(c, event,
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
    print("Preparing request for IRIS ...")

    # BK network doesn't return data when using the IRIS client.
    # this option switches to NCEDC if BK is 
    if "BK" in network:
        client_name = "NCEDC"
        print("\nWARNING. Request for BK network. Switching to NCEDC client")
        c = Client("NCEDC")
    else:
        client_name = "IRIS"

    evtime = event.preferred_origin().time

    print("Download stations...")
    stations = c.get_stations(network=network, starttime=evtime - before,
                              endtime=evtime + after, channel=channel,
                              level="response")

    print(stations)

    sta_limit_distance(event, stations, min_dist=min_dist, max_dist=max_dist)

    print("Downloading waveforms...")
    bulk_list = make_bulk_list_from_stalist(
        stations, evtime - before, evtime + after, channel=channel)
    _stream = c.get_waveforms_bulk(bulk_list)
    print(_stream)

    # set reftime
    stream = obspy.Stream()
    stream = set_reftime(_stream, evtime)

    if ifDemean:
        stream.detrend('demean')

    if ifDetrend:
        stream.detrend('linear')

    if ifRemoveResponse:
        stream.remove_response(inventory=stations, pre_filt=pre_filt,
                               output="VEL")

    if scale_factor > 0:
        print("\n--> WARNING -- rescaling amplitudes by %f" % scale_factor)
        for tr in stream.traces:
            tr.data = tr.data * scale_factor

    stream.detrend('demean')

    print("--> Adding SAC metadata...")
    st2 = add_sac_metadata(stream, ev=event, stalist=stations)

    # 20160902 cralvizuri@alaska.edu -- 
    # see also correct_sac_tshift
    # This command overwrites SAC headers "b" and "e" with a timeshift. But
    # this is not neeeded since obspy handles this internally.
    # This command is now disabled.
    #time_shift_sac(st2, -1 * before)

    if resample_freq != 0:
        print("\n--> WARNING -- RESAMPLING")
        print("--> New sample rate = %5.1f\n" % resample_freq)
        resample(st2, freq=resample_freq)

    # Now do some QA: throw out traces with missing data
    # keep a log with status of trace extraction
    # this is mirrored in llnl_tool.py and iris_tools.py
    output_log = "data_processing_status" + "_" + client_name + ".log"
    fid = open(output_log, "w")
    fid.write("\n--------------------\n%s\n" % event.short_str())
    fid.write("evtime net sta cha starttime endtime npts length (sec)\n")
    for tr in st2:
        fid.write("\n%s %s %s %s %s %s %6s %.2f sec" % (evtime, \
                tr.stats.network, tr.stats.station, tr.stats.channel, \
                tr.stats.starttime, tr.stats.endtime, tr.stats.npts, \
                float(tr.stats.npts / tr.stats.sampling_rate)))
        if tr.stats.npts < tr.stats.sampling_rate * (before + after):
            print("WARNING. Missing data for station %s" % tr.stats.station)
            print("WARNING. consider removing this station")
            fid.write(" -- data missing.")
            # 20160912 cralvizuri@alaska.edu --
            # the original code removes waveforms that do not have the same
            # length as the requested window. 
            # Rejection is now disabled per discussion today.
            # This is also mirrored in getwaveform_llnl.py
            st2.remove(tr)

    write_stream_sac(st2, evtime)

    if ifrotate:
        rotate_and_write_stream(st2, evtime)

    if ifCapInp:
        write_cap_weights(st2, evtime, client_name)

    if ifEvInfo:
        write_ev_info(event, evtime)

    # 20160902 cralvizuri@alaska.edu -- 
    # see also time_shift_sac
    # This command overwrites SAC headers "b" and "e" with a timeshift. But
    # this is not neeeded since obspy handles this internally.
    # This command is now disabled.
    #Fix b and e sac headers
    #correct_sac_tshift(evtime.strftime('%Y%m%d%H%M%S%f')[:-3]+'/', \
    #        before, after)

