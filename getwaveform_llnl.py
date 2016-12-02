#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tools for interfacing the LLNL DB, ObsPy, and SAC input/output.
"""
from __future__ import print_function

import os

import obspy
from scipy import signal

from util_write_cap import *

def run_get_waveform(llnl_db_client, event, 
                     min_dist=20, max_dist=300, before=100, after=300, 
                     network='*', station = '*', channel='BH*',
                     resample_freq=20.0,
                     ifrotate=True, ifCapInp=True, ifRemoveResponse=True,
                     ifDetrend=True, ifDemean=True, ifEvInfo=True,
                     scale_factor=10.0**2,
                     ipre_filt = 1,
                     pre_filt=(0.005, 0.006, 10.0, 15.0), 
                     icreateNull=1,
                     ifFilter=False, fmin=0.02, fmax=1, filter_type='bandpass', 
                     zerophase=False, corners=4, 
                     iplot_response = False):
    """
    Get SAC waveforms for an event

    basic usage:
        run_get_waveform(ev)

    event -    The LLNL DB event number
    llnl_db_client - An active LLNLDB client instance.
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
    client_name = "LLNL"
    print("Preparing request for LLNL ...")

    # Database
    idb = 3

    # Get event an inventory from the LLNL DB.
    event_number = int(event.event_descriptions[0].text)
    #event = llnl_db_client.get_obspy_event(event)
    inventory = llnl_db_client.get_inventory()

    evtime = event.origins[0].time

    print("--> Total stations in LLNL DB: %i" % (
        len(inventory.get_contents()["stations"])))
    sta_limit_distance(event, inventory, min_dist=min_dist, max_dist=max_dist)
    print("--> Stations after filtering for distance: %i" % (
        len(inventory.get_contents()["stations"])))

    stations = set([sta.code for net in inventory for sta in net])

    _st = llnl_db_client.get_waveforms_for_event(event_number)
    _stream = obspy.Stream()
    for tr in _st:
        if tr.stats.station in stations:
            _stream.append(tr)

    # Create SAC objects and add headers    
    for tr in _stream:
        # This is the best way to make sac objects
        tr.write('tmppp.sac', format='SAC')
        tmptr = obspy.read('tmppp.sac').traces[0]
        tr.stats.sac = tmptr.stats.sac
        # Add units to the sac header (it should be in COUNTS before response is removed)
        tr.stats.sac['kuser0'] = 'llnl'

    # Save raw waveforms
    _stream = add_sac_metadata(_stream,idb=idb,ev=event, stalist=inventory) 
    evname_key=_stream[0].stats.sac['kevnm']
    write_stream_sac(_stream, evname_key=evname_key)
    os.rename(evname_key,'RAW')

    # set reftime
    stream = obspy.Stream()
    stream = set_reftime(_stream, evtime)

    #print("--> Extracted %i waveforms from the LLNL db database." % len(st))

    if ifDemean:
        stream.detrend('demean')

    if ifDetrend:
        stream.detrend('linear')

    if ifFilter:
        prefilter(stream, fmin, fmax, zerophase, corners, filter_type)

    if ifRemoveResponse:
        resp_plot_remove(stream, ipre_filt, iplot_response, scale_factor, stations)
    else:
        # output RAW waveforms
        decon=False
        print("WARNING -- NOT correcting for instrument response")

    if scale_factor > 0:
        amp_rescale(stream, scale_factor)

        # The response information in the LLNL database seems to have wrong
        # amplitude calibration. This correction rescales amplitudes by channel
        amp_rescale_llnl(stream, scale_factor)

    print("--> Adding SAC metadata...")
    st2 = add_sac_metadata(stream,idb=idb, ev=event, stalist=inventory)

    # Set the sac header KEVNM with event name
    # This applies to the events from the LLNL database
    # NOTE this command is needed at the time of writing files, so it has to
    # be set early
    st2, evname_key = rename_if_LLNL_event(st2, evtime)

    # 20160902 cralvizuri@alaska.edu -- 
    # see also correct_sac_tshift
    # This command overwrites SAC headers "b" and "e" with a timeshift. But
    # this is not neeeded since obspy handles this internally.
    # This command is now disabled.
    #time_shift_sac(st2, -1 * before)

    # Now do some QA: throw out traces with missing data
    # keep a log with status of trace extraction
    # this is mirrored in llnl_tool.py and iris_tools.py
    output_log = "data_processing_status" + "_" + client_name + ".log"
    fid = open(output_log, "w")
    fid.write("\n--------------------\n%s\n" % event.short_str())

    for tr in st2:
        fid.write("\n%s %s %s %s %s %s %6s %.2f sec" % (evtime, \
                tr.stats.network, tr.stats.station, tr.stats.channel, \
                tr.stats.starttime, tr.stats.endtime, tr.stats.npts, \
                float(tr.stats.npts / tr.stats.sampling_rate)))
        if (tr.stats.starttime > evtime - before) or \
                (tr.stats.endtime < evtime + after):
            print("===========")
            print("Event time:", evtime)
            print("Trace starttime:", tr.stats.starttime,
                  "Minimum required starttime:", evtime - before, "Pass:",
                  tr.stats.starttime > evtime - before)
            print("Trace endtime:", tr.stats.endtime,
                  "Maximum required endtime:", evtime + before,
                  "Pass:", tr.stats.endtime < evtime + after)
            print("Removing %s. Not in requested temporal range." % tr.id)
            print("===========")
            print("WARNING. Missing data for station %s" % tr.stats.station)
            print("WARNING. consider removing this station")
            fid.write(" -- data missing.")
            # 20160912 cralvizuri@alaska.edu --
            # the original code removes waveforms that do not have the same
            # length as the requested window. 
            # Rejection is now disabled per discussion today.
            # This is also mirrored in getwaveform_iris.py
            #st2.remove(tr)

    print("--> %i waveforms left." % len(st2))

    # fill gaps with 0
    # print(st2)    # code for debugging.
    #st2.merge(method=0,fill_value=0)
    st2.merge(fill_value='interpolate')
    # print(st2)    # code for debugging.

    fid.write("\n--------After filling the gaps------------")
    for tr in st2:
        fid.write("\n%s %s %s %s %s %s %6s %.2f sec" % (evtime, \
                tr.stats.network, tr.stats.station, tr.stats.channel, \
                tr.stats.starttime, tr.stats.endtime, tr.stats.npts, \
                float(tr.stats.npts / tr.stats.sampling_rate)))

    # Get list of unique stations + location (example: 'KDAK.00')
    stalist = []
    for tr in st2.traces:
        #stalist.append(tr.stats.station)
        stalist.append(tr.stats.network + '.' + tr.stats.station +'.'+ tr.stats.location + '.'+ tr.stats.channel[:-1])

    # Crazy way of getting a unique list of stations
    stalist = list(set(stalist))
    # print(stalist)    # for debugging.
    st2 = trim_maxstart_minend(stalist, st2)

    if not st2:
        raise Exception("No waveforms left")

    if resample_freq != 0:
        print("\n--> WARNING -- Resampling")
        print("--> New sample rate = %5.1f" % resample_freq)
        resample_cut(st2, resample_freq, evtime, before, after)

    write_stream_sac(st2, evname_key)
    # Move raw waveforms inside this directory
    os.rename('RAW',evname_key+'/RAW')

    if ifrotate:
        rotate_and_write_stream(st2, evname_key, icreateNull)

    if ifCapInp:
        write_cap_weights(st2, evname_key, client_name, event)

    if ifEvInfo:
        write_ev_info(event, evname_key)

    # 20160902 cralvizuri@alaska.edu -- 
    # see also time_shift_sac
    # This command overwrites SAC headers "b" and "e" with a timeshift. But
    # this is not neeeded since obspy handles this internally.
    # This command is now disabled.
    #Fix b and e sac headers
    #correct_sac_tshift(evtime.strftime('%Y%m%d%H%M%S%f')[:-3]+'/', \
    #        before, after)

