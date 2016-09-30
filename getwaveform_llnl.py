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
                     network='*', station = '*', channel='BH*', resample_freq=20.0,
                     ifrotate=True, ifCapInp=True, ifRemoveResponse=True,
                     ifDetrend=True, ifDemean=True, ifEvInfo=True,
                     scale_factor=10.0**2,
                     pre_filt=(0.005, 0.006, 10.0, 15.0)):
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
    _st2 = obspy.Stream()
    for tr in _st:
        if tr.stats.station in stations:
            _st2.append(tr)

    # set reftime
    st = obspy.Stream()
    st = set_reftime(_st2, evtime)

    print("--> Extracted %i waveforms from the LLNL db database." % len(st))

    if ifDemean:
        st.detrend('demean')

    if ifDetrend:
        st.detrend('linear')

    if ifRemoveResponse:
        decon=True
        passed_st = obspy.Stream()
        failed_st = []
        for tr in st:
            if tr.stats.channel[-1] not in ["Z", "N", "E"]:
                print("%s is not vertical, north, or east. Skipped." % tr.id)
                continue
            try:
                llnl_db_client.remove_response(
                        tr=tr, pre_filt=pre_filt, output="VEL")
            except Exception as e:
                print("Failed to correct %s due to: %s" % (tr.id, str(e)))
                failed_st.append(tr)
                continue
            else:
                passed_st.append(tr)
        st = passed_st
        if not st:
            raise Exception("No waveform managed to get instrument corrected")

        print("--> %i waveforms failed to instrument correct" % len(failed_st))
        print("--> %i waveforms managed to instrument correct" % len(passed_st))

    else:
        # output RAW waveforms
        decon=False
        print("WARNING -- NOT correcting for instrument response")

    if scale_factor > 0:
        print("\n--> WARNING -- rescaling amplitudes by %f" % scale_factor)
        for tr in st.traces:
            tr.data = tr.data * scale_factor

    st.detrend('demean')

    print("--> Adding SAC metadata...")
    st2 = add_sac_metadata(st, ev=event, stalist=inventory)

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

    if not st2:
        raise Exception("No waveforms left")

    if resample_freq != 0:
        print("\n--> WARNING -- Resampling")
        print("--> New sample rate = %5.1f" % resample_freq)
        resample_cut(st2, resample_freq, evtime, before, after)

    write_stream_sac(st2, evtime)

    if ifrotate:
        rotate_and_write_stream(st2, evtime)

    if ifCapInp:
        write_cap_weights(st2, evtime, client_name, event)

    if ifEvInfo:
        write_ev_info(event, evtime)

    if decon is False:
        print("WARNING waveforms NOT corrected for instrument response")

    # 20160902 cralvizuri@alaska.edu -- 
    # see also time_shift_sac
    # This command overwrites SAC headers "b" and "e" with a timeshift. But
    # this is not neeeded since obspy handles this internally.
    # This command is now disabled.
    #Fix b and e sac headers
    #correct_sac_tshift(evtime.strftime('%Y%m%d%H%M%S%f')[:-3]+'/', \
    #        before, after)

