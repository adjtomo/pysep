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

def run_get_waveform(c, event, idb,
                     min_dist=20, max_dist=300, before=100, after=300,
                     network='*', station = '*', channel='BH*', 
                     resample_freq=20, 
                     ifrotateRTZ=True, ifrotateUVW=False, ifCapInp=True, ifRemoveResponse=True,
                     ifDetrend=True, ifDemean=True, ifTaper=False, ifEvInfo=True,
                     scale_factor=10.0**2, ipre_filt = 1,
                     pre_filt=(0.005, 0.006, 10.0, 15.0),
                     icreateNull = 1,
                     ifFilter=False, fmin=.02, fmax=1, filter_type='bandpass', 
                     zerophase=False, corners=4, 
                     iplot_response = False):
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
    
    evtime = event.origins[0].time
    
    if idb==1:
        print("Preparing request for IRIS ...")
        # BK network doesn't return data when using the IRIS client.
        # this option switches to NCEDC if BK is 
        if "BK" in network:
            client_name = "NCEDC"
            print("\nWARNING. Request for BK network. Switching to NCEDC client")
            c = Client("NCEDC")
        else:
            client_name = "IRIS" 

        print("Download stations...")
        stations = c.get_stations(network=network, station=station, 
                                  channel=channel,
                                  starttime=evtime - before, endtime=evtime + after,
                                  level="response")
        inventory = stations    # so that llnl and iris scripts can be combined
        
        print(stations)
        
        sta_limit_distance(event, stations, min_dist=min_dist, max_dist=max_dist)
        
        print("Downloading waveforms...")
        bulk_list = make_bulk_list_from_stalist(
            stations, evtime - before, evtime + after, channel=channel)
        stream_raw = c.get_waveforms_bulk(bulk_list)
            
    elif idb==3:
        client_name = "LLNL"
        print("Preparing request for LLNL ...")

        # Get event an inventory from the LLNL DB.
        event_number = int(event.event_descriptions[0].text)
        # event = llnl_db_client.get_obspy_event(event)
        inventory = c.get_inventory()
        
        print("--> Total stations in LLNL DB: %i" % (
                len(inventory.get_contents()["stations"])))
        sta_limit_distance(event, inventory, min_dist=min_dist, max_dist=max_dist)
        print("--> Stations after filtering for distance: %i" % (
                len(inventory.get_contents()["stations"])))

        stations = set([sta.code for net in inventory for sta in net])
        
        _st = c.get_waveforms_for_event(event_number)
        stream_raw = obspy.Stream()
        for tr in _st:
            if tr.stats.station in stations:
                stream_raw.append(tr)

    # set reftime
    stream = obspy.Stream()
    stream = set_reftime(stream_raw, evtime)

    if ifDemean:
        stream.detrend('demean')

    if ifDetrend:
        stream.detrend('linear')

    if ifFilter:
        prefilter(stream, fmin, fmax, zerophase, corners, filter_type)

    if ifRemoveResponse:
        resp_plot_remove(stream, ipre_filt, pre_filt, iplot_response, scale_factor, stations)
    else:
        # output RAW waveforms
        decon=False
        print("WARNING -- NOT correcting for instrument response")

    if scale_factor > 0:
        amp_rescale(stream, scale_factor)
        if idb ==3:
            amp_rescale_llnl(stream, scale_factor)

    print("--> Adding SAC metadata...")
    st2 = add_sac_metadata(stream,idb=idb, ev=event, stalist=inventory)

    # Set the sac header KEVNM with event name
    # This applies to the events from the LLNL database
    # NOTE this command is needed at the time of writing files, so it has to
    # be set early
    st2, evname_key = rename_if_LLNL_event(st2, evtime)

    #if resample_freq != 0:
    #    print("\n--> WARNING -- RESAMPLING")
    #    print("--> New sample rate = %5.1f\n" % resample_freq)
        # sporadical crashes (see iex = 18)
    #    resample(st2, freq=resample_freq)
        # appears to be stable
    #    resample_cut(st2, resample_freq,evtime, before, after)

    # Do some waveform QA
    # - (disabled) Throw out traces with missing data
    # - log waveform lengths and discrepancies
    # - Fill-in missing data -- Carl request
    do_waveform_QA(st2, client_name, event, evtime, before, after)

    # Get list of unique stations + locaiton (example: 'KDAK.00')
    stalist = []
    for tr in st2.traces:
        #stalist.append(tr.stats.station)
        stalist.append(tr.stats.network + '.' + tr.stats.station +'.'+ tr.stats.location + '.'+ tr.stats.channel[:-1])

    # Crazy way of getting a unique list of stations
    stalist = list(set(stalist))

    # match start and end points for all traces
    st2 = trim_maxstart_minend(stalist, st2, client_name, event, evtime, resample_freq, before, after)

    # save raw waveforms in SAC format
    path_to_waveforms = evname_key + "/RAW"
    write_stream_sac_raw(stream_raw, path_to_waveforms, evname_key, idb, event, stations=inventory)

    # Taper waveforms (optional; Generally used when data is noisy- example: HutchisonGhosh2016)
    # https://docs.obspy.org/master/packages/autogen/obspy.core.trace.Trace.taper.html
    # To get the same results as the default taper in SAC, use max_percentage=0.05 and leave type as hann.
    if ifTaper:
        st2.taper(max_percentage=0.5, type='hann',max_length=None, side='both')

    # save processed waveforms in SAC format
    path_to_waveforms = evname_key 
    write_stream_sac(st2, path_to_waveforms, evname_key)

    if ifrotateRTZ:
        rotate_and_write_stream(st2, evname_key, icreateNull, ifrotateUVW)

    if ifCapInp:
        write_cap_weights(st2, evname_key, client_name, event)

    if ifEvInfo:
        write_ev_info(event, evname_key)

