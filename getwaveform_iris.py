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
                     network='*', station = '*', channel='BH*', 
                     resample_freq=20, 
                     ifrotate=True, ifCapInp=True, ifRemoveResponse=True,
                     ifDetrend=True, ifDemean=True, ifEvInfo=True,
                     scale_factor=10.0**2,
                     pre_filt=(0.005, 0.006, 10.0, 15.0), icreateNull=1,
                     ifFilter=False, fmin=.02, fmax=1, filt_type='bandpass', 
                     zerophase=False, corners=4):
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

    #evtime = event.preferred_origin().time
    evtime = event.origins[0].time

    print("Download stations...")
    stations = c.get_stations(network=network, station=station, 
            channel=channel,
            starttime=evtime - before, endtime=evtime + after,
            level="response")

    print(stations)

    sta_limit_distance(event, stations, min_dist=min_dist, max_dist=max_dist)

    print("Downloading waveforms...")
    bulk_list = make_bulk_list_from_stalist(
        stations, evtime - before, evtime + after, channel=channel)
    _stream = c.get_waveforms_bulk(bulk_list)
    #print(_stream)  # code for debugging. partial print stations

    # Add headers to the raw waveforms
    # Save raw waveforms
    _stream = add_sac_metadata(_stream,ev=event, stalist=stations) 
    evname_key=_stream[0].stats.sac['kevnm']
    write_stream_sac(_stream, evname_key=evname_key)
    os.rename(evname_key,'RAW')

    # set reftime
    stream = obspy.Stream()
    stream = set_reftime(_stream, evtime)

    if ifDemean:
        stream.detrend('demean')

    if ifDetrend:
        stream.detrend('linear')

    if ifFilter:
        for tr in stream:
            print('Filtering ', tr.stats.network +'.'+ tr.stats.station +'.'+ tr.stats.location +'.'+ tr.stats.channel)
            tr.filter('bandpass',freqmin=fmin,freqmax=fmax,zerophase=False,corners=corners)

    if ifRemoveResponse:
        for tr in stream:
            print('Removing instrument response from ' + tr.stats.network +'.'+ tr.stats.station +'.'+ tr.stats.location +'.'+ tr.stats.channel)
            tr.remove_response(inventory=stations, pre_filt=pre_filt,
                               output="VEL")

    if scale_factor > 0:
        print("\n--> WARNING -- rescaling amplitudes by %f" % scale_factor)
        for tr in stream.traces:
            tr.data = tr.data * scale_factor

    stream.detrend('demean')

    print("--> Adding SAC metadata...")
    st2 = add_sac_metadata(stream, ev=event, stalist=stations)

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
    fid.write("evtime net sta loc cha starttime endtime npts length (sec)\n")
    for tr in st2:
        fid.write("\n%s %s %s %s %s %s %s %6s %.2f sec" % (evtime, \
                tr.stats.network, tr.stats.station, tr.stats.location, tr.stats.channel, \
                tr.stats.starttime, tr.stats.endtime, tr.stats.npts, \
                float(tr.stats.npts / tr.stats.sampling_rate)))
        if tr.stats.npts < tr.stats.sampling_rate * (before + after):
            print("WARNING. data available < (before + after) for station " + \
                    tr.stats.network +'.'+ tr.stats.station +'.'+ tr.stats.location +'.'+ tr.stats.channel + " -- consider removing this station")
            fid.write(" -- data missing.")
            # 20160912 cralvizuri@alaska.edu --
            # the original code removes waveforms that do not have the same
            # length as the requested window. 
            # Rejection is now disabled per discussion today.
            # This is also mirrored in getwaveform_llnl.py
            # st2.remove(tr)

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

    # Get list of unique stations + locaiton (example: 'KDAK.00')
    stalist = []
    for tr in st2.traces:
        #stalist.append(tr.stats.station)
        stalist.append(tr.stats.network + '.' + tr.stats.station +'.'+ tr.stats.location + '.'+ tr.stats.channel[:-1])

    # Crazy way of getting a unique list of stations
    stalist = list(set(stalist))
    # print(stalist)    # for debugging.
    st2 = trim_maxstart_minend(stalist, st2)

    #st3 = obspy.Stream()
    ## Trim the edges in case 3 channels have different lengths
    #for stn in stalist:
    #    # split STNM.LOC
    #    tmp = stn.split('.')
    #    netw = tmp[0]
    #    station = tmp[1]
    #    location = tmp[2]
    #    chan = tmp[3] + '*'
    #    # Get 3 traces (subset based on matching station name and location code)
    #    substr = stream.select(network=netw,station=station,location=location,channel=chan)
    #    # Find max startime and min end time for stations with number of channels = 1 or 2 or 3
    #    if len(substr) == 1:
    #        max_starttime = substr[0].stats.starttime
    #        min_endtime = substr[0].stats.endtime
    #    if len(substr) == 2:
    #        max_starttime = max(substr[0].stats.starttime,substr[1].stats.starttime)
    #        min_endtime = min(substr[0].stats.endtime,substr[1].stats.endtime)
    #    if len(substr) == 3:
    #        max_starttime = max(substr[0].stats.starttime,substr[1].stats.starttime,substr[2].stats.starttime)
    #        min_endtime = min(substr[0].stats.endtime,substr[1].stats.endtime,substr[2].stats.endtime)
    #    print(substr[0].stats.station, max_starttime, min_endtime)
    #    try:
    #        substr.trim(starttime=max_starttime, endtime=min_endtime, pad=False, nearest_sample=True, fill_value=0)
    #    except:
    #        print('WARNING: stattime larger than endtime for channels of', netw, '.', station, '.', location)
    #        continue
    #    for tr in substr.traces:
    #        st3 = st3.append(tr)
    #    
    #st2=st3

    fid.write("\n--------After trimming the edges (in case the 3 channels have different lengths)------------")
    for tr in st2:
        fid.write("\n%s %s %s %s %s %s %6s %.2f sec" % (evtime, \
                tr.stats.network, tr.stats.station, tr.stats.channel, \
                tr.stats.starttime, tr.stats.endtime, tr.stats.npts, \
                float(tr.stats.npts / tr.stats.sampling_rate))) 

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

