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


def zerophase_chebychev_lowpass_filter(trace, freqmax):
    """
    Custom Chebychev type two zerophase lowpass filter useful for
    decimation filtering.
    This filter is stable up to a reduction in frequency with a factor of
    10. If more reduction is desired, simply decimate in steps.
    Partly based on a filter in ObsPy.
    :param trace: The trace to be filtered.
    :param freqmax: The desired lowpass frequency.
    Will be replaced once ObsPy has a proper decimation filter.
    """
    # rp - maximum ripple of passband, rs - attenuation of stopband
    rp, rs, order = 1, 96, 1e99
    ws = freqmax / (trace.stats.sampling_rate * 0.5)  # stop band frequency
    wp = ws  # pass band frequency

    while True:
        if order <= 12:
            break
        wp *= 0.99
        order, wn = signal.cheb2ord(wp, ws, rp, rs, analog=0)

    b, a = signal.cheby2(order, rs, wn, btype="low", analog=0, output="ba")

    # Apply twice to get rid of the phase distortion.
    trace.data = signal.filtfilt(b, a, trace.data)


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


def run_get_waveform(ev, min_dist=20, max_dist=300, before=100, after=300,
                     network='*', channel='BH*', samp_freq=20.0, rotate=True,
                     output_cap_weight_file=True, remove_response=True,
                     detrend=True, demean=True, output_event_info=True,
                     scale_factor=10.0**2,
                     pre_filt=(0.005, 0.006, 10.0, 15.0)):
    """
    Get SAC waveforms for an event

    basic usage:
        run_get_waveform(ev)

    ev -       obspy Event object
    min_dist - minimum station distance (default = 50)
    max_dist - maximum station distance (default =200)
    before -   time window length before the event time (default= 100)
    after  -   time window length after the event time (default = 300)
    network -  network codes of potential stations (default=*)
    channel -  component(s) to get, accepts comma separated (default='BH*')
    samp_freq- sampling frequency to resample files (default 20.0, 0 for
                                                     no resampling)
    rotate - Boolean, if true will output sac files rotated to baz
               unrotated sac files will also be written
    output_cap_weight_file - Boolean, make weight files for CAP
    output_event_info- Boolean, output 'ev_info.dat' containg event info (True)
    remove_response - Boolean, will remove response (True)
    detrend - Boolean, will remove linear trend from data (True)
    demean  - Boolean, will insult the data (True)
    scale_factor - scale all data by one value (10.0**2)
                    This usually puts the data in the units required by CAP
                    From m/s to cm/s
    pre_filt  - list, corner frequencies of filter to apply before deconv
                a good idea when deconvolving (remove_response=True)
    """
    # BK network doesn't return data when using the IRIS client.
    # this option switches to NCEDC if BK is 
    if "BK" in network:
        print("\nWARNING. BK network. Using NCEDC client")
        c = Client("NCEDC")
    else:
        print("\nUsing IRIS client")
        c = Client("IRIS")

    evtime = ev.preferred_origin().time

    print("Download stations...")
    stations = c.get_stations(network=network, starttime=evtime - before,
                              endtime=evtime + after, channel=channel,
                              level="response")

    print(stations)

    sta_limit_distance(ev, stations, min_dist=min_dist, max_dist=max_dist)

    print("Downloading waveforms...")
    bulk_list = make_bulk_list_from_stalist(
        stations, evtime - before, evtime + after, cmp=channel)
    stream = c.get_waveforms_bulk(bulk_list)
    print(stream)

    if demean:
        stream.detrend('demean')

    if detrend:
        stream.detrend('linear')

    if remove_response:
        stream.remove_response(inventory=stations, pre_filt=pre_filt,
                               output="VEL")

    if scale_factor > 0:
        for tr in stream.traces:
            tr.data = tr.data * scale_factor

    stream.detrend('demean')

    st2 = add_sac_metadata(stream, ev=ev, stalist=stations)

    time_shift_sac(st2, -1 * before)

    if samp_freq != 0:
        resample(st2, freq=samp_freq)

    # Now do some QA: throw out traces with missing data
    # keep a log with status of trace extraction
    # this is mirrored in llnl_tool.py and iris_tools.py
    outlog = "get_data_status_IRIS.log"
    fid = open(outlog, "a")
    fid.write("\n--------------------\n%s\n" % ev.short_str())

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

    if rotate:
        rotate_and_write_stream(st2, evtime)

    if output_cap_weight_file:
        write_cap_weights(st2, evtime)

    if output_event_info:
        write_ev_info(ev, evtime)


def rotate_and_write_stream(stream, reftime):
    """
    Rotate an obspy stream to backazimuth

    Must contain event and station metadata in the sac header
    (stream.traces[0].stats.sac)

    stream - Obspy Stream() object
    reftime- UTCDateTime() object; this will be the output directory
        (default to current directory)
    """
    if reftime == '':
        outdir = './'
    else:
        outdir = './' + reftime.strftime('%Y%m%d%H%M%S%f')[:-3] + '/'
    print(outdir)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Directory is made, now rotate
    # Sorted stream makes for structured loop
    stream.sort()

    # Add backazimuth to metadata in a particular way...
    for tr in stream.traces:
        tr.stats.back_azimuth = tr.stats.sac['baz']
        print(tr.stats.station, ' backazimuth: ', tr.stats.back_azimuth)

        # Do some qc: throw out any stations without three components
        substr = stream.select(station=tr.stats.station)
        # 20160805 cralvizuri@alaska.edu -- currently this option
        # also skips if stream is not 3 (eg there is PFO.CI and PFO.TS -- i.e.
        # 6 streams). We want all possible data, so change this to < 3
        #if len(substr) != 3:
        if len(substr) < 3:
            for subtr in substr:
                print('One or more components missing: removing ',
                      subtr.stats.station, ' ', subtr.stats.channel,
                      ' Number of traces: ', len(substr))
                stream.remove(subtr)

    # Get list of unique stations
    stalist = []
    for tr in stream.traces:
        stalist.append(tr.stats.station)

    # Crazy way of getting a unique list of stations
    stalist = list(set(stalist))
    for station in stalist:
        substr = stream.select(station=station)
        try:
            substr.rotate('NE->RT')
        except:
            "Rotation failed, skipping..."
            continue

    # stream.rotate('NE->RT') #And then boom, obspy rotates everything!
    # Fix cmpaz metadata for Radial and Transverse components
    for tr in stream.traces:
        if tr.stats.channel[-1] == 'R':
            tr.stats.sac['cmpaz'] = tr.stats.sac['az']
        elif tr.stats.channel[-1] == 'T':
            tr.stats.sac['cmpaz'] = tr.stats.sac['az']+90.0
            if tr.stats.sac['cmpaz'] > 360.0:
                tr.stats.sac['cmpaz'] += -360
        # Now Write
        outfnam = outdir + tr.stats.station + '_' + tr.stats.network + '.' + \
            tr.stats.channel[-1].lower()
        tr.write(outfnam, format='SAC')


def write_cap_weights(stream, reftime=''):
    """
    Write CAP weight files from an Obspy stream

    Assumes stream has metadata added already

    reftime - a UTCDateTime Object
    """
    if reftime == '':
        outdir = './'
    else:
        outdir = './' + reftime.strftime('%Y%m%d%H%M%S%f')[:-3] + '/'
    outform = ('%12s %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f '
               '%4.0f %4.0f\n')
    laststa = ''
    #f = open(outdir + 'weight.dat', 'w')   # original (overwrites)
    # append instead of overwrite. needed when fetching data from multiple 
    # sources (IRIS, NCEDC)
    f = open(outdir + 'weight.dat', 'a')    

    # Sort the traces by distance
    sorted_traces = sorted(stream.traces, key=lambda k: k.stats.sac['dist'])
    for tr in sorted_traces:
        s = tr.stats
        # Only write once per 3 components
        if s.station == laststa:
            continue
        f.write(outform % (s.station + '_' + s.network, s.sac['dist'],
                1, 1, 1, 1, 1, 0, 0, 0, 0, 0))
        laststa = s.station


def write_ev_info(ev, reftime):
    if reftime == '':
        outdir = './'
    else:
        outdir = './'+reftime.strftime('%Y%m%d%H%M%S%f')[:-3]+'/'
    outform = '%s %f %f %f %f'
    f = open(outdir + 'ev_info.dat', 'w')
    f.write(outform % (
        ev.preferred_origin().time, ev.preferred_origin().longitude,
        ev.preferred_origin().latitude, ev.preferred_origin().depth / 1000.0,
        ev.magnitudes[0].mag))


def add_sac_metadata(st, ev=[], stalist=[]):
    """
    Add event and station metadata to an Obspy stream
    """
    # Loop over each trace
    for tr in st.traces:
        # Write each one
        tr.write('tmppp.sac', format='SAC')
        # This is the best way to make sac objects
        tmptr = obspy.read('tmppp.sac').traces[0]
        # Loop over all the networks
        for net in stalist:
            # Find the right station
            for stan in net:
                # Hopefully there isn't more than one
                if tr.stats.station == stan.code:
                    sta = stan
        tr.stats.sac = tmptr.stats.sac

        # Station info
        tr.stats.sac['stla'] = sta.latitude
        tr.stats.sac['stlo'] = sta.longitude
        # Change to kilometers
        tr.stats.sac['stel'] = sta.elevation / 1000.0

        # Event info
        tr.stats.sac['evla'] = ev.preferred_origin().latitude
        tr.stats.sac['evlo'] = ev.preferred_origin().longitude
        tr.stats.sac['evdp'] = ev.preferred_origin().depth/1000
        tr.stats.sac['mag'] = ev.preferred_magnitude().mag

        # !!!!Weird!!!
        tr.stats.sac['kevnm'] = \
            ev.preferred_origin().time.strftime('%Y%m%d%H%M%S%f')[:-3]

        # Station-event info
        tr.stats.sac['dist'], tr.stats.sac['az'], tr.stats.sac['baz'] = \
            obspy.core.util.gps2DistAzimuth(
                tr.stats.sac['evla'], tr.stats.sac['evlo'],
                tr.stats.sac['stla'], tr.stats.sac['stlo'])

        tr.stats.back_azimuth = tr.stats.sac['baz']
        # Is this right?
        tr.stats.sac['gcarc'] = tr.stats.sac.dist * 111.19

        # Kilometers
        tr.stats.sac['dist'] = tr.stats.sac['dist'] / 1000

        # Now add component info. Fucking obspy should do this but doesn't
        tmp = tr.stats.channel[2]
        tr.stats.sac['cmpinc'] = 90.0
        tr.stats.sac['cmpaz'] = 0.0

        if tmp == 'E':
            tr.stats.sac['cmpaz'] = 90.0
        elif tmp == 'Z':
            tr.stats.sac['cmpinc'] = 0.0

        # Now some extra things
        tr.stats.sac['o'] = 0
        # Component isn't left-handed?
        tr.stats.sac['lpspol'] = 0
        # I don't even know
        tr.stats.sac['lcalda'] = 1

        print(tr.stats.station, tr.stats.sac['evlo'], tr.stats.sac['evla'],
              tr.stats.sac['stlo'], tr.stats.sac['stla'], tr.stats.sac['dist'])
    return st


def time_shift_sac(st, tshift=0):
    """
    Shift the b and e values in the sac header
    """
    for tr in st.traces:
        tr.stats.sac['b'] = tr.stats.sac['b'] + tshift
        tr.stats.sac['e'] = tr.stats.sac['e'] + tshift


class Stalist(list):
    """
    A class that is just a list of station names for now
    """
    def __init__(self, fnam):
        # fnam = a filename for an ascii file with a list of station names
        fid = open(fnam, 'r')
        a = fid.readlines()
        for line in a:
            # Assumes one name per line
            self.append(line.strip().split()[0])

    def make_bulk_list(self, t1, t2, net="*", loc="*", cmp="BH*"):
        # Create a list for client.Client.get_stations_bulk or
        # get_waveforms_bulk
        # t1 and t2 are in UTCDateTime format
        self.bulk_list = []
        for sta in self:
            self.bulk_list.append((net, sta, loc, cmp, t1, t2))


def make_bulk_list_from_stalist(stations, t1, t2, loc='*', cmp='BH*'):
    bulk_list = []
    for net in stations:
        for sta in net:
            bulk_list.append((net.code, sta.code, loc, cmp, t1, t2))
    return bulk_list


def sta_limit_distance(ev, stations, min_dist=0, max_dist=100000,
                       ifverbose=False):
    # Remove stations greater than a certain distance from event
    elat = ev.preferred_origin().latitude
    elon = ev.preferred_origin().longitude
    remlist = []

    # Loop over network and stations
    for net in stations:
        for sta in net:
            dist, az, baz = obspy.core.util.gps2DistAzimuth(
                elat, elon, sta.latitude, sta.longitude)

            # Discard station if too far
            if dist / 1000.0 < min_dist or dist / 1000.0 > max_dist:
                # Keep a list of what to discard.
                remlist.append(sta)
        # Loop over discard list
        for sta in remlist:
            # Not sure why this is needed, but it is.
            if sta in net.stations:
                net.stations.remove(sta)
    if ifverbose:
        for net in stations:
            for sta in net:
                dist, az, baz = obspy.core.util.gps2DistAzimuth(
                    elat, elon, sta. latitude, sta.longitude)
                print(sta.code, elon, elat, sta.longitude, sta.latitude,
                      dist / 1000)


def write_stream_sac(st, reftime='', odir='./', use_sta_as_dirname=False):
    """
    Writes out all of the traces in a stream to sac files

    Naming convention - DATE.NET.STA.CMP.SAC
    """
    for tr in st.traces:
        if reftime == '':
            usetime = tr.stats.starttime
        else:
            usetime = reftime
        outdir = odir+usetime.strftime('%Y%m%d%H%M%S%f')[:-3] + '/'

        if use_sta_as_dirname:
            outdir = odir + tr.stats.station + '/'

        print(outdir)

        if not os.path.exists(outdir):
            os.makedirs(outdir)
            print(outdir, tr.stats.station)

        outfnam = outdir + usetime.strftime('%Y%m%d%H%M%S%f')[:-3] + '.' + \
            tr.stats.station + '.' + tr.stats.channel + '.' + \
            tr.stats.network + '.sac'
        tr.write(outfnam, format='SAC')
