#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Utilities to prepare files for CAP
For reference see versions prior to Aug 25, 2016 for:
    getwaveform_iris.py
    getwaveform_llnl.py

20160825 cralvizuri <cralvizuri@alaska.edu>
"""

import obspy
from obspy.io.sac import SACTrace
import obspy.signal.rotate as rt
import os
from scipy import signal
import numpy as np
import util_helpers

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

def rotate_and_write_stream(stream, evname_key):
    """
    Rotate an obspy stream to backazimuth

    Must contain event and station metadata in the sac header
    (stream.traces[0].stats.sac)

    stream - Obspy Stream() object
    reftime- UTCDateTime() object; this will be the output directory
        (default to current directory)
    """

    # pending removing
    #if reftime == '':
    #    outdir = './'
    #else:
    #    outdir = './' + reftime.strftime('%Y%m%d%H%M%S%f')[:-3] + '/'
    #print(outdir)
    #

    # get name key for output directory and files
    outdir = evname_key                  

    if not os.path.exists(outdir):
       os.makedirs(outdir)

    # Directory is made, now rotate
    # Sorted stream makes for structured loop
    stream.sort()

    # Add backazimuth to metadata in a particular way...
    for tr in stream.traces:
        tr.stats.back_azimuth = tr.stats.sac['baz']
        # code for debugging. print station and backazimuth
        #print(tr.stats.station, ' backazimuth: ', tr.stats.back_azimuth)   

        # Do some qc: throw out any stations without three components
        substr = stream.select(station=tr.stats.station)

        # Select the alphabetically lowest band code.
        band_code = sorted([tr.stats.channel[0] for tr in substr])[0]
        substr = substr.select(channel=band_code + "*")

#        # 20160805 cralvizuri@alaska.edu -- currently this option
#        # also skips if stream is not 3 (eg there is PFO.CI and PFO.TS -- i.e.
#        # 6 streams). We want all possible data, so change this to < 3
#        #if len(substr) != 3:
        if len(substr) < 3:
            for subtr in substr:
                print('One or more components missing: removing ',
                      subtr.stats.station, ' ', subtr.stats.channel,
                      ' Number of traces: ', len(substr))
                #stream.remove(subtr)

    # Get list of unique stations + locaiton (example: 'KDAK.00')
    stalist = []
    for tr in stream.traces:
        #stalist.append(tr.stats.station)
        stalist.append(tr.stats.network + '.' + tr.stats.station +'.'+ tr.stats.location + '.'+ tr.stats.channel[:-1])
        
    # Initialize stream object
    # For storing extra traces in case there are less than 3 compnents   
    st_new = obspy.Stream()

    # Crazy way of getting a unique list of stations
    stalist = list(set(stalist))
    #print(stalist)
    for stn in stalist:
        # split STNM.LOC
        tmp = stn.split('.')
        netw = tmp[0]
        station = tmp[1]
        location = tmp[2]
        chan = tmp[3] + '*'
        # Get 3 traces (subset based on matching station name and location code)
        substr = stream.select(network=netw,station=station,location=location,channel=chan)
        #print(substr)   # code for debugging. print stream information
        if len(substr) < 3:
            print('WARNING:', len(substr), 'traces available for rotation. Adding NULL traces - ', \
                      netw + '.' + station + '.' + location + '.' + chan)
            d1 = substr[0].data
            dip1 = substr[0].stats.sac['cmpinc']
            az1 = substr[0].stats.sac['cmpaz']
            for itr in range(len(substr),3):
                tmp = substr[0].copy()
                #print(tmp.data)
                substr.append(tmp)
                substr[itr].data = np.zeros(substr[0].stats.npts)
                substr[itr].stats.sac['cmpinc'] = dip1 + 90.0
                tmp = None
                if itr == 1:
                    substr[itr].stats.sac['cmpaz'] = az1
                if itr == 2:
                    substr[itr].stats.sac['cmpaz'] = az1 + 90.0

        # Rotate to NEZ first
        # Sometimes channels are not orthogonal (example: 12Z instead of NEZ)
        # Run iex = 5 (run_getwaveform.py) for one such case       
        d1 = substr[0].data
        d2 = substr[1].data
        d3 = substr[2].data
        az1 = substr[0].stats.sac['cmpaz']
        az2 = substr[1].stats.sac['cmpaz']
        az3 = substr[2].stats.sac['cmpaz']
        dip1 = substr[0].stats.sac['cmpinc']
        dip2 = substr[1].stats.sac['cmpinc']
        dip3 = substr[2].stats.sac['cmpinc']
        #print('=R==>',substr[0].stats.channel, substr[0].stats.sac['cmpinc'], substr[0].stats.sac['cmpaz'], \
        #          substr[1].stats.channel, substr[1].stats.sac['cmpinc'],substr[1].stats.sac['cmpaz'], \
        #          substr[2].stats.channel, substr[2].stats.sac['cmpinc'], substr[2].stats.sac['cmpaz'])
        print('--> Station ' + netw + '.' + station + '.' + location + \
                  ' Rotating random orientation to NEZ.')
        data_array = rt.rotate2zne(d1, az1, dip1, d2, az2, dip2, d3, az3, dip3)
        # Rotates an arbitrarily oriented three-component vector to ZNE( [0]-Z, [1]-N, [2]-E)
        # XXX: Check 012 in correct order? 
        substr[0].data =  data_array[2]  # E
        substr[1].data =  data_array[1]  # N
        substr[2].data =  data_array[0]  # Z

        # Fix the channel names in the traces.stats
        substr[0].stats.channel = substr[0].stats.channel[0:2] + 'E'
        substr[1].stats.channel = substr[0].stats.channel[0:2] + 'N'
        substr[2].stats.channel = substr[0].stats.channel[0:2] + 'Z'
        # Fix the sac headers since the traces have been rotated now
        substr[0].stats.sac['cmpaz'] = 90.0
        substr[1].stats.sac['cmpaz'] = 0.0
        substr[2].stats.sac['cmpaz'] = 0.0
        # matlab files had following cmpinc: E = 90, N = 90, Z = 0
        # XXX: Will this cause problem??
        substr[0].stats.sac['cmpinc'] = 0.0
        substr[1].stats.sac['cmpinc'] = 0.0
        substr[2].stats.sac['cmpinc'] = -90.0
        # Fix sac headers
        substr[0].stats.sac['kcmpnm'] = substr[0].stats.channel
        substr[1].stats.sac['kcmpnm'] = substr[1].stats.channel
        substr[2].stats.sac['kcmpnm'] = substr[2].stats.channel
 
        # save NEZ waveforms
        for tr in substr:
            #outfnam = outdir + reftime.strftime('%Y%m%d%H%M%S%f')[:-3] + '.' \
            outfnam = outdir + '/' + evname_key + '.' \
                + tr.stats.network + '.' + tr.stats.station + '.' \
                + tr.stats.location + '.' + tr.stats.channel[:-1] + '.' \
                + tr.stats.channel[-1].lower()
            tr.write(outfnam, format='SAC')

        try:
            print('--> Station ' + netw + '.' + station + '.' + location + \
                      ' Rotating ENZ to RTZ.')
            substr.rotate('NE->RT')
            #print('=R==>',substr[0].stats.channel, substr[0].stats.sac['cmpinc'], substr[0].stats.sac['cmpaz'], \
            #      substr[1].stats.channel, substr[1].stats.sac['cmpinc'],substr[1].stats.sac['cmpaz'], \
            #      substr[2].stats.channel, substr[2].stats.sac['cmpinc'], substr[2].stats.sac['cmpaz'])
        except:
            "Rotation failed, skipping..."
            continue

        # append substream to the main stream
        st_new = st_new + substr

    # replace stream object
    stream = st_new

    # stream.rotate('NE->RT') #And then boom, obspy rotates everything!
    # Fix cmpaz metadata for Radial and Transverse components
    for tr in stream.traces:
        if tr.stats.channel[-1] == 'R':
            tr.stats.sac['kcmpnm'] = tr.stats.channel[0:2] + 'R'
            tr.stats.sac['cmpaz'] = tr.stats.sac['az']
        elif tr.stats.channel[-1] == 'T':
            tr.stats.sac['kcmpnm'] = tr.stats.channel[0:2] + 'T'
            tr.stats.sac['cmpaz'] = tr.stats.sac['az']+90.0
            if tr.stats.sac['cmpaz'] > 360.0:
                tr.stats.sac['cmpaz'] += -360

        # Now Write
        # 20160805 cralvizuri@alaska.edu -- some llnl stations have traces with
        # multiple channel types. The original filename does not include 
        # channel, so they're overwritten. eg KNB_BB and KNB_HF both are output
        # as KNB_[component], so one is lost. this change fixes that.
        #outfnam = outdir + tr.stats.station + '_' + tr.stats.network + '.' + \
        #    tr.stats.channel[-1].lower()
        #outfnam = outdir + reftime.strftime('%Y%m%d%H%M%S%f')[:-3] + '.' \
        outfnam = outdir + '/' + evname_key + '.' \
                + tr.stats.network + '.' + tr.stats.station + '.' \
                + tr.stats.location + '.' + tr.stats.channel[:-1] + '.' \
                + tr.stats.channel[-1].lower()
        tr.write(outfnam, format='SAC')

def write_cap_weights(stream, evname_key, client_name='', event=''):
    """
    Write CAP weight files from an Obspy stream

    Assumes stream has metadata added already

    reftime - a UTCDateTime Object
    """
    # For debugging LLNL example
    # if event != '':
    #    for pick in event.picks:
    #        print(pick.waveform_id.network_code, pick.waveform_id.station_code, pick.waveform_id.channel_code, \
    #                  pick.waveform_id.location_code, pick.time, pick.phase_hint)

    # get name key for output directory and files
    outdir = evname_key

    #if reftime == '':
    #    outdir = './'
    #else:
    #    outdir = './' + reftime.strftime('%Y%m%d%H%M%S%f')[:-3] + '/'


    outform = ('%35s %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.3f %4.0f %4.0f '
               '%4.0f %4.0f\n')

    laststa = ''
    lastloc = ''
    lastcha = ''

    #f = open(outdir + 'weight.dat', 'w')   # original (overwrites weightfile)
    # append instead of overwrite. This is needed when fetching data from
    # multiple sources (IRIS, NCEDC). Else weight files are overwritten.
    # weight.dat    - append. includes all stations from all clients
    # weight_body_XX.dat  - write a weight (body waves) file for client XX
    # weight_surf_XX.dat  - write a weight (surf waves) file for client XX
    wfile = outdir + '/' + "weight.dat"
    wfile_body = outdir + '/' + "weight_body" + ".dat"
    wfile_surf = outdir + '/' + "weight_surf" + ".dat"
    wfile_body_client = outdir + '/' + "weight_body" + "_" + client_name + ".dat"
    wfile_surf_client = outdir + '/' + "weight_surf" + "_" + client_name + ".dat"
    f =  open(wfile, 'a') 
    fb = open(wfile_body, 'a') 
    fs = open(wfile_surf, 'a') 
    fbc = open(wfile_body_client, 'w') 
    fsc = open(wfile_surf_client, 'w') 

    # Sort the traces by distance
    sorted_traces = sorted(stream.traces, key=lambda k: k.stats.sac['dist'])
    for tr in sorted_traces:
        current_sta = tr.stats
        # Only write once per 3 components
        if (current_sta.station == laststa) and (current_sta.channel[:-1] == lastcha) and (current_sta.location == lastloc):
            continue

        Ptime = 0
        for pick in event.picks:
            if (pick.waveform_id.network_code == tr.stats.network and pick.waveform_id.station_code == tr.stats.station \
                    and pick.waveform_id.channel_code[2].upper() == 'Z' and pick.waveform_id.location_code == tr.stats.location and pick.phase_hint == 'Pn'):
                # For debugging
                # print(pick.waveform_id.network_code, pick.waveform_id.station_code, pick.waveform_id.channel_code, \
                #          pick.waveform_id.location_code, pick.time, pick.phase_hint,tr.stats.channel, tr.stats.location )
                Ptime = pick.time - event.origins[0].time
        
        #outfnam = reftime.strftime('%Y%m%d%H%M%S%f')[:-3] + '.' \
        outfnam = outdir + '/' + evname_key + '.' \
                + tr.stats.network + '.' + tr.stats.station + '.' \
                + tr.stats.location + '.' + tr.stats.channel[:-1]

        f.write(outform % (outfnam, current_sta.sac['dist'], 
            1, 1, 1, 1, 1, Ptime, 0, 0, 0, 0))
        fb.write(outform % (outfnam, current_sta.sac['dist'], 
            1, 1, 0, 0, 0, Ptime, 0, 0, 0, 0))
        fs.write(outform % (outfnam, current_sta.sac['dist'], 
            0, 0, 1, 1, 1, Ptime, 0, 0, 0, 0))
        fbc.write(outform % (outfnam, current_sta.sac['dist'], 
            1, 1, 0, 0, 0, Ptime, 0, 0, 0, 0))
        fsc.write(outform % (outfnam, current_sta.sac['dist'], 
            0, 0, 1, 1, 1, Ptime, 0, 0, 0, 0))

        laststa = current_sta.station
        lastloc = current_sta.location
        lastcha = current_sta.channel[:-1]

def write_ev_info(ev, evname_key):
    outdir = evname_key
    #if reftime == '':
    #    outdir = './'
    #else:
    #    outdir = './'+reftime.strftime('%Y%m%d%H%M%S%f')[:-3]+'/'

    fout_event_info = outdir + '/' + evname_key + "_event_info.dat"
    outform = '%s %f %f %f %f'

    f = open(fout_event_info, 'w')
    f.write(outform % (
        ev.origins[0].time, ev.origins[0].longitude,
        ev.origins[0].latitude, ev.origins[0].depth / 1000.0,
        ev.magnitudes[0].mag))

def add_sac_metadata(st, ev=[], stalist=[]):
    """
    Add event and station metadata to an Obspy stream
    """
    fid = open('traces_inventory_log', "w")
    out_form = ('%s %s %s %s %s %s %s %s %s')
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
        tr.stats.sac['evla'] = ev.origins[0].latitude
        tr.stats.sac['evlo'] = ev.origins[0].longitude
        tr.stats.sac['evdp'] = ev.origins[0].depth/1000
        m = ev.preferred_magnitude() or ev.magnitudes[0]
        tr.stats.sac['mag'] = m.mag

        # Add P arrival time
        for pick in ev.picks:
            if (pick.waveform_id.network_code == tr.stats.network and pick.waveform_id.station_code == tr.stats.station \
                    and pick.waveform_id.channel_code[2].upper() == 'Z' and pick.waveform_id.location_code == tr.stats.location and pick.phase_hint == 'Pn'):
                Ptime = pick.time - ev.origins[0].time
                tr.stats.sac['a'] = Ptime

        # !!!!Weird!!!
        tr.stats.sac['kevnm'] = \
            ev.origins[0].time.strftime('%Y%m%d%H%M%S%f')[:-3]

        # Station-event info
        tr.stats.sac['dist'], tr.stats.sac['az'], tr.stats.sac['baz'] = \
            obspy.geodetics.gps2dist_azimuth(
                tr.stats.sac['evla'], tr.stats.sac['evlo'],
                tr.stats.sac['stla'], tr.stats.sac['stlo'])

        tr.stats.back_azimuth = tr.stats.sac['baz']
        # Is this right?
        tr.stats.sac['gcarc'] = tr.stats.sac.dist * 111.19

        # Kilometers
        tr.stats.sac['dist'] = tr.stats.sac['dist'] / 1000

        # Now add component info. CMPAZ and CMPINC info is in the station inventor
        # and not in the trace.stats
        tmp = tr.stats.channel
        
        # match trace station info with the station inventory info
        for net in stalist:
            for sta in net.stations:
                for ch in sta.channels:
                    # code for debugging. print stations names from traces and inventory
                    # fid.write(out_form % ('\n--->', tr.stats.channel, ch.code, tr.stats.location, ch.location_code, tr.stats.station, sta.code, tr.stats.network, net.code))
                    if tr.stats.channel == ch.code.upper() and \
                            tr.stats.location == ch.location_code and \
                            tr.stats.station == sta.code and \
                            tr.stats.network == net.code:
                        # fid.write(out_form % ('\n--->', tr.stats.channel, ch.code, tr.stats.location, ch.location_code, tr.stats.station, sta.code, tr.stats.network, net.code))
                        # code for debugging. print azimuth and dip
                        # print('--->', net.code, sta.code, ch.location_code, ch.code, 'Azimuth:', ch.azimuth, 'Dip:', ch.dip) 
                        tr.stats.sac['cmpinc'] = ch.dip
                        tr.stats.sac['cmpaz'] = ch.azimuth
                
        # obspy has cmpinc for Z component as -90
        #if tmp == 'Z':
        #    tr.stats.sac['cmpinc'] = 0.0
      
        #-------- obsolete (remove?) ------
        # We should not be assuming the component azimuth=0 and cmpinc=90
        #tr.stats.sac['cmpinc'] = 90.0
        #tr.stats.sac['cmpaz'] = 0.0
            
        #if tmp == 'E':
        #    tr.stats.sac['cmpaz'] = 90.0
        #elif tmp == 'Z':
        #    tr.stats.sac['cmpinc'] = 0.0
        #--------------------------------

        # Now some extra things
        tr.stats.sac['o'] = 0
        # Component isn't left-handed?
        tr.stats.sac['lpspol'] = 0
        # I don't even know
        tr.stats.sac['lcalda'] = 1

        # code for debugging.
        #print(tr.stats.station, tr.stats.sac['evlo'], tr.stats.sac['evla'],
        #      tr.stats.sac['stlo'], tr.stats.sac['stla'], tr.stats.sac['dist'])
    return st

def check_if_LLNL_event(event_time):
    """
    Update SAC header kevnm if current event is in the Ford paper
    """
    sec_threshold = 5

    # event times and names are from the LLNL database
    event_time_name_LLNL = {
            # nuclear tests
            "1988-02-15T18:10:00.090": "KERNVILLE", 
            "1989-06-27T15:30:00.020": "AMARILLO", 
            "1989-09-14T15:00:00.100": "DISKO_ELM", 
            "1989-10-31T15:30:00.090": "HORNITOS", 
            "1989-12-08T15:00:00.090": "BARNWELL", 
            "1990-03-10T16:00:00.080": "METROPOLIS", 
            "1990-06-13T16:00:00.010": "BULLION", 
            "1990-06-21T18:15:00.040": "AUSTIN", 
            "1990-11-14T19:17:00.710": "HOUSTON", 
            "1991-03-08T21:02:45.080": "COSO_SILVER", 
            "1991-04-04T19:00:00.000": "BEXAR", 
            "1991-09-14T19:00:00.050": "HOYA", 
            "1991-10-18T19:12:00.000": "LUBBOCK", 
            "1991-11-26T18:34:59.670": "BRISTOL", 
            "1992-03-26T16:30:00.000": "JUNCTION", 
            "1992-09-18T17:00:00.010": "HUNTERS_TROPHY", 
            "1992-09-23T15:04:00.000": "DIVIDER",
            # mine collapses
            "1982-08-05T14:00:00.090": "ATRISCO_Hole",
            "1995-02-03T15:26:12.870": "Trona_Mine_1",
            "2000-01-30T14:46:53.330": "Trona_Mine_2",
            # earthquakes
            "1992-06-29T10:14:20.060": "Little_Skull_Main",
            "1992-07-05T06:54:13.080": "Little_Skull_Aftershock",
            "1995-07-31T12:34:47.050": "Timber_Mountain",
            "1997-04-26T01:49:35.440": "Groom_Pass",
            "1997-09-12T13:36:54.650": "Calico_Fan",
            "1998-12-12T01:41:31.980": "Warm_Springs",
            "1999-01-23T03:00:32.370": "Frenchman_Flat1",
            "1999-01-27T10:44:23.300": "Frenchman_Flat2",
            "2002-06-14T12:40:44.450": "Little_Skull",
            # earthquakes in Ford but not in the LLNL database
            "1996-09-05T08:16:56.09": "Amargosa",
            "1997-06-14T19:48:19.93": "Indian_Springs",
            "2007-01-24T11:30:16.10": "Raslton"
            }

    _event_time = obspy.UTCDateTime(event_time)
    evname_key = util_helpers.otime2eid(event_time)     # object
    is_an_llnl_event = False
    for llnl_evtime, evname in event_time_name_LLNL.items():
        # update all headers kevnm if this is an LLNL event
        if abs(_event_time - obspy.UTCDateTime(llnl_evtime)) <= sec_threshold:
            evname_key = evname         # string
            is_an_llnl_event = True
            break

    return evname_key, is_an_llnl_event

def rename_if_LLNL_event(st, event_time):
    """
    Update SAC header kevnm if current event is in the Ford paper
    """
    evname_key, is_an_llnl_event = check_if_LLNL_event(event_time)

    if(is_an_llnl_event):
        print("--> WARNING. This is an LLNL event. " +\
                "New event name: " + evname_key)
        for tr in st.traces:
            tr.stats.sac['kevnm'] = evname_key
    else:
        evname_key = st[0].stats.sac['kevnm']

    return st, evname_key

def time_shift_sac(st, tshift=0):
    """
    Shift the b and e values in the sac header
    """
    for tr in st.traces:
        tr.stats.sac['b'] = tr.stats.sac['b'] + tshift
        tr.stats.sac['e'] = tr.stats.sac['e'] + tshift

def correct_sac_tshift(targetdir,before=100.0,after=300.0):
    """
    Change the b and e header values of every sac file by 
    before and after (floats)

    NOTE
    - requires SAC
    """
    import os
    os.chdir(targetdir)
    fnam='sac_cmd'
    ff=open(fnam,'w')
    ff.write('r *.r *.t *.z *.sac;\n')
    ff.write('ch b '+str(-1*before)+' e '+str(after)+'\n')
    #ff.write('w over\nquit')
    ff.write('w over;\n')
    ff.write('quit;\n')
    ff.close()
    os.system('sac < sac_cmd;')
    os.chdir('../')

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

    def make_bulk_list(self, t1, t2, net="*", loc="*", channel="BH*"):
        # Create a list for client.Client.get_stations_bulk or
        # get_waveforms_bulk
        # t1 and t2 are in UTCDateTime format
        self.bulk_list = []
        for sta in self:
            self.bulk_list.append((net, sta, loc, channel, t1, t2))

def make_bulk_list_from_stalist(stations, t1, t2, loc='*', channel='BH*'):
    bulk_list = []
    for net in stations:
        for sta in net:
            bulk_list.append((net.code, sta.code, loc, channel, t1, t2))
    return bulk_list

def sta_limit_distance(ev, stations, min_dist=0, max_dist=100000,
                       ifverbose=True):
    # Remove stations greater than a certain distance from event
    elat = ev.origins[0].latitude
    elon = ev.origins[0].longitude
    reftime = ev.origins[0].time
    remlist = []
    rename_if_LLNL_event

    #outdir = './' + reftime.strftime('%Y%m%d%H%M%S%f')[:-3] + '/'
    evname_key, is_an_llnl_event = check_if_LLNL_event(reftime)
    outdir = evname_key
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outfile = outdir + '/' + evname_key + "_station_list_ALL.dat"
    f = open(outfile, 'w') 
    outform = '%s %s %f %f %f %f\n'

    # Loop over network and stations
    for net in stations:
        for sta in net:
            dist, az, baz = obspy.geodetics.gps2dist_azimuth(
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

    # write the full list of stations to file
    if ifverbose:
        for net in stations:
            for sta in net:
                dist, az, baz = obspy.geodetics.gps2dist_azimuth(
                    elat, elon, sta.latitude, sta.longitude)
                print(sta.code, elon, elat, sta.longitude, sta.latitude,
                      dist / 1000)
                f.write(outform % (sta.code, net.code, sta.latitude, sta.longitude, dist / 1000, az))

#def write_stream_sac(st, reftime='', odir='./', use_sta_as_dirname=False):
def write_stream_sac(st, evname_key):
    """
    Writes out all of the traces in a stream to sac files

    Naming convention - DATE.NET.STA.CMP.SAC
    """

    #if reftime == '':
    #    usetime = st[0].stats.starttime
    #else:
    #    usetime = reftime

    # get name key for output directory and files
    outdir = evname_key

    # this smells like debug code. disabling for now.
    #if use_sta_as_dirname:
    #    outdir = odir + tr.stats.station + '/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        #print(outdir, tr.stats.station)

    # write sac waveforms
    for tr in st.traces:
        #outfnam = outdir + usetime.strftime('%Y%m%d%H%M%S%f')[:-3] + '.' \
        outfnam = outdir + "/" + evname_key + '.' \
                + tr.stats.network + '.' + tr.stats.station + '.' \
                + tr.stats.location + '.' + tr.stats.channel + '.sac'
        tr.write(outfnam, format='SAC')

def set_reftime(stream, evtime):
    """
    Set a reftime for all traces. For CAP we currently set the reftime as 
    the event time.
    """
    stream2 = obspy.Stream()
    for tr in stream:
        sac = SACTrace.from_obspy_trace(tr)
        sac.reftime = evtime
        sac.o = 0
        tr = sac.to_obspy_trace()
        stream2.append(tr)
    return(stream2)

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
        try:
            tr.detrend("linear")
        except:
            print("WARNING. Unable to detrend for " + tr.stats.network \
                    + '.' + tr.stats.station + '.' + tr.stats.channel + '. '\
                    "Waveforms may contain NaN or Inf values")
            # XXX REMOVE TRACE? TR.REMOVE(...
            continue
        tr.taper(max_percentage=0.02, type="hann")
        tr.data = np.require(tr.data, requirements=["C"])
        tr.interpolate(sampling_rate=freq, method="lanczos", a=8,
                       window="blackman")

def resample_cut(st, freq, evtime, before, after):
    """
    Custom resampling with a very sharp zerophase filter.
    """
    new_nyquist = 0.5 * freq
    for tr in st:
        current_nyquist = 0.5 * tr.stats.sampling_rate
        # Filter if necessary.
        if new_nyquist < current_nyquist:
            zerophase_chebychev_lowpass_filter(trace=tr, freqmax=new_nyquist)
        try:
            tr.detrend("linear")
        except:
            print("WARNING. Rejecting station %s. Unable to detrend. " % \
                    tr.stats.station)
            print("This may be due to NaN and Inf values in the data.")
            # XXX REMOVE TRACE? TR.REMOVE(...
            continue
        tr.taper(max_percentage=0.02, type="hann")
        tr.data = np.require(tr.data, requirements=["C"])
        if (tr.stats.starttime > evtime - before):
            print("WARNING. trace starttime > otime-before for station %s" % \
                    tr.stats.station)
            before = evtime - tr.stats.starttime # assuming starttime > evtime!
            print("Setting before = evtime - starttime = %f" % before)
        if (tr.stats.endtime < evtime + after):
            print("WARNING. trace endtime < otime+after for station %s" % \
                    tr.stats.station)
            print("Setting after = evtime + endtime = %f" % after)

            # Needed (eg LLNL station KNB). see iex=10 in run_getwaveform.py
            try:
                after = evtime + tr.stats.endtime
            except:
                print("Removing %s: Unable to calculate endtime. Continuing"\
                        % tr.stats.station)
                st.remove(tr)
                continue

        starttime = evtime - before
        npts = (before + after) * freq
        tr.interpolate(sampling_rate=freq, method="lanczos", a=8,
                       window="blackman", starttime = starttime, npts = npts)

def trim_maxstart_minend(stalist, st2):
    """
    Function to cut the start and end points of a stream.
    The starttime and endtime are set to the latest starttime and earliest
    endtime from all channels (up to 3 channels).
    """
    print("---> trimming end points")
    temp_stream = obspy.Stream()

    # Trim the edges in case 3 channels have different lengths
    for stn in stalist:
        # split station codes
        split_stream = stn.split('.')
        netw = split_stream[0]
        station = split_stream[1]
        location = split_stream[2]
        chan = split_stream[3] + '*'

        # Get 3 traces (subset based on matching station name and location code)
        #temp_stream = stream.select(network=netw,station=station,location=location,channel=chan)
        select_st = st2.select(network = netw, station = station, \
                location = location, channel = chan)

        # Find max startime, min endtime for stations with 1 or 2 or 3 channels
        if len(select_st) == 1:
            max_starttime = select_st[0].stats.starttime
            min_endtime = select_st[0].stats.endtime
        if len(select_st) == 2:
            max_starttime = max(select_st[0].stats.starttime,\
                    select_st[1].stats.starttime)
            min_endtime = min(select_st[0].stats.endtime,\
                    select_st[1].stats.endtime)
        if len(select_st) == 3:
            max_starttime = max(select_st[0].stats.starttime,\
                    select_st[1].stats.starttime,\
                    select_st[2].stats.starttime)
            min_endtime = min(select_st[0].stats.endtime,\
                    select_st[1].stats.endtime,\
                    select_st[2].stats.endtime)    

        try:
            select_st.trim(starttime=max_starttime,\
                               endtime = min_endtime, pad = False, nearest_sample = True,\
                               fill_value=0)
            for tr in select_st.traces:
                temp_stream = temp_stream.append(tr)
        except:
            print('WARNING: stattime larger than endtime for channels of', netw, '.', station, '.', location)
        
    return temp_stream

def amp_rescale_llnl(st, scale_factor):
    """
    Rescale amplitudes for LLNL data. The scales are different for the 
    BB and HF channels
    """
    # scales based on tests with hoya event
    scale_factor_BB = 4.0e-4 
    scale_factor_HF = 1.8e-2

    for tr in st.traces:
        if ('BB' in tr.stats.channel):
            tr.data = tr.data * scale_factor * scale_factor_BB
            print("--> WARNING -- " + tr.stats.station + '.' + \
                    tr.stats.channel + " rescaling by " + "%f" % scale_factor_BB)
        elif ('HF' in tr.stats.channel):
            tr.data = tr.data * scale_factor * scale_factor_HF
            print("--> WARNING -- " + tr.stats.station + '.' + \
                    tr.stats.channel + " rescaling by " + "%f" % scale_factor_HF)
        else:
            # original -- apply a single factor to all
            tr.data = tr.data * scale_factor
            print("\n--> WARNING -- rescaling amplitudes by %f" % scale_factor)

    return st

