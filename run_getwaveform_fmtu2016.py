#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to get waveform data for events at the Nevada Test Site. 
This data is for the following paper

C. Alvizuri, V. Silwal, L. Krischer, and C. Tape. Estimation of full moment
tensors with uncertainties, for earthquakes, volcanic events, and nuclear
tests. Geophysics (in prep.).

Data sources:
    LLNL
    IRIS
    NCEDC

Data references:
    Walter et al (2004) (LLNL dataset)
    Ford et al (2009)   (event selection)

20160801 cralvizuri <cralvizuri@alaska.edu>
"""

import obspy
import getwaveform_iris
import getwaveform_llnl
import llnl_db_client
from obspy.clients.fdsn import Client

# get event data from iris catalogs
def get_iris_evcat(t0):
    c = Client("IRIS")  # also can use NCEDC
    cat = c.get_events(starttime = t0-10, endtime = t0+10)
    return cat[0]

# get data  
def get_data_iris_ncedc(cat0):
    print("Calling getwaveform_iris")

    # parameters for waveform request
    min_dist = 0 
    max_dist = 1200
    before = 100
    after = 1200
    channel = 'EH*'
    resample_freq = 0.0                      # 0 for no resampling

    # parameters for CAP
    rotate = True
    output_cap_weight_file = True
    remove_response = True
    detrend = True
    demean = True
    output_event_info = True
    pre_filt = (0.005, 0.006, 5.0, 10.0)
    scale_factor = 10**2    # for CAP use 10**2 (to convert m/s to cm/s)

    # BK data has to be requested through NCEDC
    clients = ["IRIS", "NCEDC"]
    #clients = ["IRIS"]
    for thisclient in clients:
        if thisclient is "IRIS":
            network='*'
            # network = 'CI,TS,II,IM,IU' 
        elif thisclient is "NCEDC":
            network = 'BK'
        else:
            print("error. client not recognized")

        client = Client(thisclient)
        try:
            getwaveform_iris.run_get_waveform(c = client, event = cat0, 
                    min_dist = min_dist, max_dist = max_dist, 
                    before = before, after = after, 
                    network = network, channel = channel, 
                    resample_freq = resample_freq, 
                    ifrotate = rotate, ifCapInp = output_cap_weight_file, 
                    ifRemoveResponse = remove_response, ifDetrend = detrend, 
                    ifDemean = demean, ifEvInfo = output_event_info, 
                    scale_factor = scale_factor, pre_filt = pre_filt)
        except Exception as msg:
            print("WARNING. There was a problem with client %s" % thisclient)
            raise(msg)
            print("Continuing with next client\n")
            pass

def get_data_llnl(evid):
    client = llnl_db_client.LLNLDBClient(
            "/store/raw/LLNL/UCRL-MI-222502/westernus.wfdisc")

    # parameters for waveform request
    min_dist = 0 
    max_dist = 1200
    before = 100
    after = 600
    network = '' 
    channel = 'BH*'
    resample_freq = 0
    resample_freq = 20

    # parameters for CAP
    rotate = True
    output_cap_weight_file = True
    remove_response = True
    detrend = True
    demean = True
    output_event_info = True
    pre_filt = (0.005, 0.006, 5.0, 10.0)
    #scale_factor = 10.0**2  # original
    scale_factor = 2e-1     # Hoya. CHECK VALUES FOR ALL OTHERS

    evid = int(evid)
    try:
        getwaveform_llnl.run_get_waveform(llnl_db_client = client, event = evid, 
                network = network, channel = channel, 
                min_dist = min_dist, max_dist = max_dist, 
                before = before, after = after, 
                resample_freq = resample_freq,
                ifrotate = rotate, ifCapInp = output_cap_weight_file, 
                ifRemoveResponse = remove_response, ifDetrend = detrend, 
                ifDemean = demean, ifEvInfo = output_event_info, 
                scale_factor = scale_factor, pre_filt = pre_filt)
    except:
        print("WARNING. No waveforms returned with LLNL client")
        print("Continuing\n")
        pass

# List of events and their EVIDs from Ford (2009).
# First I matched event times in Ford to those in the LLNL database. 
# Then wrote their matching EVID.
# Using EVIDs I wrote a script to output event times from the LLNL database.
evids_times_explosions = {
        "584497":  "1988-02-15T18:10:00.000000Z", 
        "602554":  "1989-06-27T15:30:00.000000Z", 
        "605660":  "1989-09-14T15:00:00.100000Z", 
        "607806":  "1989-10-31T15:30:00.000000Z", 
        "609251":  "1989-12-08T15:00:00.000000Z",
        "612817":  "1990-03-10T16:00:00.000000Z", 
        "616762":  "1990-06-13T16:00:00.000000Z", 
        "617078":  "1990-06-21T18:15:00.000000Z", 
        "623055":  "1990-11-14T19:17:00.700000Z", 
        "627879":  "1991-03-08T21:02:45.000000Z",
        "628994":  "1991-04-04T19:00:00.000000Z", # no waveforms for NCEDC
        "635527":  "1991-09-14T19:00:00.000000Z", # Hoya
        "636899":  "1991-10-18T19:12:00.000000Z",
        "638595":  "1991-11-26T18:35:00.000000Z",
        "643767":  "1992-03-26T16:30:00.000000Z", 
        "653134":  "1992-09-18T17:00:00.000000Z", 
        "653332":  "1992-09-23T15:04:00.000000Z", 
        }

evids_times_quakes = {
        "648766": "1992-06-29T10:14:22.480000Z", 
        "649220": "1992-07-05T06:54:13.560000Z", 
        "706312": "1995-07-31T12:34:46.860000Z", 
        "737983": "1997-04-26T01:49:35.410000Z", 
        "743984": "1997-09-12T13:36:55.420000Z", 
        "768593": "1998-12-12T01:41:31.370000Z", 
        "770646": "1999-01-23T03:00:33.200000Z", 
        "770868": "1999-01-27T10:44:23.310000Z", 
        "1592802":"2002-06-14T12:40:44.450000Z" 
        }

evids_times_collapses = {
        "522227": "1982-08-05T14:00:00.000000Z",
        "697661": "1995-02-03T15:26:10.660000Z",
       "1324942": "2000-01-30T14:46:51.310000Z"
        }

event_list = [evids_times_explosions, evids_times_quakes,
        evids_times_collapses]

for evid, otime in evids_times_collapses.items():
    # get catalog data
    t0 = obspy.UTCDateTime(otime)
    cat0 = get_iris_evcat(t0)
    print(cat0)

    # download waveforms
    get_data_iris_ncedc(cat0)
    #get_data_llnl(evid)    # LLNL data. currently broken ...
