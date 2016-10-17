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
    sec_before_after = 2
    cat = c.get_events(starttime = t0-sec_before_after, endtime = t0+sec_before_after)
    if(len(cat) > 1):
        print("WARNING -- multiple events found. Using the first in the list.")
    return cat[0]

# get data  
def get_data_iris_ncedc(cat0):
    print("Calling getwaveform_iris")

    # parameters for waveform request
    min_dist = 0 
    max_dist = 1200
    before = 100
    after = 1200
    network = '*'
    station = '*'
    channel = 'BH*,LH*'
    resample_freq = 20      # 0 for no resampling

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
                    network = network, station = station, channel = channel,
                    resample_freq = resample_freq, 
                    ifrotate = rotate, ifCapInp = output_cap_weight_file, 
                    ifRemoveResponse = remove_response, ifDetrend = detrend, 
                    ifDemean = demean, ifEvInfo = output_event_info, 
                    scale_factor = scale_factor, pre_filt = pre_filt)
        except Exception as msg:
            print("WARNING. There was a problem with client %s" % thisclient)
            print(msg)
            print("Continuing ... \n")
            pass
    print("Done")

def get_data_llnl(evid):
    client = llnl_db_client.LLNLDBClient(
            "/store/raw/LLNL/UCRL-MI-222502/westernus.wfdisc")

    # parameters for waveform request
    min_dist = 0 
    max_dist = 1200
    before = 100
    after = 600
    network = '*' 
    station = '*' 
    channel = 'BH*,LH*,EH*'
    resample_freq = 20

    # parameters for CAP
    rotate = True
    output_cap_weight_file = True
    remove_response = True
    detrend = True
    demean = True
    output_event_info = True
    pre_filt = (0.005, 0.006, 5.0, 10.0)
    scale_factor = 10**2    # for CAP use 10**2 (to convert m/s to cm/s)

    #evid = int(evid)
    sec_before_after_event = 10  # time window to search for a target event in a catalog
    otime = obspy.UTCDateTime(evid)
    # get event time and event ID

    cat = client.get_catalog()
    mintime_str = "time > %s" % (otime - sec_before_after_event)
    maxtime_str = "time < %s" % (otime + sec_before_after_event)
    print(mintime_str, maxtime_str)
    cat0 = cat.filter(mintime_str, maxtime_str)[0]
    print(cat0)

    try:
        getwaveform_llnl.run_get_waveform(llnl_db_client = client, event = cat0, 
                network = network, station = station, channel = channel, 
                min_dist = min_dist, max_dist = max_dist, 
                before = before, after = after, 
                resample_freq = resample_freq,
                ifrotate = rotate, ifCapInp = output_cap_weight_file, 
                ifRemoveResponse = remove_response, ifDetrend = detrend, 
                ifDemean = demean, ifEvInfo = output_event_info, 
                scale_factor = scale_factor, pre_filt = pre_filt)
    except Exception as msg:
        print("\n==========================================================\n")
        print("WARNING. Event ", cat0)
        print("There was a problem with LLNL client:")
        print(msg)
        print("Continuing to next event...")
        print("\n==========================================================\n")
        pass
    print("Done")

def getdata_iris_llnl(dataset, llnl=False, iris=False):
    """
    get data from LLNL and IRIS + NCEDC
    """
    for evid, otime in dataset.items():
        # get catalog data
        t0 = obspy.UTCDateTime(otime)
        cat0 = get_iris_evcat(t0)
        print(cat0)

        # download waveforms
        if(llnl): get_data_llnl(otime)
        if(iris): get_data_iris_ncedc(cat0)

# List of events and their EVIDs from Ford (2009).
# First I matched event times in Ford to those in the LLNL database. 
# Then wrote their matching EVID.
# Using EVIDs I wrote a script to output event times from the LLNL database.
evids_times_explosions = {
        "584497":  "1988-02-15T18:10:00.000000Z",
        "602554":  "1989-06-27T15:30:00.000000Z",
        "605660":  "1989-09-14T15:00:00.100000Z",  # KNB, MNV. problem applying rotate2zne
        "607806":  "1989-10-31T15:30:00.000000Z",
        "609251":  "1989-12-08T15:00:00.000000Z",
        "612817":  "1990-03-10T16:00:00.000000Z",
        "616762":  "1990-06-13T16:00:00.000000Z",
        "617078":  "1990-06-21T18:15:00.000000Z",
        "623055":  "1990-11-14T19:17:00.700000Z", 
        "627879":  "1991-03-08T21:02:45.000000Z", 
        "628994":  "1991-04-04T19:00:00.000000Z", # no waveforms from NCEDC.
        "635527":  "1991-09-14T19:00:00.000000Z", # ok (Hoya)
        "636899":  "1991-10-18T19:12:00.000000Z",
        "638595":  "1991-11-26T18:35:00.000000Z", 
        "643767":  "1992-03-26T16:30:00.000000Z", 
        "653134":  "1992-09-18T17:00:00.000000Z",  # KNB, LAC, MNV are very short traces (useless)
        "653332":  "1992-09-23T15:04:00.000000Z"  
        }

evids_times_quakes = {
        "648766": "1992-06-29T10:14:22.480000Z", 
        "649220": "1992-07-05T06:54:13.560000Z", 
        "706312": "1995-07-31T12:34:46.860000Z", # No waveform managed to get instrument corrected (LLNL)
        "737983": "1997-04-26T01:49:35.410000Z",
        "743984": "1997-09-12T13:36:55.420000Z", 
        "768593": "1998-12-12T01:41:31.370000Z", 
        "770646": "1999-01-23T03:00:33.200000Z",
        "770868": "1999-01-27T10:44:23.310000Z", 
        "1592802":"2002-06-14T12:40:44.450000Z",  
        #"-99999":  "1996-09-05T08:16:56.090000Z",# IRIS returns deeper event (7.5km, not 5km) 1996-09-05T08:16:55.780000Z | +36.745, -116.282 | 3.2 mL
        #"-99999": "1997-06-14T19:48:19.930000Z", # IRIS returns shallower event (5km, not 7km) 1997-06-14T19:48:21.200000Z | +36.806, -115.847 | 3.4 mL
        #"-99999": "2007-01-24T11:30:16.100000Z"  # IRIS returns smaller event (not 4.09) 2007-01-24T11:30:16.650000Z | +37.424, -117.064 | 2.9 mb
        }

evids_times_collapses = {
        "522227": "1982-08-05T14:00:00.000000Z", # No waveform managed to get instrument corrected (LLNL)
        "697661": "1995-02-03T15:26:10.660000Z",
        "1324942": "2000-01-30T14:46:51.310000Z" # No waveform managed to get instrument corrected (LLNL)
        }

# get the waveforms
getdata_iris_llnl(evids_times_explosions, llnl=True, iris=True)
getdata_iris_llnl(evids_times_quakes, llnl=True, iris=True)
getdata_iris_llnl(evids_times_collapses, llnl=True, iris=True)

