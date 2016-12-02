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
    #station = '*'
    station = '*,-PURD,-NV33,-GPO'
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
origins_explosions = {
        "KERNVILLE     ": ["1988-02-15T18:10:00.09", "-116.472", "37.314", "542", "Ford2009", "5.30", "ml", "NCSN"], 
        "AMARILLO      ": ["1989-06-27T15:30:00.02", "-116.354", "37.275", "640", "Ford2009", "4.90", "ml", "NCSN"], 
        "DISKO_ELM     ": ["1989-09-14T15:00:00.10", "-116.164", "37.236", "261", "Ford2009", "4.40", "ml", "NCSN"], 
        "HORNITOS      ": ["1989-10-31T15:30:00.09", "-116.492", "37.263", "564", "Ford2009", "5.40", "ml", "NCSN"], 
        "BARNWELL      ": ["1989-12-08T15:00:00.09", "-116.410", "37.231", "601", "Ford2009", "5.30", "ml", "NCSN"], 
        "METROPOLIS    ": ["1990-03-10T16:00:00.08", "-116.056", "37.112", "469", "Ford2009", "4.94", "md", "NCSN"], 
        "BULLION       ": ["1990-06-13T16:00:00.09", "-116.421", "37.262", "674", "Ford2009", "5.34", "md", "NCSN"], 
        "AUSTIN        ": ["1990-06-21T18:15:00.00", "-116.005", "36.993", "350", "Ford2009", "4.11", "md", "NCSN"], 
        "HOUSTON       ": ["1990-11-14T19:17:00.07", "-116.372", "37.227", "594", "Ford2009", "4.86", "md", "NCSN"], 
        "COSO          ": ["1991-03-08T21:02:45.08", "-116.075", "37.104", "417", "Ford2009", "4.50", "ml", "NCSN"], 
        "BEXAR         ": ["1991-04-04T19:00:00.00", "-116.314", "37.296", "629", "Ford2009", "5.08", "md", "NCSN"], 
        "HOYA          ": ["1991-09-14T19:00:00.08", "-116.429", "37.226", "658", "Ford2009", "5.40", "ml", "NCSN"], 
        "LUBBOCK       ": ["1991-10-18T19:12:00.00", "-116.046", "37.063", "457", "Ford2009", "4.75", "md", "NCSN"], 
        "BRISTOL       ": ["1991-11-26T18:35:00.07", "-116.070", "37.096", "457", "Ford2009", "4.80", "ml", "NCSN"], 
        "JUNCTION      ": ["1992-03-26T16:30:00.00", "-116.361", "37.272", "622", "Ford2009", "4.82", "ml", "NCSN"], 
        "HUNTERS_TROPHY": ["1992-09-18T17:00:00.08", "-116.211", "37.207", "385", "Ford2009", "3.87", "md", "NCSN"], 
        "DIVIDER       ": ["1992-09-23T15:04:00.00", "-115.989", "37.021", "340", "Ford2009", "4.13", "md", "NCSN"] 
        }

origins_quakes = {
"      Little_Skull_Main": ["1992-06-29T10:14:23.18 ", "-116.289 ", "36.698  ", "9070.0", " CI", "5.4 ", "ms ", "CI "], 
"Little_Skull_Aftershock": ["1992-07-05T06:54:13.52 ", "-116.276 ", "36.685  ", "5070.0", " CI", "4.44", "ml ", "CI "], 
"        Timber_Mountain": ["1995-07-31T12:34:47.35 ", "-116.415 ", "37.107  ", "4487.0", " CI", "4.0 ", "ml ", "CI "], 
"               Amargosa": ["1996-09-05T08:16:55.40 ", "-116.266 ", "36.729  ", "9070.0", " CI", "3.7 ", "ml ", "CI "], 
"             Groom_Pass": ["1997-04-26T01:49:35.20 ", "-115.937 ", "37.16   ", "5000.0", " US", "4.3 ", "ml ", "US "], 
"         Indian_Springs": ["1997-06-14T19:48:19.45 ", "-115.801 ", "36.648  ", "5793.0", " CI", "3.81", "ml ", "CI "], 
"             Calico_Fan": ["1997-09-12T13:36:54.94 ", "-116.277 ", "36.887  ", "6037.0", " CI", "4.05", "ml ", "CI "], 
"           Warm_Springs": ["1998-12-12T01:41:32.00 ", "-116.29  ", "37.51   ", "0.0   ", "REN", "4.1 ", "mb ", "US "], 
"       Frenchman_Flat_1": ["1999-01-23T03:00:32.00 ", "-115.92  ", "36.82   ", "0.0   ", "REN", "3.7 ", "ml ", "US "], 
"       Frenchman_Flat_2": ["1999-01-27T10:44:23.30 ", "-115.989 ", "36.816  ", "5000.0", "US ", "4.8 ", "mwr", "BRK"], 
"           Little_Skull": ["2002-06-14T12:40:45.36 ", "-116.293 ", "36.7163333", "9653.0", "CI ", "4.58", "mw ", "CI "], 
"                Ralston": ["2007-01-24T11:30:16.099", "-117.0986", "37.4133 ", "6100.0", "NN ", "4.1 ", "ml ", "NN "] 
}

origins_collapses = {
        "ATRISCO_Hole": ["1982-08-05T14:21:38.000", "-116.0065", "37.0842", " 320.0", "Walter2009", "3.5", "ms", "Walter2009"], 
        "Trona_Mine_1": ["1995-02-03T15:26:10.690", "-109.640 ", "41.529 ", "1000.0", "  US      ", "5.3", "mb", "US"], 
        "Trona_Mine_2": ["2000-01-30T14:46:51.310", "-109.679 ", "41.464 ", "1000.0", "  US      ", "4.4", "mb", "US"] 
        }

# get the waveforms
getdata_iris_llnl(origins_explosions, llnl=True, iris=True)
getdata_iris_llnl(origins_quakes, llnl=True, iris=True)
getdata_iris_llnl(origins_collapses, llnl=True, iris=True)

