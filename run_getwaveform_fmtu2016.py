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

import os
import shutil
import util_helpers
import getwaveform as gw
import llnl_db_client
import obspy
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.core.event import Event, Origin, Magnitude

def getwf_iris_ncedc_llnl(origin0, client_pick):

    print("Running function getwf_iris_ncedc_llnl")

    # parameters for waveform request
    tbefore_sec = 100
    tafter_sec = 600

    # DEFAULT SETTINGS (see getwaveform_iris.py)
    rotateRTZ = True
    rotateUVW = False   # works only if 'rotateRTZ = True'
    output_cap_weight_file = True
    detrend = True
    demean = True
    output_event_info = True
    taper = False
    ifplot_spectrogram = False

    # for CAP all waveforms need to have the same sample rate
    resample_freq = 20.0         # 0 for no resampling
    scale_factor = 10**2         # for CAP use 10**2  (to convert m/s to cm/s)

    # event parameters
    sec_before_after_event = 10  # time window to search for a target event in a catalog
    min_dist = 0 
    max_dist = 1200

    # station parameters
    network = '*'                # all networks
    station = '*,-PURD,-NV33,-GPO'  # all stations
    channel = 'BH?,LH?'

    overwrite_ddir = 1           # 1 = delete data directory if it already exists
    icreateNull = 0              # create Null traces so that rotation can work (obsby stream.rotate require 3 traces)

    # filter
    # set ipre_filt = 0 to prevent extra filtering
    ifFilter = False 
    filter_type = 'bandpass'
    # LLNL filter 10-50 sec is most stable. Some waveforms >= 100 sec show
    # oscillations. Use factor of 2 limit = 200 sec for prefilter
    f1 = 1/200
    f2 = 1/10
    zerophase = True             # False = causal, True = acausal
    corners = 4                  # Is corner in Obspy same as Pole in SAC?
    remove_response = True
    iplot_response = False
    ipre_filt = 2                # 0 No pre_filter
                                 # 1 default pre_filter (see getwaveform_iris.py)
                                 # 2 user-defined pre_filter
    f0 = 0.5 * f1
    f3 = 2.0 * f2
    # The following are for the FMTU paper
    f0 = 0.005; f1 = 0.006; f3 = 10; f4 = 15
    pre_filt = (f0, f1, f2, f3)  # applies for ipre_filt = 2 only

    # NOTE event data from user-defined catalog!
    # initialize objects
    ev = Event()
    org = Origin()
    mag = Magnitude()

    # build objects
    org.time        = UTCDateTime(origin0[0])
    org.longitude   = origin0[1]
    org.latitude    = origin0[2]
    org.depth       = origin0[3]
    mag.mag         = origin0[5]
    mag.magnitude_type = origin0[6]    # Mw, ml, mb, ...

    ev.origins.append(org)
    ev.magnitudes.append(mag)

    # Delete existing data directory
    eid = util_helpers.otime2eid(ev.origins[0].time)
    ddir = './'+ eid
    if os.path.exists('RAW'):
        print("WARNING. %s already exists. Deleting ..." % ddir)
        shutil.rmtree('RAW')
    if overwrite_ddir and os.path.exists(ddir):
        print("WARNING. %s already exists. Deleting ..." % ddir)
        shutil.rmtree(ddir)

    if client_pick is "IRIS":
        print('Using client %s' % client_pick)
        idb = 1 
        client_list = ["IRIS", "NCEDC"]
        print("WARNING using event data from user-defined catalog")

    # LLNL
    if client_pick is "LLNL":
        print('Using client %s' % client_pick)
        idb = 3
        client_list = ["LLNL"]

        client = llnl_db_client.LLNLDBClient(
                "/store/raw/LLNL/UCRL-MI-222502/westernus.wfdisc")
        # get event time and event ID
        otime = obspy.UTCDateTime(origin0[0])
        cat = client.get_catalog()
        mintime_str = "time > %s" % (otime - sec_before_after_event)
        maxtime_str = "time < %s" % (otime + sec_before_after_event)
        print(mintime_str + "\n" + maxtime_str)
        ev = cat.filter(mintime_str, maxtime_str)
        nev = len(ev)
        if nev == 1:
            ev = ev[0]  # doesn't have magnitude (ATRISCO)
        elif nev > 1:
            ev = ev[1]  # [0] may not include magnitude. [1] may (LLNL)
        else:
            print("No events in the catalog for the given time period. Stop.")

    # The IRIS requests include BK data, but it has to be requested through 
    # the NCEDC client
    for iclient in client_list:
        if iclient is "IRIS":
            network = network
            client = Client(iclient)
        elif iclient is "NCEDC":
            network = 'BK'
            station = '*'   # doesn't like "-staX"
            client = Client(iclient)

        gw.run_get_waveform(c = client, event = ev, idb = idb, 
                min_dist = min_dist, max_dist = max_dist, 
                before = tbefore_sec, after = tafter_sec, 
                network = network, station = station, channel = channel, 
                resample_freq = resample_freq,
                ifrotateRTZ = rotateRTZ, ifrotateUVW = rotateUVW,
                ifCapInp = output_cap_weight_file, 
                ifRemoveResponse = remove_response,
                ifDetrend = detrend, ifDemean = demean, ifTaper = taper,
                ifEvInfo = output_event_info,
                scale_factor = scale_factor,
                icreateNull = icreateNull,
                ipre_filt = ipre_filt, pre_filt = pre_filt, ifFilter = ifFilter, 
                fmin = f1, fmax = f2, filter_type = filter_type, 
                zerophase = zerophase, corners = corners, 
                iplot_response = iplot_response, 
                ifplot_spectrogram = ifplot_spectrogram)

def call_getwf(dataset, llnl=False, iris=False):
    """
    get data from LLNL and IRIS + NCEDC
    """
    for evname, origin in dataset.items():
        print("Processing event %s %s" % (evname, origin))
        # download waveforms
        if(iris):
            getwf_iris_ncedc_llnl(origin, client_pick="IRIS")
        if(llnl):
            getwf_iris_ncedc_llnl(origin, client_pick="LLNL")

# List of events and their EVIDs from Ford (2009).
# First I matched event times in Ford to those in the LLNL database. 
# Then wrote their matching EVID.
# Using EVIDs I wrote a script to output event times from the LLNL database.
origins_explosions = {
        "KERNVILLE     ": ["1988-02-15T18:10:00.09", "-116.472", "37.314", "542", "Ford2009", "5.30", "ml", "NCSN"], 
        "AMARILLO      ": ["1989-06-27T15:30:00.02", "-116.354", "37.275", "640", "Ford2009", "4.90", "ml", "NCSN"], 
        "DISKO_ELM     ": ["1989-09-14T15:00:00.10", "-116.164", "37.236", "261", "Ford2009", "4.40", "ml", "NCSN"], 
        "HORNITOS      ": ["1989-10-31T15:30:00.09", "-116.492", "37.263", "564", "Ford2009", "5.40", "ml", "NCSN"], 
        "BAMWELL       ": ["1989-12-08T15:00:00.09", "-116.410", "37.231", "601", "Ford2009", "5.30", "ml", "NCSN"], 
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
call_getwf(origins_explosions, llnl=True, iris=False)
call_getwf(origins_quakes,     llnl=False, iris=False)
call_getwf(origins_collapses,  llnl=False, iris=False)

