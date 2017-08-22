#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

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

This script is based on script run_getwaveform_fmtu2016.py

20170816 cralvizuri <celso.alvizuri@gmail.com>
"""

import obspy
from getwaveform import *
def get_ev_info(ev_info, iex):
    #from getwaveform import *
    # import sys
    # 
    def client2ev_info(iev_info, iclient):
        """
        call run_get_waveform
        """

        # prepare waveform requests
        # IRIS and BK
        iev_info.idb = 1
        if "IRIS" in iclient:
            print("Working data client", iclient, "event", iev_info.evname)
            iev_info.network = '*,-XK,-XM,-XS,-XC,-XU,-XT,-XE'
            iev_info.station = '*,-PURD,-NV33,-GPO'
            iev_info.channel = 'BH?,LH?'
        elif "BK" in iclient:
            print("Working data client", iclient, "event", iev_info.evname)
            iev_info.network = 'BK'
            iev_info.station = '*'  # BK doesn't filter, use '*'
        elif "LLNL" in iclient:
            print("Working data client", iclient, "event", iev_info.evname)
            iev_info.idb = 3
            iev_info.station = '*'  # all stations
            iev_info.channel = '*'

        pass

    #events_file="/home/vipul/REPOSITORIES/manuscripts/alvizuri/papers/2016fmtu/data/event_info_llnl.txt"
    #events_file="/home/vipul/REPOSITORIES/GEOTOOLS/python_util/util_data_syn/event_info_llnl_test.txt"
    events_file="test_data/event_info_llnl.txt"
    events_file="test_data/event_info_llnl2.txt"
    # KERNVILLE              , 1988-02-15T18:10:00.09, -116.472,  37.314,   542, Ford2009,      5.30, ml, NCSN 
    # AMARILLO               , 1989-06-27T15:30:00.02, -116.354,  37.275,   640, Ford2009,      4.90, ml, NCSN 
    # DISKO_ELM              , 1989-09-14T15:00:00.10, -116.164,  37.236,   261, Ford2009,      4.40, ml, NCSN 

    fid = open(events_file, "r")
    #data = fid.readlines()[0:iex]
    data = fid.readlines()
    fid.close()

    ev_info = getwaveform()
    ev_info.overwrite_ddir = 1 
    ev_info.use_catalog = 0          # do not use event catalog for source parameters
    ev_info.min_dist = 0
    ev_info.max_dist = 1200
    ev_info.tbefore_sec = 100
    ev_info.tafter_sec = 600
    ev_info.scale_factor = 100 
    ev_info.resample_TF = False
    ev_info.resample_freq = 20        
    ev_info.f1 = 1/200  # fmin
    ev_info.f2 = 1/10  # fmax
    ev_info.ifsave_stationxml = False

    dblist = [1, 3]
    ev_info_list = []
    for row in data:
        # event objects for each event
        line = row.split()
        iev_info = ev_info.copy()
        iev_info.otime = obspy.UTCDateTime(line[1])
        iev_info.elon = line[2]
        iev_info.elat = line[3]
        iev_info.edep = line[4]
        iev_info.emag = line[6]
        #iev_info.evname = util_helpers.otime2eid(iev_info.otime)
        iev_info.evname = line[0]
        iev_info.get_events_client()

        client2ev_info(iev_info, 'IRIS'); ev_info_list.append(iev_info)
        client2ev_info(iev_info, 'LLNL'); ev_info_list.append(iev_info)
        client2ev_info(iev_info, 'BK'); ev_info_list.append(iev_info)

    return(ev_info_list)
