#=============================================================
# run_getwaveform.py
# 
# This script will fetch seismic waveforms, then process them, then write sac files.
# Used heavily within the UAF seismology group.
#
# This script contains a large number of examples in two categories:
# A. examples that target a current or previous bug
# B. examples of important events for modeling efforts that others may want to reproduce
#
# In the future, we will try to automatically run these examples for each code update.
# For now, we will simply try to regularly re-run the examples.
#
# contributors: Celso Alvizuri, Lion Krischer, Vipul Silwal, Carl Tape
# 
# To run this script:
# python run_getwaveform.py
#
# TO DO
# + filetags for the case iFilter = True (like lp10, bp10_40, etc)
# + provide better options and handling for data gaps (like: "toss waveform if gaps ar 0.01*length")
# + the KEVNM header cannot store the time to the 3rd millisecond character
#   probably the best approach is to write the spill-over character into another field
#   (or reconstruct the EID from the origin time, if that is your convention)
#
#=============================================================

import obspy
import copy
import util_helpers
import shutil   # only used for deleting data directory
import os
import sys
import getwaveform_new
import run_getwaveform_input

# EXAMPLES (choose one)
iex = 200
print("Running example iex =", iex)

#================================================================
# username and password for accessing embargoed data from IRIS
# Register here: http://ds.iris.edu/ds/nodes/dmc/forms/restricted-data-registration/
# Run example iex = 4 to check
user = ''
password = ''
#================================================================

# Get event info
ev_info = run_getwaveform_input.getwaveform_input()
ev_info.get_extraction_info(iex)

#================================================================
# fetch and process waveforms
# IRIS
if ev_info.idb == 1:
    # import functions to access waveforms
    #import getwaveform_iris
    from obspy.clients.fdsn import Client
    from obspy.core.event import Event, Origin, Magnitude
    if not user and not password:
        client = Client("IRIS")
    else:
        client = Client("IRIS",user=user,password=password)
    # will only work for events in the 'IRIS' catalog
    # (future: for Alaska events, read the AEC catalog)
    if ev_info.use_catalog==1:
        print("WARNING using event data from the IRIS catalog")
        cat = client.get_events(starttime = ev_info.otime - ev_info.sec_before_after_event,\
                                endtime = ev_info.otime + ev_info.sec_before_after_event)
        ev = cat[0]
        
        ref_time_place = ev
        if ev_info.rlat != ev_info.dummyval:
            ref_time_place.origins[0].latitude = ev_info.rlat
            ref_time_place.origins[0].longitude = ev_info.rlon
            ref_time_place.origins[0].time = ev_info.rtime 
    else:
        print("WARNING using event data from user-defined catalog")
        ev = Event()
        org = Origin()
        org.latitude = ev_info.elat
        org.longitude = ev_info.elon
        org.depth = ev_info.edep
        org.time = ev_info.otime
        mag = Magnitude()
        mag.mag = ev_info.emag
        mag.magnitude_type = "Mw"
        ev.origins.append(org)
        ev.magnitudes.append(mag)

        if ev_info.rlat == ev_info.dummyval:
            # By default this should be the event time and location unless we want to grab stations centered at another location
            ev_info.rlat = ev_info.elat
            ev_info.rlon = ev_info.elon
            ev_info.rtime = ev_info.otime
        
        ref_time_place = Event()
        ref_org = Origin()
        ref_org.latitude = ev_info.rlat
        ref_org.longitude = ev_info.rlon
        ref_org.time = ev_info.rtime
        ref_org.depth = 0 # dummy value
        ref_time_place.origins.append(ref_org)
        ref_time_place.magnitudes.append(mag) # more dummies

# LLNL
if ev_info.idb == 3:
    import llnl_db_client
    #import getwaveform_llnl
    client = llnl_db_client.LLNLDBClient(
            "/store/raw/LLNL/UCRL-MI-222502/westernus.wfdisc")

    # get event time and event ID
    # XXX this needs update
    cat = client.get_catalog()
    mintime_str = "time > %s" % (ev_info.otime - ev_info.sec_before_after_event)
    maxtime_str = "time < %s" % (ev_info.otime + ev_info.sec_before_after_event)
    print(mintime_str + "\n" + maxtime_str)
    #ev = cat.filter(mintime_str, maxtime_str)[0]
    ev = cat.filter(mintime_str, maxtime_str)
    
    if len(ev) > 0:
        ev = ev[0]
        # Nothing happens here.  We can change later
        ref_time_place = ev
        print(len(ev))
    else:
        print("No events in the catalog for the given time period. Stop.")
        sys.exit(0)

# Delete existing data directory
eid = util_helpers.otime2eid(ev.origins[0].time)
ddir = './'+ eid
#if os.path.exists('RAW'):
#    print("WARNING. %s already exists. Deleting ..." % ddir)
#    shutil.rmtree('RAW')
if ev_info.overwrite_ddir and os.path.exists(ddir):
    print("WARNING. %s already exists. Deleting ..." % ddir)
    shutil.rmtree(ddir)


# track git commit
os.system('git log | head -12 > ./' + eid + '_last_2git_commits.txt')

# Extract waveforms, IRIS
getwaveform_new.run_get_waveform(c = client, event = ev, ref_time_place = ref_time_place, ev_info = ev_info)
