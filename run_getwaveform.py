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
#================================================================

import obspy
import copy
import util_helpers
import shutil   # only used for deleting data directory
import os
import sys
from getwaveform import *

# EXAMPLES (choose one)
# event_input             -- 10 default examples (including CAP examples)
# event_input_testing     -- examples used to isolate bugs or special cases
# event_input_scak        -- southern Alaska
# event_input_flats       -- Minto Flats
iproject = 'event_input'   # this is the name of file containing event info (See run_getwaveform_input.py for example)
iex = 0                    # example number within iproject.py script
                           # iproject = 'event_input_flats', iex = 6  (for looping over multiple events)

#iproject = 'gw_fmtu'    # llnl paper. NOTE iex doesnt work here

# Or parse command line input arguments
if len(sys.argv)==3:
    iproject = str(sys.argv[1])
    iex = int(sys.argv[2])

print("Running example iex =", iex)

#================================================================
# Get extraction info
# see getwaveform_new.py for input parameters dedcription
ev_info = getwaveform()              # create event object (contains default extraction parameters)
iproject = __import__(iproject)      # file containing all the examples
ev_info = iproject.get_ev_info(ev_info, iex)   # update default extraction parameters with inputted event extarction parameters
# For looping over events
# create list if ev_info is a getwaveform object
if type(ev_info ) == getwaveform:
    ev_info_list = [ev_info]
else:
    ev_info_list  = ev_info

#================================================================
# fetch and process waveforms
# IRIS

nev = len(ev_info_list)
ndir=0
dir_exists=False
for ii in range(nev):
    ev_info = ev_info_list[ii]
    
    # Get event and client info
    # create event objects and reference origin object (same type as event object) for station selection
    ev_info.get_events_client()

    # Delete existing data directory
    ev_info.evname = util_helpers.otime2eid(ev_info.ref_time_place.origins[0].time)
    ddir = './'+ ev_info.evname
    if ev_info.overwrite_ddir and os.path.exists(ddir):
        # this command does not allow to keep data from multiple sources.
        # replacing it with Warning message
        #shutil.rmtree(ddir)

        print("\n*** WARNING *** %s. This directory already exists.", ddir)
        print("Consider deleting it and rerunning this script\n")
        dir_exists=True
        ndir+=1


    # KEY COMMAND: Extract waveforms, IRIS
    # use try to recover nicely when there are multiple event requests
    try:
        ev_info.run_get_waveform()
    except Exception as e:
        print("FATAL. Unable to get data for event", ev_info.evname)
        print(e)
        print("Continuing with next event\n")
        continue

    # save extraction info (git commit and filenames)
    ev_info.save_extraction_info()

if dir_exists:
    print("WARNING. There were %d pre-existing directories." % ndir)

