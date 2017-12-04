'''
# ============================================================================
run_getwaveform.py
This script will fetch seismic waveforms, then process them, then write sac files.
Used heavily within the UAF seismology group.

This script contains a large number of examples in two categories:
A. examples that target a current or previous bug
B. examples of important events for modeling efforts that others
   may want to reproduce

In the future, we will try to automatically run these examples for
each code update.
For now, we will simply try to regularly re-run the examples.
contributors: Celso Alvizuri, Lion Krischer, Vipul Silwal, Carl Tape

To run this script:
Option A: from the command line
python run_getwaveform.py NAME_OF_EVENT_INPUT_FILE EXAMPLE_INDEX
python run_getwaveform.py event_input 3
event input file needs to be linked locally:
   ln -s event_input/event_input_flats.py .

Option B: within a bash script
check_getwaveform.bash will run all examples

TO DO
+ filetags for the case iFilter = True (like lp10, bp10_40, etc)
+ provide better options and handling for data gaps (like: "toss
  waveform if gapar 0.01*length")
+ the KEVNM header cannot store the time to the 3rd millisecond
  character probably the best approach is to write the spill-over
  character into another field (or reconstruct the EID from the
  origin time, if that is your convention)
# =============================================================================
'''

import util_helpers
import shutil
import os
import sys
import argparse
from getwaveform import getwaveform

'''
EXAMPLES (choose one)
event_input               -- default examples (including CAP examples)
event_input_testing       -- examples used to isolate bugs or special cases
event_input_scak          -- southern Alaska
event_input_flats         -- Minto Flats (iex=3; multiple events)
gw_fmtu                   -- llnl paper. NOTE iex doesnt work here
event_input_mtuq          -- Test events for MTUQ project
'''



def getargs():
    parser = argparse.ArgumentParser()

    parser.add_argument('iproject', nargs='?',
        help='name of file containing event info',
        default='event_input')

    parser.add_argument('iex', type=int, nargs='?',
        help='example number within iproject.py script',
        default='0')

    return parser.parse_args()



if __name__ == '__main__':
    args = getargs()

    print("Running example iex =", args.iex)

    # =============================================================================
    # Get extraction info
    # see getwaveform_new.py for input parameters dedcription
    ev_info = getwaveform()              # create getwaveforms object
    iproject = __import__(args.iproject)      # file containing all the examples

    # update default extraction parameters with inputted parameters
    ev_info = iproject.get_ev_info(ev_info, args.iex)

    # For looping over events
    # create list if ev_info is a getwaveform object
    if type(ev_info) not in [list, tuple]: 
        ev_info_list = [ev_info]
    else:
        ev_info_list = ev_info

    # =============================================================================
    # fetch and process waveforms

    nev = len(ev_info_list)

    for ii in range(nev):
        ev_info = ev_info_list[ii]

        # Get event and client info
        # create event objects and reference origin object
        # (same type as event object) for station selection
        ev_info.get_events_client()

        # Delete existing data directory
        # set ev_info.overwrite_ddir=false to extract data from multiple clients
        ev_info.evname = util_helpers.otime2eid(
            ev_info.ref_time_place.origins[0].time)
        ddir = './' + ev_info.evname
        if ev_info.overwrite_ddir and os.path.exists(ddir):
            print("\n*** WARNING *** . Deleted existing directory", ddir)
            shutil.rmtree(ddir)

        # KEY: Extract waveforms, IRIS
        # use try to recover nicely when there are multiple event requests
        ev_info.run_get_waveform()
        '''
        try:
            ev_info.run_get_waveform()
        except Exception as e:
            print("FATAL. Unable to get data for event", ev_info.evname)
            print(e)
            print("Continuing with next event\n")
            continue
        # save extraction info (git commit and filenames)
        ev_info.save_extraction_info()
    '''
