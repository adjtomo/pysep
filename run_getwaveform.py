"""
Run script for getwaveform.py

This script will fetch seismic waveforms, process them, then write sac files.
Used heavily within the UAF seismology group.

.. rubric::
    $ python run_getwaveform.py NAME_OF_EVENT_INPUT_FILE EXAMPLE_INDEX
    $ python run_getwaveform.py event_input 3

.. note::
    event input file needs to be linked locally:
    $ ln -s event_input/event_input_flats.py .

Option B: within a bash script
check_getwaveform.bash will run all examples


In the future, we will try to automatically run these examples for
each code update. For now, we will simply try to regularly re-run the examples.

.. note:: Examples
    event_input               -- default examples (including CAP examples)
    event_input_testing       -- examples used to isolate bugs or special cases
    event_input_scak          -- southern Alaska
    event_input_flats         -- Minto Flats (iex=3; multiple events)
    gw_fmtu                   -- llnl paper. NOTE iex doesnt work here
    event_input_mtuq          -- Test events for MTUQ project

.. note:: To Do
    + filetags for the case iFilter = True (like lp10, bp10_40, etc)
    + provide better options and handling for data gaps (like: "toss
      waveform if gapar 0.01*length")
    + the KEVNM header cannot store the time to the 3rd millisecond
      character probably the best approach is to write the spill-over
      character into another field (or reconstruct the EID from the
      origin time, if that is your convention)
"""
import os
import shutil
import warnings
import argparse
from getwaveform import GetWaveform
import util_helpers


def getargs():
    """
    Get command line arguments
    :return:
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("iproject", nargs="?",
                        help="name of file containing event info",
                        default="event_input")

    parser.add_argument("iex", type=int, nargs="?",
                        help="example number within iproject.py script",
                        default="0")

    return parser.parse_args()


if __name__ == '__main__':
    args = getargs()
    print(f"Running example iex = {args.iex}")

    ev_info = GetWaveform()
    # iproject defines the event_input information cointained in a .py file
    iproject = __import__(args.iproject)

    # Update default extraction parameters with inputted parameters
    ev_info = iproject.get_ev_info(ev_info, args.iex)

    # create list if ev_info is a getwaveform object
    if type(ev_info) not in [list, tuple]:
        ev_info_list = [ev_info]
    else:
        ev_info_list = ev_info

    # Fetch and process waveforms
    for ii in range(len(ev_info_list)):
        ev_info = ev_info_list[ii]

        # Get event and client info, create event objects and reference
        # origin object (same type as event object) for station selection
        ev_info.get_events_client()

        # Delete existing data directory
        # set ev_info.overwrite_ddir=false to extract data from multiple clients
        ev_info.evname = util_helpers.otime2eid(
            ev_info.ref_time_place.origins[0].time)
        ddir = os.path.join(os.getcwd(), ev_info.evname)
        if ev_info.overwrite_ddir and os.path.exists(ddir):
            shutil.rmtree(ddir)
            warnings.warn("\nDeleted existing directory:\n{ddir}\n")


        # Extract waveforms from IRIS and save output information
        ev_info.run_get_waveform()
        ev_info.save_extraction_info()
        
