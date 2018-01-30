import obspy
from getwaveform import *

def get_ev_info(ev_info,iex):
    if iex is 0:
        # Virginia earthquake
        # https://earthquake.usgs.gov/earthquakes/events/2011virginia/overview.php
        print('Extracting data for SALMON EQ')
        ev_info.idb = 1   # 
        ev_info.otime = obspy.UTCDateTime("2017-05-16T04:20:48.991")
        ev_info.min_lat = 58.0 
        ev_info.max_lat = 61.5
        ev_info.min_lon = -156.0
        ev_info.max_lon = -150.0
        ev_info.tbefore_sec = 20
        ev_info.tafter_sec = 120
        ev_info.network = 'AK,AT,AV,CN,II,IU,US,XM,TA,XE,XR,XZ,YV,XV,ZE,XG'
        ev_info.channel = 'BH?,HH?'
        ev_info.remove_response = False
        ev_info.resample_TF = False
        #ev_info.use_catalog = 1 
        # Shorter run
        ev_info.phase_window = True
        ev_info.phases = ["P","S"] # Choosing time period with respect to P & S
        #ev_info.phases = ["P","P"]
        ev_info.write_sac_phase = True
        ev_info.ipre_filt = 1 # 0 = no prefiter.  1 = default prefilter (see getwaveform_iris.py) where is this file?
        ev_info.output_cap_weight_file = False
        ev_info.ifsave_stationxml = False
            # For SALMON data
        ev_info.user = "csrichards2@alaska.edu"
        ev_info.password = "Pimd2jg8hVbD"
           # ev_info_temp.user = None
           # ev_info_temp.password = None

    return ev_info
