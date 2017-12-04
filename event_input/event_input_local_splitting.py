import obspy
import read_event_obspy_file as reof
from getwaveform import *

def get_ev_info(ev_info,iex):

# Loop over multiple event       
    if iex == 1:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        
        #events_file = "/home/ksmith/REPOSITORIES/manuscripts/kyle/papers/basinamp/data/basinamp_obspy.txt"
        events_file = "/home/crichards/LOCAL_SPLITTING/salmon_redoubt_obspy.txt" #output from rs_redoubt.m
        eids,otimes,elons,elats,edeps,emags = reof.read_events_obspy_file(events_file)

        ev_info_list = []
        for xx in range(7,8):
            ev_info_temp = ev_info.copy()
            ev_info_temp.otime = obspy.UTCDateTime(otimes[xx])
            ev_info_temp.elat = elats[xx]
            ev_info_temp.elon = elons[xx]
            ev_info_temp.edep = edeps[xx]
            ev_info_temp.emag = emags[xx]
            ev_info_temp.eid = eids[xx]
            
            # subset of stations -- use lon-lat subset region
            #ev_info_temp.min_dist = 0
            #ev_info_temp.max_dist = 200
            ev_info_temp.min_lat = 59.0
            ev_info_temp.max_lat = 63.0
            ev_info_temp.min_lon = -154.0
            ev_info_temp.max_lon = -148.0
            #ev_info_temp.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
            ev_info_temp.network = 'AK,AT,AV,CN,II,IU,US,XM,TA,XE,XR,XZ,YV,XV,ZE,XG'
            ev_info_temp.channel = 'BH?,HH?'

            ev_info_temp.tbefore_sec = 20
            ev_info_temp.tafter_sec = 120
            #ev_info_temp.resample_freq = 40
            #ev_info_temp.scale_factor = 100
            ev_info_temp.resample_TF = False

            #ev_info_temp.rotateRTZ = True
            ev_info_temp.remove_response = False

            # For SALMON data
            ev_info_temp.user = None
            ev_info_temp.password = None

            # append getwaveform objects
            ev_info_list.append(ev_info_temp)
        
        # always return ev_info
        ev_info = ev_info_list

    
    return(ev_info)
#=================================================================================
# END EXAMPLES
#=================================================================================
