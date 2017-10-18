import obspy
import read_event_obspy_file as reof
from getwaveform import *

def get_ev_info(ev_info,iex):
# ===============================================================

# tomoDD events for Melissa (and local splitting?)     
    if iex == 0:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        ev_info.resample_TF = False
        ev_info.scale_factor = 1         # no scale factor
        ev_info.output_cap_weight_file = False
        ev_info.ifsave_stationxml = False

        events_file = "/home/carltape/PROJECTS/SALMON/data/salmon_tomodd_obspy_sub.txt"
        eids,otimes,elons,elats,edeps,emags = reof.read_events_obspy_file(events_file)
        
        ev_info.ifFilter = True
        ev_info.filter_type = 'highpass'
        ev_info.f1 = 1.0

        # subset of stations
        #ev_info.min_dist = 0
        #ev_info.max_dist = 200
        ev_info.min_lat = 59
        ev_info.max_lat = 62
        ev_info.min_lon = -156
        ev_info.max_lon = -148

        ev_info.tbefore_sec = 20
        ev_info.tafter_sec = 140
        ev_info.network = 'AK,AT,AV,CN,II,IU,US,XM,TA,XE,XR,XZ,YV,XV,ZE,XG'
        #ev_info.network = "AK"
        ev_info.channel = 'BH?,HH?'

        ev_info_list = []
        #for xx in range(5,10):
        for xx in range(len(eids)):
            ev_info_temp = ev_info.copy()
            ev_info_temp.otime = obspy.UTCDateTime(otimes[xx])
            ev_info_temp.elat = elats[xx]
            ev_info_temp.elon = elons[xx]
            ev_info_temp.edep = edeps[xx]
            ev_info_temp.emag = emags[xx]
            ev_info_temp.eid = eids[xx]
            # append getwaveform objects
            ev_info_list.append(ev_info_temp)
        
        # always return ev_info
        ev_info = ev_info_list
