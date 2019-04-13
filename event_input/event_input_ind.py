import obspy
import read_event_obspy_file as reof
from getwaveform import *

def get_ev_info(ev_info,iex):
# ===============================================================
# Example Himalayan event
    if iex == 0:
        ev_info.use_catalog = 0
        ev_info.otime = obspy.UTCDateTime("2009-09-21T09:43:52.300")
        ev_info.min_dist = 0 
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300

        # List of  networks in Himalayas
        ev_info.network = 'Y2'
        #ev_info.network = 'ZO'
        ev_info.client_name = 'IRIS'
        ev_info.station = '*'
        ev_info.channel = 'H??,B??'
        ev_info.use_catalog = 0 
        ev_info.elat = 30.83
        ev_info.elon = 78.984
        ev_info.edep = 53000
        ev_info.emag = 4.7
        ev_info.ifverbose = True
        ev_info.resample_freq = 50
        ev_info.scale_factor = 100
        #ev_info.phase_window = False 
        #-------for specfem------------
        #ev_info.tbefore_sec = 0
        #ev_info.resample_TF = False
        #ev_info.scale_factor = 1
        #ev_info.outformat = 'DISP'
        #------------------------------
    return(ev_info)
