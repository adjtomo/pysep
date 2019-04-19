import obspy
import read_event_obspy_file as reof
from getwaveform import *

def get_ev_info(ev_info,iex):
# ===============================================================
# Example Himalayan event
    if iex == 0:
        ev_info.use_catalog = 0
        ev_info.otime = obspy.UTCDateTime("2009-09-21T09:43:52.300")
        ev_info.otime = obspy.UTCDateTime("2010-06-22T23:14:11.930")
        ev_info.otime = obspy.UTCDateTime("2010-07-06T19:08:23.700")
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
        #----------------------
        ev_info.elat = 30.83
        ev_info.elon = 78.984
        ev_info.edep = 53000
        ev_info.emag = 4.7
        #----------------------
        ev_info.elat = 29.88
        ev_info.elon = 80.47
        ev_info.edep = 25000
        ev_info.emag = 5.3
        #----------------------
        ev_info.elat = 29.55
        ev_info.elon = 80.33
        ev_info.edep = 14000
        ev_info.emag = 4.9
        #----------------------
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
        
    if iex == 1:
        ev_info.use_catalog = 0
        ev_info.otime = obspy.UTCDateTime("2015-12-18T22:16:56.460")
        ev_info.min_dist = 0 
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300

        # List of  networks in Himalayas
        ev_info.network = 'Y2'
        ev_info.network = 'ZO'
        ev_info.client_name = 'ETH'
        ev_info.station = '*'
        ev_info.channel = 'H??,B??'
        ev_info.use_catalog = 0 
        ev_info.elat = 29.33
        ev_info.elon =81.64
        ev_info.edep = 16000
        ev_info.emag = 5.25
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
