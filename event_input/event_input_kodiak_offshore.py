import obspy
import read_event_obspy_file as reof
from getwaveform import *

def get_ev_info(ev_info,iex):
# ===============================================================
# Mw 8 Kodiak, Jan 23 2018
    if iex == 0:
        ev_info.use_catalog = 1
        ev_info.otime = obspy.UTCDateTime("2018-01-23 09:31:42")
        ev_info.min_dist = 0 
        ev_info.max_dist = 1000
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = 'AK,CN,II,IU,US,XM,XV,XZ,YV,TA'  # note: cannot use '*' because of IM
        ev_info.channel = 'BH?,HH?'
        #ev_info.channel = 'BN?,HN?,EN?'    # strong motion
        ev_info.use_catalog = 0 
        ev_info.elat = 56.046
        ev_info.elon = -149.073
        ev_info.edep = 25000
        ev_info.emag = 7.9
        ev_info.resample_freq = 50         # for CAP
        ev_info.scale_factor = 100         # for CAP

        # parameters for examining step response (causal low-pass filter on raw waveforms)
        # delete AUQ, SPCP, SPBG
        #ev_info.ifFilter = True
        #ev_info.filter_type = 'lowpass'
        #ev_info.f1 = 1/4
        #ev_info.zerophase = False
        #ev_info.remove_response = False
        #ev_info.ipre_filt = 0
        #ev_info.demean = False
        #ev_info.detrend = False
        #ev_info.resample_TF = False        # if False then resample_freq is taken from SAC files
        #ev_info.scale_factor = 1           # no scale factor
        #ev_info.rotateRTZ = False
        #ev_info.rotateUVW = True
        #ev_info.isave_raw = True

#-----------------------------------------------------------------
# Kodiak Offshore earthquake aftershocks
# XXX Replace it with a file (Extract multiple events)
    if iex == 1:
        ev_info.use_catalog = 1
        ev_info.otime = obspy.UTCDateTime("2018-01-26 01:09:30")
        ev_info.min_dist = 0
        ev_info.max_dist = 1000
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = 'AK,CN,II,IU,US,XM,XV,XZ,YV,TA'  # note: cannot use '*' because of IM
        ev_info.channel = 'BH?,HH?'
        #ev_info.channel = 'BN?,HN?,EN?'    # strong motion
        ev_info.use_catalog = 0 
        ev_info.elat = 56.0382
        ev_info.elon = -149.2686
        ev_info.edep = 11000
        ev_info.emag = 5.3
        ev_info.resample_freq = 50         # for CAP
        ev_info.scale_factor = 100         # for CAP

#---------------------------------------------------------------
# A cluster of events occurred on the Castle Mountain fault after the Kodiak event
# Largest event in this cluster is M4.3 (depth 3km)
    if iex == 100:
        ev_info.use_catalog = 1
        ev_info.otime = obspy.UTCDateTime("2018-01-26 21:10:04")
        ev_info.min_dist = 0
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 50
        ev_info.tafter_sec = 300
        ev_info.network = 'AK,AT,CN,II,IU,US,XM,XV,XZ,YV,TA'  # note: cannot use '*' because of IM
        ev_info.channel = 'BH?,HH?'
        #ev_info.channel = 'BN?,HN?,EN?'    # strong motion
        ev_info.use_catalog = 0 
        ev_info.elat = 61.7511
        ev_info.elon = -148.1374
        ev_info.edep = 3500
        ev_info.emag = 4.3
        ev_info.resample_freq = 50         # for CAP
        ev_info.scale_factor = 100         # for CAP

        # parameters for examining step response
        #ev_info.ifFilter = True
        #ev_info.filter_type = 'lowpass'
        #ev_info.f1 = 1/10
        #ev_info.zerophase = False
        #ev_info.remove_response = False
        #ev_info.ipre_filt = 0
        #ev_info.demean = False
        #ev_info.detrend = False
        #ev_info.resample_TF = False        # if False then resample_freq is taken from SAC files
        #ev_info.scale_factor = 1           # no scale factor
        #ev_info.rotateRTZ = False
        #ev_info.rotateUVW = True
        #ev_info.isave_raw = True

    return(ev_info)
