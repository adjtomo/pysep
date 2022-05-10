import obspy
import read_event_obspy_file as reof
from getwaveform import *

def get_ev_info(ev_info,iex):
# ===============================================================
# Vipul Silwal: 6 moment tensors in Tape et al. (2018 Nature Geoscience) 
# python run_getwaveform.py event_input_nennuc 0
    if iex == 0:
        ev_info.use_catalog = 0
        # -------------------------------------------------
        # 2012 VLFE + EQ (from EPSL 2013 paper)
        # Note: No PS stations when extracting through IRIS
        #ev_info.otime = obspy.UTCDateTime("2012-04-11T09:21:57.444")
        #ev_info.elat = 64.9222
        #ev_info.elon = -148.9461
        #ev_info.edep = 16000
        #ev_info.emag = 3.80
        # -------------------------------------------------
        # 2013 VLFE
        #ev_info.otime = obspy.UTCDateTime("2013-03-12T07:39:50.214")
        #ev_info.elat = 64.7161
        #ev_info.elon = -148.9505
        #ev_info.edep = 23000
        #ev_info.emag = 3.5
        # -------------------------------------------------
        # 2015 VLFE
        #ev_info.otime = obspy.UTCDateTime("2015-09-12T03:25:12.711")
        #ev_info.elat = 65.1207
        #ev_info.elon = -148.6646
        #ev_info.edep = 21000
        #ev_info.emag = 3.80
        # -------------------------------------------------
        # 2015 EQ
        #ev_info.otime = obspy.UTCDateTime("2015-10-22T13:16:15.794")
        #ev_info.elat = 64.7334
        #ev_info.elon = -149.0388
        #ev_info.edep = 18000
        #ev_info.emag = 2.6
        # -------------------------------------------------
        # 2015 EQ
        #ev_info.otime = obspy.UTCDateTime("2015-10-31T02:56:35.572")
        #ev_info.elat = 64.4285
        #ev_info.elon = -149.6969
        #ev_info.edep = 25000
        #ev_info.emag = 3.40
        # -------------------------------------------------
        # 2016 VLFE + EQ
        ev_info.otime = obspy.UTCDateTime("2016-01-14T19:04:10.727")
        ev_info.elat = 64.6827
        ev_info.elon = -149.2479
        ev_info.edep = 17000
        ev_info.emag = 3.70
        # Output pole-zero file (Chen Ji needed this for source estimation)
        #ev_info.ifsave_sacpaz = True
        # -------------------------------------------------
        
        ev_info.min_dist = 0
        ev_info.max_dist = 200
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = 'AV,CN,AT,TA,AK,XV,II,IU,US'
        # ev_info.network = 'XV,AK'
        ev_info.channel = 'BH?,HH?'
        # for CAP
        ev_info.scale_factor = 100
        ev_info.resample_freq = 50

# 2016 pre-earthquake VLFE in Tape et al. (2018 Nature Geoscience) 
# python run_getwaveform.py event_input_nennuc 1
    if iex == 1:
        ev_info.overwrite_ddir = 1 
        ev_info.ifverbose = True    # debugging output

        ev_info.use_catalog = 0
        ev_info.otime = obspy.UTCDateTime("2016-01-14T19:04:10.727")
        ev_info.elat = 64.6827
        ev_info.elon = -149.2479
        ev_info.edep = 17000
        ev_info.emag = 3.70
        
        ev_info.min_dist = 0
        ev_info.max_dist = 30
        ev_info.tbefore_sec = 2000
        ev_info.tafter_sec = 2000
        ev_info.network = 'AV,CN,AT,TA,AK,XV,II,IU,US'
        ev_info.channel = 'BH?,HH?'
        #ev_info.network = 'XV'
        #ev_info.station = 'FTGH'
        ev_info.resample_TF = False
        ev_info.scale_factor = 0

        # to investigate step response
        ev_info.isave_raw = True
        ev_info.isave_raw_processed = True  # for comparison with Figs S10-S12
        ev_info.isave_ENZ = False
        ev_info.rotateENZ = False    
        ev_info.rotateRTZ = False
        #ev_info.rotateUVW = True
        
        # filter options
        ev_info.ifFilter = True
        ev_info.zerophase = False    # causal
        ev_info.filter_type = 'bandpass'
        ev_info.f1 = 1/100    # fmin
        ev_info.f2 = 1/20     # fmax
        ev_info.corners = 4   # number of poles in sac
        ev_info.remove_response = True
        #ev_info.water_level = 100000   # defaut = 60
        ev_info.outformat = 'VEL'
        #ev_info.ifsave_sacpaz = True
        ev_info.ipre_filt = 1
        ev_info.demean = True
        ev_info.detrend = True
        ev_info.taper = 0.05    # 0.05 default in process_data_new.pl

        # seismograms for CAP
        if 0:
            #ev_info.otime = obspy.UTCDateTime("2016-01-14T19:03:50.727")  # 20 s earlier
            ev_info.isave_ENZ = False
            ev_info.rotateENZ = True    
            ev_info.rotateRTZ = True
            ev_info.resample_TF = True
            ev_info.scale_factor = 100
            ev_info.resample_freq = 50
            ev_info.isave_raw_processed = True 
            ev_info.output_cap_weight_file = False

    return(ev_info)
