import obspy
import read_event_obspy_file as reof
from getwaveform import *

def get_ev_info(ev_info,iex):
# ===============================================================
# SilwalTape2016 example event (Anchorage) -- python run_getwaveform.py event_input_mtuq2022 1
    if iex == 1:
        ev_info.overwrite_ddir = 1
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        ev_info.otime = obspy.UTCDateTime("2009-04-07T20:12:55.351")
        ev_info.min_dist = 0
        ev_info.max_dist = 300
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300

        # RAW and ENZ files can be used when checking if you are receiving all possible data (example station: SOLD)
        ev_info.isave_raw = False
        ev_info.isave_raw_processed = False
        ev_info.isave_ENZ = False

        # Networks of interest (All within distance limitations)
        ev_info.network = 'AK,YV'
        ev_info.channel = 'BH?'

        # Event information
        ev_info.elat = 61.45420
        ev_info.elon = -149.7428
        ev_info.edep = 33033.60
        ev_info.emag = 4.6

        # scaling and resampling
        ev_info.resample_TF = True
        ev_info.resample_freq = 50

        # Scaling depends on units (This assumes velocity in units cm/s)
        ev_info.scale_factor = 100

# 2020 Southern California event
    if iex == 2:
        ev_info.overwrite_ddir = 1
        ev_info.use_catalog = 0
        ev_info.otime = obspy.UTCDateTime("2020-04-04T01:53:18.920")
        ev_info.min_dist = 0
        ev_info.max_dist = 300
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300

        # RAW and ENZ files can be used when checking if you are receiving all possible data
        ev_info.isave_raw = False
        ev_info.isave_raw_processed = False
        ev_info.isave_ENZ = False

        # Only receive CI network stations
        ev_info.network = 'CI'
        ev_info.channel = 'BH?'
        # CI.SWS causes script to crash (Station problems?)
        ev_info.station = 'BAR,IKP,PLM,GLA,BC3,PDM,IRM,DAN,GMR,TUQ,HEC,GSC,RRX,BBR,SCI2,CIA,SDD,VTV,ADO,ARV,DGR,SVD,DJJ,FMP,-SWS' #Receive subset of stations
        
        # Event specific information
        ev_info.elat = 33.490
        ev_info.elon = -116.506
        ev_info.edep =  10500.0
        ev_info.emag = 4.9

        # scaling and resampling
        ev_info.resample_TF = True
        ev_info.resample_freq = 50
        # See iex == 1 for more info
        ev_info.scale_factor = 100

# 2017 North Korea Event
    if iex == 3:
        ev_info.overwrite_ddir = 1
        ev_info.use_catalog = 0
        ev_info.otime = obspy.UTCDateTime("2017-09-03 03:30:01.760")
        ev_info.min_dist = 0
        ev_info.max_dist = 1300
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        #^^^^^^^^^^^^^^^^^^^^^^^^^^

        # RAW and ENZ files can be used when checking if you are receiving all possible data
        ev_info.isave_raw = False
        ev_info.isave_raw_processed = False
        ev_info.isave_ENZ = False

        # Network and Channel requests CHANGE THIS vvvvvvvv
        ev_info.network = 'IC,IU,G,JP'
        ev_info.channel = 'BH?'
        ev_info.station = 'MDJ,INCN,JTU,MAJO,INU,BJT,YSS'
        ev_info.location = '00'
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
        # Event specific information CHANGE THIS vvvvvvvvvvvv
        ev_info.elat = 41.3324
        ev_info.elon = 129.0297
        ev_info.edep =  1000.0
        ev_info.emag = 5.18
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        # scaling and resampling
        ev_info.resample_TF = True
        ev_info.resample_freq = 5
        # See iex == 1 for more info
        ev_info.scale_factor = 100

# Iceland Event
    if iex == 4:
        ev_info.overwrite_ddir = 1
        ev_info.use_catalog = 0
        #CHANGE THIS vvvvvvvvv
        #
        ev_info.otime = obspy.UTCDateTime("2014-08-25T16:19:03.0")
        ev_info.min_dist = 20
        ev_info.max_dist = 300
        ev_info.tbefore_sec = 60
        ev_info.tafter_sec = 360
        #^^^^^^^^^^^^^^^^^^^^^^^^^^

        # RAW and ENZ files can be used when checking if you are receiving all possible data
        ev_info.isave_raw = False
        ev_info.isave_raw_processed = False
        ev_info.isave_ENZ = False

        # Network and Channel requests CHANGE THIS vvvvvvvv
        ev_info.network = 'Z7'
        ev_info.channel = 'HH?'
        ev_info.station = 'RODG,DYSA,LIND,LOKT,LAUF,KALF,HELI,FAG,SVIN,K250' #Receive all stations except CI.SWS (if not needed just delete)
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
        # Event specific information CHANGE THIS vvvvvvvvvvvv
        ev_info.elat = 64.612
        ev_info.elon = -17.472
        ev_info.edep =  5000.0
        ev_info.emag = 4.6
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        # scaling and resampling
        ev_info.resample_TF = True
        ev_info.resample_freq = 50
        # See iex == 1 for more info
        ev_info.scale_factor = 100

    return(ev_info)
