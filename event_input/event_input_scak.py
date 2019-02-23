import obspy
import read_event_obspy_file as reof
from getwaveform import *

def get_ev_info(ev_info,iex):
# ===============================================================
# SilwalTape2016 example event (Anchorage)
    if iex == 0:
        ev_info.use_catalog = 1
        ev_info.otime = obspy.UTCDateTime("2009-04-07T20:12:55.351")
        ev_info.min_dist = 0 
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300
        ev_info.network = 'AK,AT,AV,CN,II,IU,US,XM,XV,XZ,YV'  # note: cannot use '*' because of IM
        ev_info.channel = 'BH?'
        ev_info.use_catalog = 0 
        ev_info.elat = 61.45420
        ev_info.elon = -149.7428
        ev_info.edep = 33033.60
        # ev_info.rlat = 61.45420
        # ev_info.rlon = -149.7428
        # ev_info.rtime = obspy.UTCDateTime("2009-04-07T20:12:55.351")
        ev_info.emag = 4.6
        ev_info.resample_freq = 50
        ev_info.scale_factor = 100         # no scale factor

# Iniskin earthquake
# NOTE: must enter username and password above to get SALMON (ZE) stations
    if iex == 1:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2016-01-24T10:30:37.400") 
        ev_info.elat = 59.75
        ev_info.elon = -153.27
        ev_info.edep = 110700
        ev_info.emag = 7.1
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 800
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US'  # IM will probably crash it
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 0        # no resampling
        ev_info.scale_factor = 1         # no scale factor
        
        # parameters for examining step response (causal low-pass filter on raw waveforms)
        # delete AUQ, SPCP, SPBG
        ev_info.rotateRTZ = False
        ev_info.ifFilter = True
        ev_info.filter_type = 'lowpass'
        ev_info.f1 = 1/4
        ev_info.zerophase = False
        ev_info.remove_response = False
        ev_info.ipre_filt = 0
        ev_info.demean = False
        ev_info.detrend = False

        ev_info.user = None
        ev_info.password = None

# Iniskin earthquake - all strong motion
    if iex == 2:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2016-01-24T10:30:37.400") 
        ev_info.elat = 59.75
        ev_info.elon = -153.27
        ev_info.edep = 110700
        ev_info.emag = 7.1
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 800
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = '*'
        ev_info.channel = 'BN?,HN?,EN?'
        ev_info.resample_freq = 0        # no resampling
        ev_info.scale_factor = 1         # no scale factor

# Totschunda fault
    if iex == 3:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2017-02-13T07:17:12.000") 
        ev_info.elat = 62.5154
        ev_info.elon = -142.7485
        ev_info.edep = 7564
        ev_info.emag = 5.3
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 50
        ev_info.tafter_sec = 300
        ev_info.network = 'AV'
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor

# Klukwan earthquakes
    if iex == 4:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2017-05-01T12:31:53.000") 
        ev_info.elat = 59.8522
        ev_info.elon = -136.6618
        ev_info.edep = 4000
        ev_info.emag = 6.2
        ev_info.otime = obspy.UTCDateTime("2017-05-01T14:18:14.000") 
        ev_info.elat = 59.8184
        ev_info.elon = -136.7163
        ev_info.edep = 1000
        ev_info.emag = 6.0
        
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 500
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor

# Kantishna earthquakes
    if iex == 5:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2017-04-29T11:15:48.000") 
        ev_info.elat = 63.1296
        ev_info.elon = -151.1517
        ev_info.edep = 11000
        ev_info.emag = 5.2
        ev_info.otime = obspy.UTCDateTime("2017-01-31T09:38:37.000") 
        ev_info.elat = 63.0817
        ev_info.elon = -150.9427
        ev_info.edep = 135000
        ev_info.emag = 5.2
        
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 500
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor

# Cook Inlet earthquake
    if iex == 6:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2017-05-07T04:25:19.000") 
        ev_info.elat = 60.1828
        ev_info.elon = -151.6803
        ev_info.edep = 41400
        ev_info.emag = 5.2
    
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 500
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor
        ev_info.user = ''
        ev_info.password = ''
            
# EXAMPLE TEMPLATE -- USE THIS ONLY FOR QUICK TESTING
# This is a template for testing before creating an example.
# We can delete this if it creates more issues.
    if iex == 7:
        ev_info.otime = obspy.UTCDateTime("2009-04-07T20:12:55")
        ev_info.min_dist = 0 
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300
        ev_info.network = 'II,IU'
        ev_info.station = 'KDAK,COLA'
        ev_info.channel = '*'

# Cook Inlet earthquake
    if iex == 8:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2017-05-07T04:25:19.000") 
        ev_info.elat = 60.1945
        ev_info.elon = -151.6743
        ev_info.edep = 64000
        ev_info.emag = 5.2
    
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 500
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor
      
# South-east unusual event
    if iex == 9:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2016-12-30T18:27:10.406") 
        ev_info.elat = 58.6716
        ev_info.elon = -134.7691
        ev_info.edep = 274
        ev_info.emag = 3.43
    
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 500
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor
        ev_info.ifFilter = True
        ev_info.filter_type = 'bandpass'
        ev_info.f1 = 1/100  # fmin
        ev_info.f2 = 1/10  # fmax
        ev_info.corners = 4

# Prince William Sound event
    if iex == 10:
        ev_info.use_catalog = 1
        ev_info.otime = obspy.UTCDateTime("2017-11-27T22:18:30.000")
        ev_info.min_dist = 0 
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300
        ev_info.network = 'AK,AT,AV,CN,II,IU,US,XM,XV,XZ,YV'  # note: cannot use '*' because of IM
        ev_info.channel = 'BH?'
        ev_info.use_catalog = 0 
        ev_info.elat = 60.5634
        ev_info.elon = -147.4192
        ev_info.edep = 16000
        ev_info.emag = 5.3
        ev_info.resample_freq = 50
        ev_info.scale_factor = 100         # no scale factor
    
    if iex ==11:
        ev_info.use_catalog = 1
        ev_info.otime = obspy.UTCDateTime("2017-11-27T22:18:30.000")
        ev_info.min_dist = 0 
        ev_info.max_dist = 300
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 200
        ev_info.network = 'AK,AT,AV,CN,II,IU,US,XM,XV,XZ,YV'  
        ev_info.channel = 'BH?'
        ev_info.use_catalog = 0 
        ev_info.elat = 60.5634
        ev_info.elon = -147.4192
        ev_info.edep = 16000
        ev_info.emag = 5.3
        ev_info.resample_freq = 50
        ev_info.scale_factor = 100         # no scale factor
        
        ev_info.phase_window = True
        ev_info.phases = ["P","S"] # Choosing time period with respect to P & S
        #ev_info.phases = ["P","P"]
        ev_info.write_sac_phase = True
        ev_info.taupmodel = "ak_scak" 

    if iex == 12:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2018-11-30T17:29:29.000") 
        ev_info.elat = 61.3234
        ev_info.elon = -149.9234
        ev_info.edep = 43000 
        ev_info.emag = 7.0
    
        # subset of stations

        ev_info.min_dist = 0
        #ev_info.min_dist = 300
        ev_info.max_dist = 500 
        ev_info.tbefore_sec = 100 
        ev_info.tafter_sec = 1000 
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US,DE' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor
        ev_info.user = ''
        ev_info.password = ''
    return(ev_info)
#=================================================================================
# END EXAMPLES
#=================================================================================
