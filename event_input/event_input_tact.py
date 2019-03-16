import obspy
import read_event_obspy_file as reof
from getwaveform import *

def get_ev_info(ev_info,iex):
# ===============================================================
# southern Minto Flats fault
    if iex == 0:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 1          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2016-11-06T09:29:10.000") 
        #ev_info.elat = 60.1945
        #ev_info.elon = -151.6743
        #ev_info.edep = 64000
        #ev_info.emag = 5.2
    
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 500
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor

    if iex == 1:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 1          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        #ev_info.otime = obspy.UTCDateTime("2008-09-25T19:34:10.05") 
        ev_info.otime = obspy.UTCDateTime("2007-02-17T14:51:40.16") 
    
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 300
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 500
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor

# November 8, 2017 Minto Flats earthquake M3.7
    if iex == 2:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2017-11-08T06:49:11") # AEC prelim
        ev_info.elat = 64.8632
        ev_info.elon = -148.6255
        ev_info.edep = 16000
        ev_info.emag = 3.6
    
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 300
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 500
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor

    if iex == 3:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2017-12-30T11:43:16") # AEC prelim
        ev_info.elat = 63.8200
        ev_info.elon = -149.0497
        ev_info.edep = 6000
        ev_info.emag = 4.1
    
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 300
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 500
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor        
   
    # M5.1 W. of Tanana 
    if iex == 4:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2019-01-06T03:45:34") # AEC prelim
        ev_info.elat = 65.4070
        ev_info.elon = -153.2799
        ev_info.edep = 17000
        ev_info.emag = 5.1 
    
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 300
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 500

        # For DE (Nanometrics) data
        ev_info.user = 'kksmith7@alaska.edu' 
        ev_info.password = 'wmpo3NmqTcRm' 
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US,DE' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor        
    
    # M3.7 North Pole 
    if iex == 5:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2019-03-09T23:39:58") # AEC prelim
        ev_info.elat = 64.5830 
        ev_info.elon = -147.7238
        ev_info.edep = 20600 
        ev_info.emag = 3.7 
    
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 200  
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 500

        # For DE (Nanometrics) data
        #ev_info.user = 'kksmith7@alaska.edu' 
        #ev_info.password = 'wmpo3NmqTcRm' 
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US,DE' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor        

    return(ev_info)
#=================================================================================
# END EXAMPLES
#=================================================================================
