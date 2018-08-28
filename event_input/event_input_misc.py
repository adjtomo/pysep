import obspy
import read_event_obspy_file as reof
from getwaveform import *

def get_ev_info(ev_info,iex):
#-----------------------------------------------------------
# Event from HutchisonGhosh2016
# Crashes when using all networks [i.e. network = '*']
    if iex == 0:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2014-12-12T00:13:17") 
        ev_info.elat = 49.10
        ev_info.elon = -124.05
        ev_info.edep = 60000
        ev_info.emag = 4.10
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 250
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300
        ev_info.network = '*'
        ev_info.channel = 'LH?,BH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor
        
        # For plotting filtered waveforms
        ev_info.tbefore_sec = 500
        ev_info.tafter_sec = 500
        ev_info.resample_freq = 0 
        ev_info.scale_factor = 1
        ev_info.ifFilter = True
        ev_info.zerophase = True
        ev_info.filter_type = 'bandpass'
        ev_info.f1 = 1/50
        ev_info.f2 = 1/20
        ev_info.remove_response = True
        ev_info.ipre_filt = 2
        ev_info.demean = True
        ev_info.detrend = True
        ev_info.taper = True

# Kashmir earthquake
    if iex == 1:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 1          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2005-10-08T03:50:39") 
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 1000
        ev_info.tbefore_sec = 300
        ev_info.tafter_sec = 500
        ev_info.network = 'II'
        ev_info.channel = 'LH?,BH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor

# Gulf of Alaska M3.0
    if iex == 2:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2017-10-16T22:00:31") 
        ev_info.elat = 57.6160
        ev_info.elon = -142.7699
        ev_info.edep = 6000
        ev_info.emag = 3.0
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 50
        ev_info.tafter_sec = 600
        ev_info.network = 'AK'
        ev_info.network = 'AK,AT,AV,CN,II,IU,US,XM,TA,XE,XR,XZ,YV,XV,ZE,XG'
        ev_info.channel = 'BH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100        # no scale factor

# Earthquake in Fairbanks
    if iex == 3:
        ev_info.ifverbose = False
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters

        # subset waveforms
        ev_info.min_dist = 0 
        ev_info.max_dist = 200
        ev_info.tbefore_sec = 50
        ev_info.tafter_sec = 300
        ev_info.resample_freq = 50
        ev_info.scale_factor = 100
        ev_info.network = 'AK,AT,CN,II,IU,US,XM,TA,XE,XR,XZ,YV,XV,ZE,XG'
        ev_info.channel = 'BH?,HH?'
        ev_info.otime = obspy.UTCDateTime("2018-07-24T10:32:05") 
        ev_info.elat = 64.8471
        ev_info.elon = -147.7885
        ev_info.edep = 8000
        ev_info.emag = 3.2
        ev_info.rotateUVW = True 

        ev_info.user = None
        ev_info.password = None

# Earthquake on North Slope
    if iex == 4:
        ev_info.ifverbose = False
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters

        # subset waveforms
        ev_info.min_dist = 0 
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 50
        ev_info.tafter_sec = 300
        ev_info.resample_freq = 50
        ev_info.scale_factor = 100
        ev_info.network = 'AK,AT,CN,II,IU,US,XM,TA,XE,XR,XZ,YV,XV,ZE,XG'
        ev_info.channel = 'BH?,HH?'
        ev_info.otime = obspy.UTCDateTime("2018-08-12T14:58:54") 
        ev_info.elat = 69.6239
        ev_info.elon = -145.2468
        ev_info.edep = 10000
        ev_info.emag = 6.4
        ev_info.rotateUVW = True 

    return(ev_info)
#=================================================================================
# END EXAMPLES
#=================================================================================
