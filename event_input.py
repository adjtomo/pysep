import obspy
import read_event_obspy_file as reof
from getwaveform_new import *

def get_ev_info(ev_info,iex):
# ===============================================================
# EXAMPLE TEMPLATE -- USE THIS ONLY FOR QUICK TESTING
# This is a template for testing before creating an example.
# We can delete this if it creates more issues.
    if iex == 0:
        ev_info.otime = obspy.UTCDateTime("2009-04-07T20:12:55")
        ev_info.min_dist = 0 
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300
        ev_info.network = 'II,IU'
        ev_info.station = 'KDAK,COLA'
        ev_info.channel = '*'
        
#=================================================================================
# CATEGORY B EXAMPLES: important events
# 1XX: southern Alaska
# 2XX: central Alaska
# 9XX: other
#=================================================================================
# SilwalTape2016 example event (Anchorage)
    if iex == 100:
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

# Iniskin earthquake
# NOTE: must enter username and password above to get SALMON (ZE) stations
    if iex == 101:
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
    if iex == 103:
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
    if iex == 104:
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
    if iex == 105:
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
    if iex == 106:
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
    if iex == 107:
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
            
# MFFZ earthquakes for investigating the step response
# LISTED IN CHRONOLOGICAL ORDER
    if iex == 200:
        ev_info.idb = 1
        ev_info.use_catalog = 0
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2014-08-31T03:06:57.111")
        ev_info.elat = 65.1526
        ev_info.elon = -149.0398
        ev_info.edep = 16614.7
        ev_info.emag = 5.20
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2014-10-21T00:36:58.333")
        ev_info.elat = 65.1489
        ev_info.elon = -149.0413
        ev_info.edep = 13134.8
        ev_info.emag = 4.90
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2014-10-23T16:30:23.968")
        ev_info.elat = 65.1644
        ev_info.elon = -149.0523
        ev_info.edep = 20066.5
        ev_info.emag = 5.00
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2015-09-12T03:25:12.711")
        ev_info.elat = 65.1207
        ev_info.elon = -148.6646
        ev_info.edep = 15556.8
        ev_info.emag = 2.63
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2015-10-31T02:56:35.572")
        ev_info.elat = 64.4285
        ev_info.elon = -149.6969
        ev_info.edep = 23852.1
        ev_info.emag = 3.47
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2016-01-14T19:04:10.727")
        ev_info.elat = 64.6827
        ev_info.elon = -149.2479
        ev_info.edep = 22663.7
        ev_info.emag = 3.80
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2016-11-06T9:29:10.579")
        ev_info.elat = 64.1639
        ev_info.elon = -150.0626
        ev_info.edep = 23190.0
        ev_info.emag = 4.00
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2016-12-08T10:18:13.868")
        ev_info.elat = 64.1937
        ev_info.elon = -150.0376
        ev_info.edep = 24522.1
        ev_info.emag = 4.30
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2017-05-08T05:09:02.000") 
        ev_info.elat = 65.2643
        ev_info.elon = -146.922
        ev_info.edep = 9000   # AEC/Vipul
        ev_info.emag = 3.60   # Vipul
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2017-06-28T12:58:52.000") 
        ev_info.elat = 64.7569
        ev_info.elon = -148.8883
        ev_info.edep = 18000  # Vipul
        ev_info.emag = 3.50   # Vipul
        # -------------------------------------------------
        # VIPUL: WHAT ARE THESE? OTHER EVENTS?
        # ev_info.otime = obspy.UTCDateTime("2015-03-30T12:33:19.000")
        # ev_info.otime = obspy.UTCDateTime("2015-10-20T19:14:16.000")
        # ev_info.otime = obspy.UTCDateTime("2011-12-21T16:28:41.000")
        # -------------------------------------------------
        
        ev_info.min_dist = 0
        ev_info.max_dist = 300
        ev_info.tbefore_sec = 200
        ev_info.tafter_sec = 600
        ev_info.network = 'AV,CN,AT,TA,AK,XV,II,IU,US'
        # ev_info.network = 'XV,AK'
        ev_info.channel = 'BH?,HH?'
        ev_info.scale_factor = 1
        ev_info.resample_freq = 0
        # for CAP
        # ev_info.scale_factor = 10.0**2
        # ev_info.resample_freq = 50 
        
        # to investigate step response
        ev_info.rotateRTZ = True
        ev_info.rotateUVW = True
        ev_info.ifFilter = True
        ev_info.zerophase = False    # causal
        # filter_type = 'lowpass'
        ev_info.filter_type = 'bandpass'
        ev_info.f1 = 1/100  # fmin
        ev_info.f2 = 1/10  # fmax
        ev_info.corners = 4
        # ev_info.remove_response = False
        ev_info.remove_response = True
        ev_info.ipre_filt = 1
        ev_info.demean = True
        ev_info.detrend = True
        ev_info.output_cap_weight_file = False
        # ev_info.outformat = 'DISP'
        ev_info.ifsave_sacpaz = True
        ev_info.taper = 0.2
        
# NENNUC event (from Steve)
    if iex == 201:
        ev_info.idb = 1
        ev_info.use_catalog = 0          # manually enter AEC catalog parameters
        ev_info.otime = obspy.UTCDateTime("2016-01-14T19:04:10.727")
        ev_info.elat = 64.6827
        ev_info.elon = -149.2479
        ev_info.edep = 22663.7
        ev_info.emag = 3.8
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2015-09-12T03:25:12.711")
        ev_info.elat = 65.1207
        ev_info.elon = -148.6646
        ev_info.edep = 15556.8
        ev_info.emag = 2.6
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2013-03-12T07:39:50.214")
        ev_info.elat = 64.7161
        ev_info.elon = -148.9505
        ev_info.edep = 20000
        ev_info.emag = 2.1
            
        # For CAP
        ev_info.min_dist = 0
        ev_info.max_dist = 200
        ev_info.tbefore_sec = 500
        ev_info.tafter_sec = 2000
        ev_info.taper = 0.1
        ev_info.network = 'AV,CN,AT,TA,AK,XV,II,IU,US'
        # ev_info.network = 'AK,TA,II,IU,US'
        # ev_info.station = "-RAG"  # RAG has choppy data; gives error: Exception: Can't merge traces with same ids but differing sampling rates!
        ev_info.channel = 'BH?,HH?'
        ev_info.pre_filt = (0.0025, 0.003, 10.0, 15.0) # BH default
        ev_info.ipre_filt = 1

        # to investigate step response
        # ev_info.ev_info.tbefore_sec = 2000
        # ev_info.tafter_sec = 2000
        # ev_info.resample_freq = 0 
        # ev_info.scale_factor = 1
        # ev_info.ifFilter = True
        # ev_info.zerophase = False
        # ev_info.filter_type = 'bandpass'
        # ev_info.f1 = 1/100
        # ev_info.f2 = 1/20
        # ev_info.remove_response = True
        # ev_info.ipre_filt = 0
        # ev_info.demean = True
        # ev_info.detrend = True
        
# Kyle Nenana basin earthquakes for basin amplification
    if iex == 210:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # AEC source parameters
        ev_info.otime = obspy.UTCDateTime("2015-11-20T10:53:48.168") 
        ev_info.elat = 64.6210
        ev_info.elon = -149.4024
        ev_info.edep = 17113.4
        ev_info.emag = 2.67
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 200
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 400
        ev_info.network = 'AK,AT,II,IU,US,XM,XV,XZ,TA'  # no CN,AV,YV,ZE
        ev_info.channel = 'BH?,HH?'
        # ev_info.resample_freq = 0        # no resampling
        # ev_info.scale_factor = 1         # no scale factor
        # For CAP moment tensor
        ev_info.resample_freq = 50         # same as greens function 
        ev_info.scale_factor = 100         # change from m/s to cm/s
        # ev_info.ipre_filt = 0
        ev_info.remove_response = True
        # ev_info.demean = False
        # ev_info.detrend = False

# gets waveforms from M 8.3 Chile event with stations centered in Minto
    if iex == 211:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0
        ev_info.otime = obspy.UTCDateTime("2015-09-16T22:54:33.000")
        ev_info.elat = -31.5695
        ev_info.elon = -71.6543
        ev_info.edep = 22400.0
        ev_info.emag = 8.30
        ev_info.rlat = 64.7716
        ev_info.rlon = -149.1465
        ev_info.rtime = obspy.UTCDateTime("2015-09-16T23:09:15.000")
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 200
        ev_info.min_dist = 0
        ev_info.max_dist = 100
        ev_info.network = 'AK,AT,II,IU,US,XM,XV,XZ,TA'  # no CN,AV,YV,ZE
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 0        
        ev_info.scale_factor = 1        
        ev_info.remove_response = True

# gets a ???
    if iex == 212:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0
        ev_info.otime = obspy.UTCDateTime("2016-06-06T00:00:00.000")
        ev_info.elat = 64.6130
        ev_info.elon = -149.0992
        ev_info.edep = 0
        ev_info.emag = 0.00
        # ev_info.rlat = 64.7716
        # ev_info.rlon = -149.1465
        # ev_info.rtime = obspy.UTCDateTime("2015-09-16T23:09:15.000")
        ev_info.tbefore_sec = 0
        ev_info.tafter_sec = 3600
        ev_info.min_dist = 0
        ev_info.max_dist = 100
        ev_info.network = 'XV,AK,TA'  # no CN,AV,YV,ZE
        ev_info.channel = 'HH?,BH?'
        ev_info.resample_freq = 0        
        ev_info.scale_factor = 1        
        ev_info.remove_response = True
        ev_info.rotateRTZ = False
        # ev_info.pre_filt=(f0*0.001, f1*0.001, f2*1000, f3*1000)
        # ev_info.ipre_filt = 2
        ev_info.ipre_filt = 0

# Chatnika earthquake
    if iex == 213:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2017-05-08T05:09:02.000") 
        ev_info.elat = 65.2643
        ev_info.elon = -146.922
        ev_info.edep = 9000
        ev_info.emag = 3.8
        
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 500
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor
            
# NE Nenana earthquake
    if iex == 214:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2017-06-28T12:58:52") 
        ev_info.elat = 64.7569
        ev_info.elon = -148.8883
        ev_info.edep = 18000
        ev_info.emag = 3.5
        
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 500
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        
        ev_info.scale_factor = 100 

# Loop over multiple event       
    if iex == 215:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        
        events_file = "/home/ksmith/REPOSITORIES/manuscripts/kyle/papers/basinamp/data/basinamp_obspy.txt"
        eids,otimes,elons,elats,edeps,emags = reof.read_events_obspy_file(events_file)

        ev_info_list = []
        for xx in range(0,3):
            ev_info_temp = ev_info.copy()
            ev_info_temp.otime = obspy.UTCDateTime(otimes[xx])
            ev_info_temp.elat = elats[xx]
            ev_info_temp.elon = elons[xx]
            ev_info_temp.edep = edeps[xx]
            ev_info_temp.emag = emags[xx]
            ev_info_temp.eid = eids[xx]
            
            # subset of stations
            ev_info_temp.min_dist = 0
            ev_info_temp.max_dist = 500
            ev_info_temp.tbefore_sec = 100
            ev_info_temp.tafter_sec = 500
            ev_info_temp.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
            ev_info_temp.channel = 'BH?,HH?'
            ev_info_temp.resample_freq = 50        
            ev_info_temp.scale_factor = 100 

            # append getwaveform objects
            ev_info_list.append(ev_info_temp)
        
        # always return ev_info
        ev_info = ev_info_list
#----------------------------------------------------------

#-----------------------------------------------------------
# Event from HutchisonGhosh2016
# Crashes when using all networks [i.e. network = '*']
    if iex == 900:
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

    return(ev_info)
#=================================================================================
# END EXAMPLES
#=================================================================================
