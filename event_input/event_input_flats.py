import obspy
import read_event_obspy_file as reof
from getwaveform import *

def get_ev_info(ev_info,iex):
# ===============================================================            
# MFFZ earthquakes for investigating the step response
# python run_getwaveform.py event_input_flats 0
    if iex == 0:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        
        ievent = 13
        events_file = "/home/carltape/REPOSITORIES/manuscripts/carltape/papers/nennuc/clipping/data/MFFZ_step_response_obspy.txt"
        eids,otimes,elons,elats,edeps,emags = reof.read_events_obspy_file(events_file)

        ev_info_list = []
        for xx in range(ievent-1,ievent):
        #for xx in range(len(eids)):
            ev_info_temp = ev_info.copy()
            ev_info_temp.otime = obspy.UTCDateTime(otimes[xx])
            ev_info_temp.elat = elats[xx]
            ev_info_temp.elon = elons[xx]
            ev_info_temp.edep = edeps[xx]
            ev_info_temp.emag = emags[xx]
            ev_info_temp.eid = eids[xx]
            
            ev_info_temp.min_dist = 0
            ev_info_temp.max_dist = 300
            ev_info_temp.tbefore_sec = 200
            ev_info_temp.tafter_sec = 600
            #ev_info_temp.network = 'AV,CN,AT,TA,AK,XV,II,IU,US,DE'
            #ev_info_temp.network = 'AV,CN,AT,TA,AK,XV,II,IU,US'
            # ev_info_temp.network = 'XV,AK'
            ev_info_temp.network = 'XV,DE'
            ev_info_temp.channel = 'BH?,HH?'
            #ev_info_temp.scale_factor = 1
            ev_info_temp.resample_TF = False
            # for CAP
            # ev_info_temp.scale_factor = 10.0**2
            # ev_info_temp.resample_freq = 50 
        
            # For DE (Nanometrics) data
            ev_info_temp.user = 'kksmith7@alaska.edu' 
            ev_info_temp.password = 'wmpo3NmqTcRm' 

            # to investigate step response
            ev_info_temp.isave_raw = True
            ev_info_temp.isave_raw_processed  = True
            ev_info_temp.rotateENZ = False
            ev_info_temp.rotateRTZ = False
            ev_info_temp.rotateUVW = True
            # filter options
            ev_info_temp.ifFilter = True
            ev_info_temp.zerophase = False    # causal
            # filter_type = 'lowpass'
            ev_info_temp.filter_type = 'bandpass'
            ev_info_temp.f1 = 1/100  # fmin
            ev_info_temp.f2 = 1/10  # fmax
            ev_info_temp.corners = 4
            # ev_info_temp.remove_response = False
            ev_info_temp.remove_response = True
            ev_info_temp.ipre_filt = 1
            ev_info_temp.demean = True
            ev_info_temp.detrend = True
            ev_info_temp.output_cap_weight_file = False
            # ev_info_temp.outformat = 'DISP'
            ev_info_temp.ifsave_sacpaz = True
            ev_info_temp.taper = 0.2

            # append getwaveform objects
            ev_info_list.append(ev_info_temp)
        
        # always return ev_info
        ev_info = ev_info_list
        
# NENNUC event (from Steve)
    if iex == 1:
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
    if iex == 2:
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
    if iex == 3:
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
    if iex == 4:
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
    if iex == 5:
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
    if iex == 6:
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

# Loop over multiple event for MT 
    if iex == 7:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        
        events_file = "/home/ksmith/REPOSITORIES/manuscripts/kyle/papers/basinresp/data/basinamp_obspy.txt"
        eids,otimes,elons,elats,edeps,emags = reof.read_events_obspy_file(events_file)

        ev_info_list = []
        for xx in range(26,27):
            ev_info_temp = ev_info.copy()
            ev_info_temp.otime = obspy.UTCDateTime(otimes[xx])
            ev_info_temp.elat = elats[xx]
            ev_info_temp.elon = elons[xx]
            ev_info_temp.edep = edeps[xx]
            ev_info_temp.emag = emags[xx]
            ev_info_temp.eid = eids[xx]
            
            # subset of stations
            ev_info_temp.min_dist = 0
            ev_info_temp.max_dist = 200 #500
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

    if iex == 8:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        
        events_file = "/home/ksmith/REPOSITORIES/manuscripts/kyle/misc/F2TN_events_near_out.txt"
        eids,otimes,elons,elats,edeps,emags = reof.read_events_obspy_file(events_file)
        ev_info_temp.tafter_sec = 500
        ev_info_temp.network = 'AK,XV' 
        ev_info_temp.channel = 'BH?,HH?'
        ev_info_temp.resample_freq = 50        
        ev_info_temp.scale_factor = 100 

        # append getwaveform objects
        ev_info_list.append(ev_info_temp)
        
        # always return ev_info
        ev_info = ev_info_list

    if iex == 9:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2018-05-27T10:40:00") 
        ev_info.elat = 64.613
        ev_info.elon = -149.0992
        ev_info.edep = 000
        ev_info.emag = 0
        ev_info.rotateRTZ = False
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 15
        ev_info.tbefore_sec = 1800
        ev_info.tafter_sec = 3600
        ev_info.network = 'XV' 
        ev_info.channel = 'HH?'
        ev_info.resample_freq = 50        
        ev_info.scale_factor = 100 

    if iex == 10:     
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
	# GCMT source parameters
	# the otime is the centroid time and accounts for tshift
        
	# FNN1 outage 
	#ev_info.otime = obspy.UTCDateTime("2018-07-15T19:30:00") 
        #ev_info.elat = 64.5716
        #ev_info.elon = -149.2179

	# FTGH outage
        #ev_info.otime = obspy.UTCDateTime("2018-07-23T23:30:00")
        #ev_info.elat = 64.6917
        #ev_info.elon = -148.8278
      
	# FNN1 out again 
        #ev_info.otime = obspy.UTCDateTime("2018-07-18T10:00:00")  
        #ev_info.elat = 64.5716
        #ev_info.elon = -149.2179

        # F3TN radio out and on
        #ev_info.otime = obspy.UTCDateTime("2018-08-27T20:30:00")  
        #ev_info.otime = obspy.UTCDateTime("2018-09-02T16:30:00")  
        #ev_info.elat = 64.7713
        #ev_info.elon = -149.1466


        # F4TN radio out and on
        #ev_info.otime = obspy.UTCDateTime("2018-08-04T16:00:00")  
        #ev_info.otime = obspy.UTCDateTime("2018-08-07T01:30:00")  
        #ev_info.elat = 64.8340
        #ev_info.elon = -149.1528
 
        
        # FAPT radio out and on
        #ev_info.otime = obspy.UTCDateTime("2018-05-31T22:30:00")  
        #ev_info.otime = obspy.UTCDateTime("2018-05-20T06:30:00")  
        #ev_info.elat = 64.5498
        #ev_info.elon = -149.0831

	# FPAP out and on
        #ev_info.otime = obspy.UTCDateTime("2018-09-28T12:00:00")  
        #ev_info.otime = obspy.UTCDateTime("2018-10-04T00:00:00")  
        #ev_info.elat = 64.6130
        #ev_info.elon = -149.0991
	
	# F2TN out and on
        ev_info.otime = obspy.UTCDateTime("2018-08-03T18:30:00")  
        ev_info.elat = 64.7091
        ev_info.elon = -149.1327
        
        ev_info.edep = 0 
        ev_info.emag = 0
        ev_info.rotateRTZ = False
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 5 
        ev_info.tbefore_sec = 0 
        ev_info.tafter_sec = 3600 
        ev_info.network = 'XV' 
        ev_info.channel = 'LH?'
        ev_info.resample_freq = 50        
        ev_info.scale_factor = 100 

    if iex == 11:     
        
        # Manley event
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
	
        ev_info.otime = obspy.UTCDateTime("2018-08-28T15:18:44")  
        ev_info.elat = 65.1778
        ev_info.elon = -150.4964
        
        ev_info.edep = 16000
        ev_info.emag = 4.8 
        ev_info.rotateRTZ = True 
    
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 500
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        
        ev_info.scale_factor = 100

    if iex == 12: 
       
        # Manley event
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
	
        ev_info.otime = obspy.UTCDateTime("2018-10-03T03:29:37.544")  
        ev_info.elat = 64.8979 
        ev_info.elon = -148.9191
        
        ev_info.edep = 19690 
        ev_info.emag = 3.94 
        ev_info.rotateRTZ = True 
    
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 500 
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 500
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        
        ev_info.scale_factor = 100


    if iex == 13:     
        
        # Nenana EQ 
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
	
	# from AEC website 
        #ev_info.otime = obspy.UTCDateTime("2018-08-25T18:15:51") 
        #ev_info.elat = 64.6181
        #ev_info.elon = -149.2233
        #ev_info.edep = 18900

        # Reviewed time and location 
        ev_info.otime = obspy.UTCDateTime("2018-08-25T18:15:51.595") 
        ev_info.elat = 64.6224
        ev_info.elon = -149.2107
        ev_info.edep = 21166 
        ev_info.emag = 3.2 
        ev_info.rotateRTZ = True 
    
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 300 
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        
        ev_info.scale_factor = 100

    if iex == 14:     
        
        # Tanana event
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
	
        ev_info.otime = obspy.UTCDateTime("2018-10-27T16:57:27")  
        ev_info.elat = 65.2234
        ev_info.elon = -151.6636
        
        ev_info.edep = 17000
        ev_info.emag = 5.3 
        ev_info.rotateRTZ = True 
    
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 300 
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 400 
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US,DE' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        
        ev_info.scale_factor = 100

    if iex == 16:     
        
        # Anchorage event
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        
        ievent = 1 
        events_file = "/home/ksmith/REPOSITORIES/manuscripts/kyle/posters/2018/AGU_poster.txt"
        eids,otimes,elons,elats,edeps,emags = reof.read_events_obspy_file(events_file)

        ev_info_list = []
        for xx in range(ievent-1,ievent):
        #for xx in range(len(eids)):
            ev_info_temp = ev_info.copy()
            ev_info_temp.otime = obspy.UTCDateTime(otimes[xx])
            ev_info_temp.elat = elats[xx]
            ev_info_temp.elon = elons[xx]
            ev_info_temp.edep = edeps[xx]
            ev_info_temp.emag = emags[xx]
            ev_info_temp.eid = eids[xx]
            
            ev_info_temp.min_az = 0 
            ev_info_temp.max_az = 15 
            ev_info_temp.max_az = 20 
            ev_info_temp.min_dist = 0
            ev_info_temp.max_dist = 600 
            ev_info_temp.tbefore_sec = 50 
            ev_info_temp.tafter_sec = 350 
            #ev_info_temp.network = 'AV,CN,AT,TA,AK,XV,II,IU,US,DE'
            #ev_info_temp.network = 'AV,CN,AT,TA,AK,XV,II,IU,US'
            ev_info_temp.network = 'XV,DE,AK,TA'
            ev_info_temp.channel = 'BH?,HH?'
            #ev_info_temp.scale_factor = 1
            ev_info_temp.resample_TF = False
            # for CAP
            # ev_info_temp.scale_factor = 10.0**2
            # ev_info_temp.resample_freq = 50 
        
            # For DE (Nanometrics) data
            ev_info_temp.user = 'kksmith7@alaska.edu' 
            ev_info_temp.password = 'wmpo3NmqTcRm' 

            # to investigate step response
            ev_info_temp.isave_raw = True
            ev_info_temp.isave_raw_processed  = True
            ev_info_temp.rotateENZ = False
            ev_info_temp.rotateRTZ = True 
            ev_info_temp.rotateUVW = True
            # filter options
            #ev_info_temp.ifFilter = True 
            ev_info_temp.ifFilter = False 
            ev_info_temp.zerophase = False    # causal
            # filter_type = 'lowpass'
            #ev_info_temp.filter_type = 'bandpass'
            #ev_info_temp.f1 = 1/100  # fmin
            #ev_info_temp.f2 = 1/10  # fmax
            #ev_info_temp.corners = 4
            # ev_info_temp.remove_response = False
            ev_info_temp.remove_response = True
            ev_info_temp.ipre_filt = 1
            ev_info_temp.demean = True
            ev_info_temp.detrend = True
            ev_info_temp.output_cap_weight_file = False
            # ev_info_temp.outformat = 'DISP'
            ev_info_temp.ifsave_sacpaz = True
            ev_info_temp.taper = 0.2

            # append getwaveform objects
            ev_info_list.append(ev_info_temp)
        
        # always return ev_info
        ev_info = ev_info_list

    # Loop over multiple event for metric study 
    if iex == 17:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        sec_tol = 180 
        events_file = "/home/ksmith/REPOSITORIES/manuscripts/kyle/papers/basinresp/data/basinamp_obspy.txt"
        eids,otimes,elons,elats,edeps,emags = reof.read_events_obspy_file(events_file)
	
        get_SR_wf = True 
        getpreEQnoise = True; 
        ev_info_list = []
        for xx in range(3,4):
            ev_info_temp = ev_info.copy()
            ev_info_temp.otime = obspy.UTCDateTime(otimes[xx])
            ev_info_temp.elat = elats[xx]
            ev_info_temp.elon = elons[xx]
            ev_info_temp.edep = edeps[xx]
            ev_info_temp.emag = emags[xx]
            ev_info_temp.eid = eids[xx]
            
            print(ev_info_temp.tafter_sec)
 
            # subset of stations
            ev_info_temp.min_dist = 0
            ev_info_temp.max_dist = 50 
            #ev_info_temp.max_dist = 150
            
            ev_info_temp.rlat = 64.7716
            ev_info_temp.rlon = -149.1465 
            ev_info_temp.rtime = ev_info_temp.otime
            saftcase = False
            print(saftcase) 
            # Iniskin 
            if abs(ev_info_temp.otime-obspy.UTCDateTime("2016-01-24T10:30:29"))<sec_tol: 
                ev_info_temp.tafter_sec = 600; saftcase = True
	    # Chile 
            if abs(ev_info_temp.otime-obspy.UTCDateTime("2015-09-16T22:54:33"))<sec_tol: 
                ev_info_temp.rtime = obspy.UTCDateTime("2015-09-16T23:09:15.000")
                ev_info_temp.tafter_sec = 200; saftcase = True  
	    # Mariana 
            if abs(ev_info_temp.otime-obspy.UTCDateTime("2016-07-29T21:18:26"))<sec_tol: 
                ev_info_temp.rtime = obspy.UTCDateTime("2016-07-29T21:28:19.000")
	    # North Korea 
            if abs(ev_info_temp.otime-obspy.UTCDateTime("2017-09-03T03:30:01"))<sec_tol: 
                ev_info_temp.rtime = obspy.UTCDateTime("2017-09-03T03:38:55")
                ev_info_temp.tafter_sec = 600; saftcase = True  
	    # Kodiak 
            if abs(ev_info_temp.otime-obspy.UTCDateTime("2018-01-23 09:31:42"))<sec_tol: 
                ev_info_temp.tafter_sec = 1000; saftcase = True  
            # Anchorage 
            if abs(ev_info_temp.otime-obspy.UTCDateTime("2018-11-30T17:29:29"))<sec_tol: 
                ev_info_temp.tafter_sec = 600; saftcase = True 
 
            ev_info_temp.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
            ev_info_temp.channel = 'BH?,HH?'
            ev_info_temp.resample_freq = 50        
            ev_info_temp.scale_factor = 1 
            
            if get_SR_wf:
                ev_info_temp.tbefore_sec = 0 
                if not saftcase:
                    TS_length = 500
                    print("first statement",saftcase) 
                else: 
                    TS_length = ev_info_temp.tafter_sec 
                    print("second statement",saftcase) 
                
                if getpreEQnoise:
                    ev_info_temp.tafter_sec = 0
                    ev_info_temp.tbefore_sec = TS_length
 
            # EQ metrics WFs 
            else: 
                ev_info_temp.tbefore_sec = 100
                if not saftcase:
                    ev_info_temp.tafter_sec = 400
	    
                # ev_info_temp.scale_factor = 100 
            print("before sec = ", ev_info_temp.tbefore_sec)
            print("after sec = ", ev_info_temp.tafter_sec)

            # append getwaveform objects
            ev_info_list.append(ev_info_temp)
        
        # always return ev_info
        ev_info = ev_info_list
    return(ev_info)

#=================================================================================
# END EXAMPLES
#=================================================================================
