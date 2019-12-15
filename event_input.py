import obspy
import read_event_obspy_file as reof
from getwaveform import *

def get_ev_info(ev_info,iex):
# ===============================================================
# SilwalTape2016 example event (Anchorage) -- python run_getwaveform.py event_input 0
    if iex == 0:
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        ev_info.otime = obspy.UTCDateTime("2009-04-07T20:12:55.351")
        ev_info.min_dist = 0
        ev_info.max_dist = 300
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300

        # keep stations with missing components and fill the missing component with a null trace (MPEN)
        #ev_info.icreateNull = 1

        # RAW and ENZ files can be used when checking if you are receiving all possible data (example station: SOLD)
        ev_info.isave_raw = False
        ev_info.isave_raw_processed = False
        ev_info.isave_ENZ = False

        #ev_info.min_lat = 59
        #ev_info.max_lat = 62
        #ev_info.min_lon = -152
        #ev_info.max_lon = -147

        # default list of Alaska networks
        # note 1: cannot use '*' because of IM
        # note 2: may want to exclude the mid-band AV network
        # note 3: these are temporary:
        # XE BEAAR 1999
        # XR ARCTIC 2004
        # XZ STEEP 2005
        # YV MOOS 2006
        # XV FLATS 2014
        # ZE SALMON 2015
        # XG WVF 2016
        # [7C MMEP 2015]
        # TA
        #ev_info.network = 'AK,AT,AV,CN,II,IU,US,XM,XV,XZ,YV'
        #ev_info.network = 'AK' # for testing
        ev_info.network = 'AK,AT,AV,CN,II,IU,US,XM,TA,XE,XR,XZ,YV,XV,ZE,XG'

        ev_info.channel = 'BH?'
        ev_info.use_catalog = 0
        ev_info.elat = 61.45420
        ev_info.elon = -149.7428
        ev_info.edep = 33033.60
        # ev_info.rlat = 61.45420
        # ev_info.rlon = -149.7428
        # ev_info.rtime = obspy.UTCDateTime("2009-04-07T20:12:55.351")
        ev_info.emag = 4.6
        # scaling and resampling needed for CAP
        ev_info.resample_TF = True
        ev_info.resample_freq = 50
        ev_info.scale_factor = 100
        #ev_info.phase_window = False
        #-------for specfem------------
        #ev_info.tbefore_sec = 0
        #ev_info.resample_TF = False
        #ev_info.scale_factor = 1
        #ev_info.outformat = 'DISP'
        #------------------------------

# Iniskin earthquake
    if iex == 1:
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
        #ev_info.resample_freq = 0        # no resampling
        #ev_info.resample_TF = False
        #ev_info.scale_factor = 1         # no scale factor
        
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

        # For SALMON data
        #ev_info.user = None
        #ev_info.password = None
            
# MFFZ earthquakes for investigating the step response
# copied from event_input_flats
    if iex == 2:
        ev_info.use_catalog = 0
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2014-08-31T03:06:57.111")
        ev_info.elat = 65.1526
        ev_info.elon = -149.0398
        ev_info.edep = 16614.7
        ev_info.emag = 5.20
        # -------------------------------------------------
        
        ev_info.min_dist = 0
        ev_info.max_dist = 300
        ev_info.tbefore_sec = 200
        ev_info.tafter_sec = 600
        ev_info.network = 'AV,CN,AT,TA,AK,XV,II,IU,US'
        # ev_info.network = 'XV,AK'
        ev_info.channel = 'BH?,HH?'
        #ev_info.scale_factor = 1
        #ev_info.resample_freq = 0
        #ev_info.resample_TF = False
        # for CAP
        # ev_info.resample_TF = True
        # ev_info.scale_factor = 100
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
        
# Loop over multiple event       
    if iex == 3:
        ev_info.use_catalog = 0
        # text file of source parameters
        ievent = 1
        events_file = "/home/ksmith/REPOSITORIES/manuscripts/kyle/papers/nenanabasin/data/nenanabasin_obspy.txt"
        eids,otimes,elons,elats,edeps,emags = reof.read_events_obspy_file(events_file)

        ev_info_list = []
        #for xx in range(len(eids)):
        #for xx in range(1,3):
        for xx in range(ievent-1,ievent):
            ev_info_temp = ev_info.copy()
            ev_info_temp.otime = obspy.UTCDateTime(otimes[xx])
            ev_info_temp.elat = elats[xx]
            ev_info_temp.elon = elons[xx]
            ev_info_temp.edep = edeps[xx]
            ev_info_temp.emag = emags[xx]
            ev_info_temp.eid = eids[xx]
            
            # subset of stations
            ev_info_temp.min_dist = 0
            ev_info_temp.max_dist = 200
            ev_info_temp.tbefore_sec = 100
            ev_info_temp.tafter_sec = 500
            ev_info_temp.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
            #ev_info_temp.network = 'AK,AT,AV,CN,II,IU,US,XM,TA,XE,XR,XZ,YV,XV,ZE,XG'
            ev_info_temp.channel = 'BH?,HH?'
            #ev_info_temp.resample_freq = 40
            #ev_info_temp.scale_factor = 100
            #ev_info_temp.resample_TF = False

            # append getwaveform objects
            ev_info_list.append(ev_info_temp)
        
        # always return ev_info
        ev_info = ev_info_list

# Alaska data from far away event
    if iex == 4:
        ev_info.use_catalog = 0
        # Mariana Event observed by FLATS
        ev_info.otime = obspy.UTCDateTime("2016-07-29T21:18:26.000")
        ev_info.elat = 18.5439
        ev_info.elon = 145.541
        ev_info.edep = 207620.0
        ev_info.emag = 7.7
        ev_info.rtime = obspy.UTCDateTime("2016-07-29T21:28:19.000")
        # center at F3TN
        ev_info.rlat =  64.7716
        ev_info.rlon = -149.1465
        ev_info.min_dist = 0
        ev_info.max_dist = 150
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 400
        ev_info.network = 'AK,AT,II,IU,US,XM,XV,XZ,TA' # no CN,AV,YV,ZE
        ev_info.channel = 'BH?,HH?'
 
# doublet event near Manley  
    if iex == 5:
        ev_info.rlat = 64.7716  
        ev_info.rlon = -149.1465 
        ev_info.scale_factor = 1
        ev_info.min_dist = 0
        ev_info.max_dist = 150
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 200 
        ev_info.network = 'XV,TA'
        ev_info.channel = 'BH?,HH?'
        ev_info.use_catalog = 0
        ev_info.station = 'I23K,F3TN'
        # AEC source parameters
        ev_info.otime = obspy.UTCDateTime("2015-11-06T01:20:12.712") 
        ev_info.elat = 64.7552
        ev_info.elon = -151.3103
        ev_info.edep = 1502.1
        ev_info.emag = 3.35
        ev_info.rtime = ev_info.otime

# Illinois main event
    if iex == 6:
        ev_info.use_catalog=0
        ev_info.otime = obspy.UTCDateTime("2008-04-18T09:36:59.110")
        ev_info.elon = -87.886 
        ev_info.elat = 38.452 
        ev_info.edep = 14.3
        ev_info.emag = 5.2
        ev_info.min_dist = 0 
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300
        ev_info.network = 'IU,NM'
        ev_info.station = '*'
        ev_info.channel = 'LH?,BH?'

# Uturuncu main event, AlvizuriTape2016
    if iex == 7:
        ev_info.use_catalog=0
        ev_info.otime = obspy.UTCDateTime("2010-05-16T06:34:54.464")
        ev_info.elon = -67.1856
        ev_info.elat = -22.2600
        ev_info.edep = -0.6
        ev_info.emag = 2.8
        ev_info.min_dist = 0 
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300
        ev_info.network = 'XP'
        ev_info.station = '*'
        ev_info.channel = 'HH?'

# Manley earthquake (Kyle)
    if iex == 8:
        ev_info.use_catalog=0
        ev_info.otime = obspy.UTCDateTime("2016-05-18T03:25:48.320")
        ev_info.elat = 65.2466
        ev_info.elon = -151.0651
        ev_info.edep = 15156.2
        ev_info.emag = 4.2
        ev_info.rtime = ev_info.otime
        ev_info.max_dist = 400
        ev_info.tbefore_sec = 50
        ev_info.tafter_sec = 600
        ev_info.network = 'AK,AT,AV,CN,II,IU,US,XM,TA,XE,XR,XZ,YV,XV,ZE,XG'
        ev_info.channel = 'BH?,HH?'
        # for CAP
        ev_info.resample_TF = True
        ev_info.resample_freq = 50
        ev_info.scale_factor = 100

# same as iex=11 but for the IRIS database
# GOAL: For LLNL events, we do NOT want to use the IRIS source parameters:
#       origin time, hypocenter, magnitude.
# --> THIS EXAMPLE IS NOT CURRENTLY WORKING WITH sln01
# ERROR: FileNotFoundError: [Errno 2] No such file or directory: '19910914190000000/19910914190000000_ev_info.obj'
    if iex == 9:
        ev_info.overwrite_ddir = 0    # do NOT overwrite the existing directory
        ev_info.client_name = 'IRIS'
        #ev_info.client_name = 'NCEDC'
        ev_info.use_catalog = 1
        #ev_info.idb = 1            # IRIS database
        # ev_info.resample_freq = 0  # no resampling
        # ev_info.otime = obspy.UTCDateTime("1991-09-14T19:00:00.000Z")   # Hoya actual
        ev_info.otime = obspy.UTCDateTime("1991-09-14T19:00:08.031Z")   # Hoya target
        ev_info.min_dist = 0
        ev_info.max_dist = 1200
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        # needs to be run TWICE to get BK stations and IRIS stations
        #ev_info.network = 'BK'        # BK will go to NCEDC
        ev_info.network = '*'         # * will get all at IRIS DMC
        ev_info.channel = 'BH?,LH?' 
        ev_info.overwrite_ddir = 0
        ev_info.ifsave_stationxml = False
    
    if iex == 10:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        sec_tol = 180 
        events_file = "/home/ksmith/REPOSITORIES/manuscripts/kyle/papers/nenanabasin/data/nenanabasin_obspy.txt"
        eids,otimes,elons,elats,edeps,emags = reof.read_events_obspy_file(events_file)
	
        ev_info_list = []
        for xx in range(0,1):
            ev_info_temp = ev_info.copy()
            ev_info_temp.otime = obspy.UTCDateTime(otimes[xx])
            ev_info_temp.elat = elats[xx]
            ev_info_temp.elon = elons[xx]
            ev_info_temp.edep = edeps[xx]
            ev_info_temp.emag = emags[xx]
            ev_info_temp.eid = eids[xx]
            # subset of stations
            ev_info_temp.min_dist = 0
            #ev_info_temp.max_dist = 50 
            ev_info_temp.max_dist = 150
            ev_info_temp.rlat = 64.7716
            ev_info_temp.rlon = -149.1465 
            ev_info_temp.rtime = ev_info_temp.otime
            saftcase = False
            ev_info_temp.phases = ["P","S"] 
            ev_info_temp.write_sac_phase = True # put phase information in sac files
            ev_info_temp.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
            ev_info_temp.channel = 'BH?,HH?'
            ev_info_temp.resample_TF = True # to be consistent with past past waveforms
            ev_info_temp.resample_freq = 50        
            ev_info_temp.scale_factor = 1 
            ev_info_temp.tbefore_sec = 150
            ev_info_temp.tafter_sec = 1000
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
