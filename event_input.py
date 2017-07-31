import obspy
import read_event_obspy_file as reof
from getwaveform import *

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
        
# SilwalTape2016 example event (Anchorage)
    if iex == 1:
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

        # For SALMON data
        ev_info.user = None
        ev_info.password = None

# Cook Inlet earthquake
    if iex == 3:
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
# copied from event_input_flats
    if iex == 4:
        ev_info.idb = 1
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
        
# Loop over multiple event       
    if iex == 5:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        
        events_file = "/home/ksmith/REPOSITORIES/manuscripts/kyle/papers/basinamp/data/basinamp_obspy.txt"
        eids,otimes,elons,elats,edeps,emags = reof.read_events_obspy_file(events_file)

        ev_info_list = []
        for xx in range(8,9):
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
            ev_info_temp.channel = 'BH?,HH?'
            ev_info_temp.resample_freq = 40        
            ev_info_temp.scale_factor = 100 
            ev_info_temp.resample_TF = False

            # append getwaveform objects
            ev_info_list.append(ev_info_temp)
        
        # always return ev_info
        ev_info = ev_info_list

# nuclear event: LLNL (see also iex = 7)
# GOAL: To find events in the LLNL database based on a target origin time, rather than an eid.
#       (The reference time (NZYEAR, etc) should then be assigned as the actual origin time,
#       not the target origin time.)
# DEBUGGING HELPER LINE:
#   saclst NPTS o b e NZHOUR NZMIN NZSEC NZMSEC f 19910914190000000/*.z
# (This will show clearly that the reference time is NOT the origin time.)
    if iex == 6:
        # to get the LLNL client, which is a private repo from Lion Krischer:
        # cd $REPOS
        # git clone https://GITHUBUSERNAME@github.com/krischer/llnl_db_client.git
        # then follow instructions for install
        ev_info.idb = 3              # LLNL database
        # resample_freq = 0    # no resampling and no cutting
        
        # TARGET origin time (8.031 s from actual origin time)
        ev_info.otime = obspy.UTCDateTime("1991-09-14T19:00:08.031Z")
        # evid = 635527        # Hoya event id in LLNL database
        
        ev_info.min_dist = 0 
        ev_info.max_dist = 1200
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = '*'        # note that the client will look for BK stations in the list
        ev_info.channel = 'BH?'      # ALL channels from LLNL are returned
        # ev_info.scale_factor = 10.0**2  # original
        ev_info.scale_factor = 2e-1     # Hoya  
        ev_info.overwrite_ddir = 0

# same as iex=6 but for the IRIS database
# GOAL: For LLNL events, we do NOT want to use the IRIS source parameters:
#       origin time, hypocenter, magnitude.
    if iex == 7:
        ev_info.idb = 1            # IRIS database
        # ev_info.resample_freq = 0  # no resampling
        # ev_info.otime = obspy.UTCDateTime("1991-09-14T19:00:00.000Z")   # Hoya actual
        ev_info.otime = obspy.UTCDateTime("1991-09-14T19:00:08.031Z")   # Hoya target
        ev_info.min_dist = 0 
        ev_info.max_dist = 1200
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        # needs to be run TWICE to get BK stations and IRIS stations
        # ev_info.network = 'BK'        # BK will go to NCEDC
        ev_info.network = '*'         # * will give all at IRIS DMC
        ev_info.channel = 'BH?,LH?' 
        ev_info.overwrite_ddir = 0

# Alaska data from far away event
    if iex == 8:
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
        ev_info.scale_factor = 1
        ev_info.min_dist = 0
        ev_info.max_dist = 150
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 400
        ev_info.network = 'AK,AT,II,IU,US,XM,XV,XZ,TA' # no CN,AV,YV,ZE
        ev_info.channel = 'BH?,HH?'
    
    if iex == 9:
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
        # ev_info.resample_TF = False
        # ev_info.resample_freq = 50

    return(ev_info)
#=================================================================================
# END EXAMPLES
#=================================================================================
