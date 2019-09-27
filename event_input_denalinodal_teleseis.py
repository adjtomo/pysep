import obspy
import read_event_obspy_file as reof
from getwaveform import *

def get_ev_info(ev_info,iex):
# ===============================================================
# denali nodal
# python run_getwaveform.py event_input_denalinodal 0
    if iex == 0:
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        ev_info.ifverbose = True         # debugging output
        # subset waveforms (Run 1)
        ev_info.min_dist = 0 
        ev_info.max_dist = 300
        #ev_info.min_lat = 60
        #ev_info.max_lat = 65
        #ev_info.min_lon = -151
        #ev_info.max_lon = -148
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300
        # for CAP
        ev_info.resample_freq = 50
        ev_info.scale_factor = 100

        # Run 1: non-nodal data
        ev_info.overwrite_ddir = 0       # delete data dir if it exists
        #ev_info.ifph5 = False
        #ev_info.network = 'AK,AT,AV,CN,II,IU,TA,US,XM,XV,DE,ZE'
        #ev_info.channel = 'BH?,HH?'

        # Run 2: denali nodal data
        ev_info.ifph5 = True
        ev_info.network = 'ZE'
        ev_info.channel = '*'
        ev_info.station = '*' # CAUTION: Takes a LONG time!
        #ev_info.station = '5900' # note: PH5 does not understand -XXX station code

        ev_info.isave_raw = True 
        ev_info.isave_raw_processed = True
        ev_info.rotateENZ = True
        ev_info.rotateRTZ = True 

        ev_info.user = 'ctape@alaska.edu'
        ev_info.password = 'dln3mjKtap3m9'

        # see paper_nodal_parks.m
        ievent = 1
        events_file = "/home/carltape/REPOSITORIES/GEOTOOLS/python_util/pysep_dev/input/denali_parks_obspy.txt"
        eids,otimes,elons,elats,edeps,emags = reof.read_events_obspy_file(events_file)

        ev_info_list = []
        #for xx in range(len(eids)):
        for xx in range(ievent-1,ievent):
            ev_info_temp = ev_info.copy()
            ev_info_temp.otime = obspy.UTCDateTime(otimes[xx])
            ev_info_temp.elat = elats[xx]
            ev_info_temp.elon = elons[xx]
            ev_info_temp.edep = edeps[xx]
            ev_info_temp.emag = emags[xx]
            ev_info_temp.eid = eids[xx]

            # append getwaveform objects
            ev_info_list.append(ev_info_temp)
        
        # always return ev_info
        ev_info = ev_info_list

    if iex == 1:

        ev_info.overwrite_ddir = 0       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        ev_info.ifverbose = True         # debugging output
        # subset waveforms
        #ev_info.min_dist = 0 
        #ev_info.max_dist = 300
        #ev_info.min_lat = 60
        #ev_info.max_lat = 65
        #ev_info.min_lon = -151
        #ev_info.max_lon = -148
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 3600
        #ev_info.resample_freq = 50
        #ev_info.scale_factor = 100
        ev_info.network = 'AK,AT,TA,XV,DE'
        #ev_info.network = 'AK,AT,AV,CN,II,IU,TA,US,XM,XV,DE'
        #ev_info.network = 'AK,AT,AV,CN,II,IU,TA,US,XM,XV,DE,ZE'
        #ev_info.network = 'ZE' #nodes
        #ev_info.station = '*' # CAUTION: Takes a LONG time!
        ev_info.station = 'CUT,FAPT,UAF01,BWN,MCK,RND'
        #ev_info.station = '1304,1305,154,1253,1199,1196,1197,1195,1156,1157,1158,1011,1010,1012' #nodes for comparison with permanent stations
        #ev_info.station = '5900' # PH5 does not understand -XXX station code
        ev_info.channel = 'BH?,HH?,DH?' #for permanent station       
        #ev_info.channel = '*' #for nodes
        # Rotate 
        # ERROR (R==> DP1 90.0 0.0 DP2 90.0 0.0 DPZ 90.0 0.0)
        #       ValueError: The given directions are not linearly independent, 
        #       at least within numerical precision. Determinant of the base change matrix: 0
        ev_info.isave_raw = True 
        ev_info.isave_raw_processed = True
        ev_info.rotateENZ = True
        ev_info.rotateRTZ = True 

        #ev_info.ifph5 = True #must be true for nodes and false (commented out) for permanent stations
        ev_info.user = 'ctape@alaska.edu'
        ev_info.password = 'dln3mjKtap3m9'

        ev_info.phase_window = True
        ev_info.phases = ["P","S"] # Choosing time period with respect to P & S
        ev_info.write_sac_phase = True
       
        ev_info.water_level = 600000
 
        #filter
        ev_info.ipre_filt = 2
        ev_info.filter_type = 'bandpass'
        ev_info.f1 = 1/20
        ev_info.f2 = 1/10
        

        ev_info.use_catalog = 0
        #ev_info.otime = obspy.UTCDateTime("2016-01-02T01:00:00.000")
        #ev_info.otime = obspy.UTCDateTime("2019-03-15T05:58:36.024")
        #ev_info.elat = 63.2812
        #ev_info.elon = -151.1263
        #ev_info.edep = 9000
        #ev_info.emag = 3.9
        ##
        ev_info.otime = obspy.UTCDateTime("2019-02-22T10:17:28.000")
        ev_info.elat = -2.26
        ev_info.elon = -77.09
        ev_info.edep = 121100
        ev_info.emag = 7.49 
        #ev_info.otime = obspy.UTCDateTime("2019-03-06T00:34:26.853")
        #ev_info.elat = 62.0550
        #ev_info.elon = -148.7027
        #ev_info.edep = 39600
        #ev_info.emag = 3.30

    if iex == 2:
    #for testing/comparing absolute amplitudes of nodal stations vs permanent stations
    #must run twice in order to get nodal and permanent stations data.
    #first run set ev_info.ifph5 = True and request nodal stations
    #second run set ev_info.ifph5 = False (or just comment it out) and request permanent stations
        
        ev_info.overwrite_ddir = 0       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        ev_info.ifverbose = True         # debugging output
        
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 3600 #hour of data for this teleseismic event

        ev_info.resample_TF = True
        ev_info.resample_freq = 50 #this is the default
        ev_info.scale_factor = 100 
        ev_info.network = 'AK,AT,AV,CN,II,IU,TA,US,XM,XV,DE' #for permanent stations
        #ev_info.network = 'ZE' #nodes
        ###ev_info.station = 'CUT,FAPT,UAF01,BWN,MCK,RND' #permanent stations near nodes below
        ###ev_info.station = '1304,1305,154,1253,1199,1196,1197,1195,1156,1157,1158,1011,1010,1012'
        ev_info.station = 'FAPT,UAF01' #for comparison with node 1304 (may need password for UAF01)
        #ev_info.station = '1304'#colocated with FAPT and UAF01
        ev_info.channel = 'BH?,HH?' #for permanent station       
        #ev_info.channel = '*' #for nodes
        
        # Rotate 
        # ERROR (R==> DP1 90.0 0.0 DP2 90.0 0.0 DPZ 90.0 0.0)
        #       ValueError: The given directions are not linearly independent, 
        #       at least within numerical precision. Determinant of the base change matrix: 0
        ev_info.isave_raw = True
        ev_info.isave_raw_processed = True
        ev_info.rotateENZ = True
        ev_info.rotateRTZ = True
        
        #ev_info.ifph5 = True #must be true for nodes and false (commented out) for permanent stations
        ev_info.user = 'ctape@alaska.edu'
        ev_info.password = 'dln3mjKtap3m9'
        
        #ev_info.phase_window = True
        ev_info.phases = ["P","S"] # Choosing time period with respect to P & S
        ev_info.write_sac_phase = True
        
        #Set the water level high to have the nodal data waveforms match that of the permanent stations
        ev_info.water_level = 6000  #default is 60
        
        #filter
        ev_info.ipre_filt = 2
        ev_info.filter_type = 'bandpass'
        ev_info.f1 = 1/20 #20s
        ev_info.f2 = 1/10 #10s
        
        ev_info.otime = obspy.UTCDateTime("2019-02-22T10:17:28.000")
        ev_info.elat = -2.26
        ev_info.elon = -77.09
        ev_info.edep = 121100
        ev_info.emag = 7.49
    return(ev_info)
