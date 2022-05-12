import obspy
import read_event_obspy_file as reof
from getwaveform import *

#iex 0 = Cantwell event (local earthquake basically directly below the nodal transect)
#iex 1 = Anchorage aftershock (south of the nodal transect)
#iex 2 = largest teleseismic event during Denali nodal deployment.  Used for comparing amplitudes of nodes and permanent stations and for testing water level effects on waveforms when removing instrument responce.
#iex 3 = uses obspy list generated from paper_nodal_parks.m for many local events

#notes:
#the closest node - permanent station pair is node 1304 with permanent stations XV.FAPT and DE.UAF01.  They are just a few meters apart

#To get both nodal and permanent station data in the same data directory, you must run your pull request twice.
#first run set ev_info.ifph5 = False (or just comment it out) and request permanent stations, networks, and channels.
#second run set ev_info.ifph5 = True and request nodal stations, networks, and channels.
#be sure ev_info.use_catalog = 0 on the second run or you will overwrite your first pull request.



def get_ev_info(ev_info,iex):
# ===============================================================
    if iex == 0:
        #Cantwell local event
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        ev_info.ifverbose = True         # debugging output
        # subset waveforms (Run 1)
        ev_info.min_dist = 0
        ev_info.max_dist = 300
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300
        # for CAP
        ev_info.resample_TF = True
        ev_info.resample_freq = 50
        ev_info.scale_factor = 100

        # Run 1: non-nodal data
        ev_info.overwrite_ddir = 1 # delete data dir if it exists
        ev_info.ifph5 = False
        ev_info.network = 'AK,AT,AV,CN,II,IU,TA,US,XM,XV,DE,ZE'
        ev_info.channel = 'BH?,HH?'
        ev_info.station = "-UAF01"

        # Run 2: denali nodal data
        #ev_info.overwrite_ddir = 0
        #ev_info.ifph5 = True
        #ev_info.network = 'ZE'
        #ev_info.channel = '*'
        ##ev_info.station = '*' # CAUTION: Takes a LONG time!
        #ev_info.station = '1304' # note: PH5 does not understand -XXX station code

        ev_info.isave_raw = True
        ev_info.isave_raw_processed = True
        ev_info.rotateENZ = True
        ev_info.rotateRTZ = True

        #Cantwell local event
        ev_info.otime = obspy.UTCDateTime("2019-02-25T18:22:30.906")
        ev_info.elat = 62.8002
        ev_info.elon = -149.6240
        ev_info.edep = 73700
        ev_info.emag = 3.10

        #Cantwell local event #2
        #ev_info.otime = obspy.UTCDateTime("2019-03-10T05:16:05.815")
        #ev_info.elat = 62.8800
        #ev_info.elon = -150.526200000000
        #ev_info.edep = 88700
        #ev_info.emag = 3.10



    if iex == 1:
        #Anchorage aftershock

        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        ev_info.ifverbose = True         # debugging output
        # subset waveforms (Run 1)
        ev_info.min_dist = 0
        ev_info.max_dist = 300
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300
        # for CAP
        ev_info.resample_TF = True
        ev_info.resample_freq = 50
        ev_info.scale_factor = 100

        # Run 1: non-nodal data
        ev_info.overwrite_ddir = 1 # delete data dir if it exists
        ev_info.ifph5 = False
        ev_info.network = 'AK,AT,AV,CN,II,IU,TA,US,XM,XV,DE,ZE'
        ev_info.channel = 'BH?,HH?'

        # Run 2: denali nodal data
        #ev_info.overwrite_ddir = 0
        #ev_info.ifph5 = True
        #ev_info.network = 'ZE'
        #ev_info.channel = '*'
        ##ev_info.station = '*' # CAUTION: Takes a LONG time!
        #ev_info.station = '1304' # note: PH5 does not understand -XXX station code

        ev_info.isave_raw = True
        ev_info.isave_raw_processed = True
        ev_info.rotateENZ = True
        ev_info.rotateRTZ = True

        ev_info.user = 'ctape@alaska.edu'
        ev_info.password = 'dln3mjKtap3m9'

        #Anchorage aftershock
        ev_info.otime = obspy.UTCDateTime("2019-02-18 17:02:46.710")
        ev_info.elat = 61.4682
        ev_info.elon = -149.9607
        ev_info.edep = 37600
        ev_info.emag = 4.40


    if iex == 2:
        #teleseismic event
        #for testing/comparing absolute amplitudes of nodal stations vs permanent stations

        ev_info.overwrite_ddir = 0       # 0 = do not delete data dir if it exists (must be 0 on the second run)
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        ev_info.ifverbose = True         # debugging output

        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 3600 #hour of data for this teleseismic event

        ev_info.resample_TF = False #note: permanent stations and nodes have different sample rates
        #ev_info.resample_freq = 50
        #ev_info.scale_factor = 100

        # Run 1: non-nodal data
        ev_info.overwrite_ddir = 0 # delete data dir if it exists
        ev_info.ifph5 = False
        ev_info.network = 'AK,AT,AV,CN,II,IU,TA,US,XM,XV,DE,ZE'
        ev_info.channel = 'BH?,HH?'
        ev_info.network = 'AK,AT,AV,CN,II,IU,TA,US,XM,XV,DE' #for permanent stations
        ev_info.station = 'CUT,FAPT,UAF01,BWN,MCK,RND' #permanent stations nearest the nodes listed below

        # Run 2: denali nodal data
        #ev_info.ifph5 = True
        #ev_info.network = 'ZE'
        #ev_info.channel = '*'
        ##ev_info.station = '*' # CAUTION: Takes a LONG time!
        #ev_info.station = '1304' # note: PH5 does not understand -XXX station code #'1304'#colocated with FAPT and UAF01
        #ev_info.station = '1304,1305,154,1253,1199,1196,1197,1195,1156,1157,1158,1011,1010,1012'


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

        #ev_info.phase_window = True #WARNING this will cut your waveforms near your chosen phase
        #ev_info.phases = ["P","S"] # Choosing time period with respect to P & S
        #ev_info.write_sac_phase = True #note: cap reads in picks from header (if set to true, be sure you want this)

        #Set the water level high to have the nodal data waveforms match that of the permanent stations at longer periods (10-20s)
        ev_info.water_level = 6000  #default is 60

        #largest event during nodal deployment
        ev_info.otime = obspy.UTCDateTime("2019-02-22T10:17:28.000")
        ev_info.elat = -2.26
        ev_info.elon = -77.09
        ev_info.edep = 121100
        ev_info.emag = 7.49

    if iex == 3:
        #local events from list
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
        ev_info.resample_TF = True
        ev_info.resample_freq = 50
        ev_info.scale_factor = 100

        # Run 1: non-nodal data
        ev_info.overwrite_ddir = 0       # delete data dir if it exists
        ev_info.ifph5 = False
        ev_info.network = 'AK,AT,AV,CN,II,IU,TA,US,XM,XV,DE,ZE'
        ev_info.channel = 'BH?,HH?'

        # Run 2: denali nodal data
        #ev_info.ifph5 = True
        #ev_info.network = 'ZE'
        #ev_info.channel = '*'
        ##ev_info.station = '*' # CAUTION: Takes a LONG time!
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

    return(ev_info)
