import obspy
import read_event_obspy_file as reof
from getwaveform import *

def get_ev_info(ev_info,iex):
# ===============================================================

# python run_getwaveform.py event_input_museum 0

    if iex == 0:
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
        ev_info.max_dist = 200
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = '*'
        ev_info.channel = 'BN?,HN?,EN?'
        ev_info.isave_raw = True
        ev_info.outformat = 'ACC'        # save as acceleration

        ev_info.rotateRTZ = False
        ev_info.isave_ENZ = True
        ev_info.isave_raw = True
        ev_info.resample_freq = 0        # no resampling
        ev_info.resample_TF = False
        ev_info.scale_factor = 1         # no scale factor

    if iex == 1:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        
        ievent = 8
        events_file = "/home/carltape/REPOSITORIES/pysep/akstrong_obspy.txt"
        eids,otimes,elons,elats,edeps,emags = reof.read_events_obspy_file(events_file)

        ev_info_list = []
        for xx in range(ievent-1,ievent):
            ev_info_temp = ev_info.copy()
            ev_info_temp.otime = obspy.UTCDateTime(otimes[xx])
            ev_info_temp.elat = elats[xx]
            ev_info_temp.elon = elons[xx]
            ev_info_temp.edep = edeps[xx]
            ev_info_temp.emag = emags[xx]
            ev_info_temp.eid = eids[xx]
            
            # subset of stations
            ev_info_temp.min_dist = 300
            ev_info_temp.max_dist = 600
            ev_info_temp.tbefore_sec = 100
            ev_info_temp.tafter_sec = 600
            ev_info_temp.network = '*'

            #ev_info_temp.channel = 'BN?,HN?'
            #ev_info_temp.outformat = 'ACC'       # save as ACC, VEL, or DISP

            ev_info_temp.channel = 'BH?,HH?'
            ev_info_temp.outformat = 'VEL'       # save as ACC, VEL, or DISP

            ev_info_temp.rotateRTZ = False
            ev_info_temp.isave_ENZ = True
            ev_info_temp.isave_raw = True
            ev_info_temp.resample_freq = 0        # no resampling
            ev_info_temp.resample_TF = False
            ev_info_temp.scale_factor = 1         # no scale factor

            # append getwaveform objects
            ev_info_list.append(ev_info_temp)
        
        # always return ev_info
        ev_info = ev_info_list

    if iex == 2:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2012-04-11T08:38:37.300") 
        ev_info.elat = 2.24
        ev_info.elon = 92.78
        ev_info.edep = 40.03
        ev_info.emag = 8.6
        # subset of stations
        ev_info.min_lat = 64.8
        ev_info.max_lat = 65.0
        ev_info.min_lon = -148.8
        ev_info.max_lon = -146.1
        ev_info.min_lat = 64.9
        ev_info.max_lat = 65.0
        ev_info.min_lon = -148.3
        ev_info.max_lon = -148.2

        ev_info.tbefore_sec = 0
        ev_info.tafter_sec = 6200
        ev_info.network = '*'
        ev_info.channel = 'BH?,HH?'
        ev_info.outformat = 'DISP'        # DISP or VEL

        ev_info.ifverbose = True      # log file

        # default pre-filter
        #  (0.0003225816857473734, 0.0006451633714947468, 12.5, 25.0)

        ev_info.ipre_filt = 2
        ev_info.f1 = 1/1000
        ev_info.f2 = 10.
        ev_info.f0 = 0.5*ev_info.f1
        ev_info.f3 = 2.0*ev_info.f2
        ev_info.pre_filt=(ev_info.f0, ev_info.f1, ev_info.f2, ev_info.f3) 

        ev_info.rotateRTZ = False
        ev_info.isave_ENZ = True
        ev_info.isave_raw = True
        ev_info.resample_freq = 0        # no resampling
        ev_info.resample_TF = False
        ev_info.scale_factor = 1         # no scale factor

    return(ev_info)
#=================================================================================
# END EXAMPLES
#=================================================================================
