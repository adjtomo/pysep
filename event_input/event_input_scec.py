import obspy
import read_event_obspy_file as reof
from getwaveform import *

def get_ev_info(ev_info,iex):
# ===============================================================
# python run_getwaveform.py event_input_scec 0
    if iex == 0:
        #ev_info.idb = 1           
        ev_info.overwrite_ddir = 0       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters

        # subset waveforms
        ev_info.min_dist = 0 
        ev_info.max_dist = 300           # 300 km: probably enough stations
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300

        # note: can we wildcard the network?
        # note: check IRIS gmap for additional networks
        # note: FOR SOME REASON YOU MAY NEED TO REQUEST IRIS DATA FIRST, THEN NCEDC
        # request 1: IRIS
        ev_info.client_name = 'IRIS'
        ev_info.network = 'AZ,CI,LB,TS,TA,NC,NN,SN,II,IU,US'
        # request 2: NCEDC
        #ev_info.client_name = 'NCEDC'
        #ev_info.network = 'BK'

        #ev_info.channel = 'BH?,HH?,EH?'
        ev_info.channel = 'BH?'
        ev_info.station = '*'

        ev_info.user = None
        ev_info.password = None

        # for CAP moment tensor code
        ev_info.resample_freq = 50
        ev_info.scale_factor = 100
        
        ievent = 16
        #events_file = "/home/vipul/REPOSITORIES/pysep/txtfiles/scec_sims_obspy.txt"
        events_file = "/home/carltape/PROJECTS/calif/MTs/scec_sims_obspy.txt"

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

            # append getwaveform objects
            ev_info_list.append(ev_info_temp)
        
        # always return ev_info
        ev_info = ev_info_list

    return(ev_info)
