import obspy
import read_event_obspy_file as reof
from getwaveform import *

def get_ev_info(ev_info,iex):
# ===============================================================
# SilwalTape2016 example event (Anchorage)
    if iex == 0:
        ev_info.use_catalog = 0
        ev_info.otime = obspy.UTCDateTime("2009-04-07T20:12:55.351")
        ev_info.min_dist = 0
        ev_info.max_dist = 300
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300

        #output all proccessing steps
        ev_info.ifverbose = True

        #keep stations with missing components and fill the missing component with a null trace (MPEN)
        #Be sure to set the null component to 0 in the weight file when running cap
        #ev_info.icreateNull = 1
        ev_info.icreateNull = 0

        #RAW and ENZ files can be used when checking if you are receiving all the data ($PYSEP/check_getwaveform.bash)
        ev_info.isave_raw = False
        ev_info.isave_raw_processed = False
        #ev_info.isave_raw = True
        #ev_info.isave_raw_processed = True
        ev_info.isave_ENZ = False
        #ev_info.isave_ENZ = True

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
        ev_info.resample_freq = 50
        ev_info.scale_factor = 100
        #ev_info.phase_window = False
        #-------for specfem------------
        #ev_info.tbefore_sec = 0
        #ev_info.resample_TF = False
        #ev_info.scale_factor = 1
        #ev_info.outformat = 'DISP'
        #------------------------------


    return(ev_info)
