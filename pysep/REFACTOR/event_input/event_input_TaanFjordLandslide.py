import obspy
import read_event_obspy_file as reof
from getwaveform import *

def get_ev_info(ev_info,iex):
# ===============================================================
# Taan Fjord landslide data, emulating data from Gualtieri & Ekström 2018 (GJI)
# python run_getwaveform.py event_input_TaanFjordLandslide 0
    if iex == 0:
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        ev_info.otime = obspy.UTCDateTime("2015-10-18T05:18:31.000")
        ev_info.min_dist = 0
        ev_info.max_dist = 50
        ev_info.tbefore_sec = 91
        ev_info.tafter_sec = 509

        # ENZ files can be used to compare data from Gualtieri & Ekström 2018 Broad-band seismic analysis and modeling of the 2015 Taan Fjord,Alaska landslide using Instaseis
        ev_info.isave_raw = False
        ev_info.isave_raw_processed = False
        ev_info.isave_ENZ = True

        ev_info.network = 'AK'
        ev_info.channel = 'LH?'
        ev_info.use_catalog = 0
        ev_info.elat = 60.175
        ev_info.elon = -141.187
        ev_info.edep = 0.00
        ev_info.emag = 4.9
        ev_info.resample_TF = True
        ev_info.resample_freq = 50
        ev_info.scale_factor = 1
        ev_info.filter_type = 'highpass'
        ev_info.ipre_filt = 2
        ev_info.pre_filt = (0.004, 0.005, 10.0, 15.0)
        ev_info.water_level = 10000
        ev_info.detrend = True               # detrend waveforms
        ev_info.demean = True                # demean waveforms
        ev_info.taper = 0.05

    return(ev_info)
#==============================================================================
