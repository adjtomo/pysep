import obspy
import read_event_obspy_file as reof
from getwaveform import *

def get_ev_info(ev_info,iex):
# ===============================================================
# North Korea events in Alvizuri and Tape (2018 SRL)
# python run_getwaveform.py event_input_nkfmtu 0
    if iex == 0:
        ev_info.use_catalog = 0     # do not use event catalog for source parameters
        # subset waveforms
        ev_info.min_dist = 0 
        ev_info.max_dist = 3000
        ev_info.tbefore_sec = 100   # should probably be >300, especially for closest stations
        ev_info.tafter_sec = 300    # 300 testing, 3000 final
        ev_info.resample_TF = True  # for CAP
        ev_info.resample_freq = 20  # for CAP (might want 50 for body waves)
        ev_info.scale_factor = 100  # for CAP (cm/s)

        # NOTE network PN is disabled since a station with missing deconv info crashes the scripts (when doing a full-station search for available data)
        # To see this add "print(station)" as line 238 in /home/calvizur/UTILS/anaconda3/envs/seis/lib/python3.6/site-packages/obspy/io/stationxml/core.py 
        #ev_info.channel = '*,-HH?,-BH?,-LH?'
        #ev_info.network = 'IU,-PN,-IM' # 
        #ev_info.channel = 'HH?,BH?,LH?'
        #ev_info.network = '*,-PN'

        # 2017-12-20 ISSUES when requesting all channels ('*'):
        #   No data available for request (nk06 event): IM JP KG SY
        #   div/0 error: network XG 
        #   event 2017-11-15 station MG04.XL: Could not find a valid Response Stage Type 
        #   event 2017-11-15: Exception: Can't merge traces with same ids but differing sampling rates! -- station BUS2, channel BH (S. korea)
        # ev_info.network = 'XG'
        # ev_info.station = 'STZ'
        #ev_info.channel = '*' 
        #ev_info.station = '*,-BUS2,-MG04' # event NK 2017-11-15 avoid BUS2. see exception above
        ev_info.channel = 'HH?,BH?,LH?'
        ev_info.network = '*,-PN,-XL'
        
        # 20191120 testing with python 3.7.4, obspy 1.1.0
        # 20061009013528000 -- 132 three-component seismograms (processed and rotated)
        # 20090525005443124 -- 178 Z-component seismograms
        #   rotation error: Station YM.38..BH* Rotating random orientation to NEZ. (also check YC.BEZH)
        # 20130212025751273 --  72 three-component seismograms (processed and rotated)
        # 20160106013000964 --  90 three-component seismograms (processed and rotated)
        # 20160909003001386 -- 108 three-component seismograms (processed and rotated)
        # 20170903033001760 -- 108 three-component seismograms (processed and rotated)
        # 20170903033831810 -- 105 three-component seismograms (processed and rotated)
        # 20160912113255770 -- 138 Z-component seismograms
        #   rotation error: Station MI.ANSV..BH* Rotating random orientation to NEZ. (also check IU.ULN)
        # 20171115052932820 -- 122 Z-component seismograms
        #   rotation error: Station MI.GOLF..BH* Rotating random orientation to NEZ. (also check RM.SLV)

        ievent = 6    # ievent = 6: 2017-09-03 nuke
        events_file = "./event_input/input/nkfmtu_obspy.txt"        
        eids,otimes,elons,elats,edeps,emags = reof.read_events_obspy_file(events_file)

        ev_info_list = []
        # run subset of events
        #for xx in range(len(eids)):
        for xx in range(ievent-1,ievent):
        #for xx in range(2,ievent):
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
