# To run this script:
# > python run_getwaveform.py
import obspy
import getwaveform_iris
import copy
from obspy.clients.fdsn import Client

c = Client("IRIS")

otime = obspy.UTCDateTime("2016-01-24T10:30:29")

cat = c.get_events(starttime = otime-10, endtime = otime+10)
network='AK,AT,YV,PS,AV,IU,II,XZ,XM'
getwaveform_iris.run_get_waveform(ev=cat[0], min_dist=0, max_dist=500, 
        before = 20, after = 600,
        network=network, channel='BH*', samp_freq=20.0, rotate=True,
        output_cap_weight_file=True, remove_response=True,
        detrend=True, demean=True, output_event_info=True,
                     pre_filt=(0.005, 0.006, 5.0, 10.0)
                     )
