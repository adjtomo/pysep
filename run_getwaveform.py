# To run this script:
# > python run_getwaveform.py
import obspy
import getwaveform_iris
import copy
from obspy.clients.fdsn import Client

c = Client("IRIS")

# EXAMPLES:
iex = 2

# Iniskin earthquake
if iex == 1:
    otime = obspy.UTCDateTime("2016-01-24T10:30:29")
    min_dist = 0 
    max_dist = 500
    before = 20
    after = 600
    network = 'AK,AT,YV,PS,AV,IU,II,XZ,XM'
    channel = 'BH*'
    samp_freq = 20.0
    rotate = True
    output_cap_weight_file = True
    remove_response = True
    detrend = True
    demean = True
    output_event_info = True
    pre_filt = (0.005, 0.006, 5.0, 10.0)

# SilwalTape2016 example event (Anchorage)
if iex == 2:
    otime = obspy.UTCDateTime("2009-04-07T20:12:55")
    min_dist = 0 
    max_dist = 500
    before = 100
    after = 300
    network = 'AK,AT,YV,PS,AV,IU,II,XZ,XM,IM'
    channel = 'BH*'
    samp_freq = 20.0
    rotate = True
    output_cap_weight_file = True
    remove_response = True
    detrend = True
    demean = True
    output_event_info = True
    pre_filt = (0.005, 0.006, 5.0, 10.0)

# Example for showing the script crashes when requesting for IM network stations
if iex == 3:
    otime = obspy.UTCDateTime("2009-04-07T20:12:55")
    min_dist = 300 
    max_dist = 500
    before = 100
    after = 300
    #network = 'AK,YV,IM' # Crashes (Error because of the IM network)
    network = 'AK,YV'   # works fine 
    channel = 'BH*'
    samp_freq = 20.0
    rotate = True
    output_cap_weight_file = True
    remove_response = True
    detrend = True
    demean = True
    output_event_info = True
    pre_filt = (0.005, 0.006, 5.0, 10.0)

cat = c.get_events(starttime = otime-10, endtime = otime+10)

# Extract waveforms 
getwaveform_iris.run_get_waveform(ev = cat[0], min_dist = min_dist, max_dist = max_dist, 
        before = before, after = after, network = network, channel = channel, 
        samp_freq = samp_freq, rotate = rotate,
        output_cap_weight_file = output_cap_weight_file, remove_response = remove_response,
        detrend = detrend, demean = demean, output_event_info = output_event_info, pre_filt = pre_filt)
