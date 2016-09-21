# To run this script:
# > python getwaveform_silwal2016.py
import obspy
import copy
import getwaveform_iris
from obspy.clients.fdsn import Client
from obspy.core.event import Event, Origin, Magnitude
import numpy as np
import util_helpers as uh

# Catalog file
filename = '/home/carltape/REPOSITORIES/manuscripts/vipul/papers/2014mt/data/SCAK_mech.txt'

# Select events
iex = [80]               # default = 8 (for example event = 2009-04-07T20:12:55.351000Z)
# Line index for 21 events
#iex = [1, 3, 6, 7, 12, 29, 30, 48, 51, 66, 69, 71, 73, 77, 80, 81, 83, 88, 93, 96, 101];
header_lines = 21        # number of header lines to skip
line_indx = np.array(iex) + header_lines

#--------------------------------------------------------------------------------
# DEFAULT SETTINGS (see getwaveform_iris.py)
# Pre-processing (manily for CAP)
rotate = True
output_cap_weight_file = True
remove_response = True
detrend = True
demean = True
output_event_info = True
# pre-filter for deconvolution
# https://ds.iris.edu/files/sac-manual/commands/transfer.html
# fmaxc should be based on sampling rate (desired channels)
# fminc should be based on the request length of the time series
fminc = 1/200
fmaxc = 10
pre_filt=(0.5*fminc, fminc, fmaxc, 2.0*fmaxc)    
#pre_filt=(0.005, 0.006, 10.0, 15.0) # BH default
# for CAP all waveforms need to have the same sample rate
resample_freq = 20.0         # =0 for no resampling
scale_factor = 0             # 10**2 to convert m/s to cm/s
# event parameters
sec_before_after_event = 10  # time window to search for a target event in a catalog

# Input parameter for extracting waveforms (common for all 21 events)
min_dist = 0 
max_dist = 500
tbefore_sec = 100
tafter_sec = 300
network = 'AK,AT,AV,CN,II,IU,XM,XV,XZ,YV'  # note: cannot use '*' because of IM
station = '*'
channel = 'BH*'
use_catalog = 0                            # To get (lat,lon, etime, dep, mag) from some catalog = 1 OR use defined = 0 (see iex=9)

#--------------------------------------------------------------------------------
# Read catalog file
for indx in line_indx:
    f = open(filename,'r')
    line = f.readlines()[indx]
    line_elements = line.split()
    
    # Get event info (otime, lat, lon, dep, mag) from the catalog
    eid = line_elements[-1]   # Fix to take microseconds
    otime = uh.eid2otime(eid)
    elat = line_elements[7]
    elon = line_elements[6]
    edep = float(line_elements[8]) * 1000.0 # in meters
    emag = line_elements[16] 
    #for line in lines:
    print(otime, elat, elon, edep, emag)
    
    # Create event object
    client = Client("IRIS")
    ev = Event()
    org = Origin()
    org.latitude = elat
    org.longitude = elon
    org.depth = edep
    org.time = otime
    mag = Magnitude()
    mag.mag = emag
    mag.magnitude_type = "Mw"
    ev.origins.append(org)
    ev.magnitudes.append(mag)

    # Extract waveforms 
    getwaveform_iris.run_get_waveform(c = client, event = ev, 
                                  min_dist = min_dist, max_dist = max_dist, 
                                  before = tbefore_sec, after = tafter_sec, 
                                  network = network, station = station, channel = channel, 
                                  resample_freq = resample_freq, ifrotate = rotate,
                                  ifCapInp = output_cap_weight_file, 
                                  ifRemoveResponse = remove_response,
                                  ifDetrend = detrend, ifDemean = demean, 
                                  ifEvInfo = output_event_info, 
                                  scale_factor = scale_factor, pre_filt = pre_filt)
