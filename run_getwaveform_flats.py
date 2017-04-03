#=============================================================
# child of run_getwaveform.py except it focuses on notes_basin.tex earthquakes
# 
# This script will fetch seismic waveforms, then process them, then write sac files.
# Used heavily within the UAF seismology group.
#
# This script contains a large number of examples in two categories:
# A. examples that target a current or previous bug
# B. examples of important events for modeling efforts that others may want to reproduce
#
# In the future, we will try to automatically run these examples for each code update.
# For now, we will simply try to regularly re-run the examples.
#
# contributors: Celso Alvizuri, Lion Krischer, Vipul Silwal, Carl Tape
# 
# To run this script:
# python run_getwaveform.py
#
# TO DO
# + filetags for the case iFilter = True (like lp10, bp10_40, etc)
#
#=============================================================

import obspy
import copy
import util_helpers
import shutil   # only used for deleting data directory
import os
import sys
import getwaveform

# DEFAULT SETTINGS (see getwaveform_iris.py)
idb = 1    # default: =1-IRIS; =2-AEC; =3-LLNL
# Pre-processing (manily for CAP)
rotateRTZ = True
rotateUVW = False   # This option works only if 'rotateRTZ = True'
output_cap_weight_file = True
detrend = True
demean = True
taper = False
output_event_info = True
outformat = 'VEL'            # Intrument response removed waveforms could be saved as 'VEL' 'DISP' 'ACC'
ifsave_sacpaz = False        # save sac pole zero (needed as input for MouseTrap module)

# for CAP all waveforms need to have the same sample rate
resample_freq = 50.0         #
         
# event parameters

sec_before_after_event = 10  # time window to search for a target event in a catalog
min_dist = 0 
# station parameters
network = '*'                # all networks
station = '*,-PURD,-NV33,-GPO'  # all stations
overwrite_ddir = 1           # 1 = delete data directory if it already exists
icreateNull = 1              # create Null traces so that rotation can work (obsby stream.rotate require 3 traces)

#------Filter--------------
# for 'bandpass' both f1 and f2 are used
# for 'lowpass' only f2 is used
# for 'highpass' only f1 is used
#
# EXAMPLES
#                                   ifFilter  zerophase  remove_response  ipre_filt
# A. CAP-ready waveforms [DEFAULT]: False     NA         True             1
# B. plot-ready waveforms:          True      True       True             2
# C. plot-ready, causal waveforms:  True      False      True             0
# D. possible sensor issues:        True      False      False            NA
#
# Perhaps you want to set ipre_filt = 0 to prevent extra filtering
ifFilter = False 
filter_type = 'bandpass'
# f1 should consider the requested length of the time series
# f2 should consider the sampling rate for the desired channels
f1 = 1/40 # fmin - highpass will keep frequencies larger than fmin
f2 = 1/10 # fmax - lowpass will keep frequencies lower than fmax
zerophase = True             # = False (causal), = True (acausal)
# 4 pole filter is more sharper at the edges than 2 pole
corners = 4                  # Is corner in Obspy same as Pole in SAC?

# pre-filter for deconvolution
# https://ds.iris.edu/files/sac-manual/commands/transfer.html
# Pre-filter will not be applied if remove_response = False 
remove_response = True
iplot_response = False
ifplot_spectrogram = False
ipre_filt = 1                # =0 No pre_filter
                             # =1 default pre_filter (see getwaveform_iris.py)
                             # =2 user-defined pre_filter (use this if you are using bandpass filter)
# For tapering down the pre-filter
f0 = 0.5*f1
f3 = 2.0*f2
pre_filt=(f0, f1, f2, f3)    # applies for ipre_filt = 2 only
#pre_filt = (0.005, 0.006, 10.0, 15.0) # BH default
idb = 1
overwrite_ddir = 1       # delete data dir if it exists
use_catalog = 0          #  # use an existing catalog (=1) or specify your own source parameters (=0)
# subset of stations
min_dist = 0
max_dist = 100
tbefore_sec = 100
tafter_sec = 400
network = 'AK,AT,II,IU,US,XM,XV,XZ,TA' # no CN,AV,YV,ZE
channel = 'BH?,HH?'
resample_freq = 0        #  =0 for no resampling
scale_factor = 1         # no scale factor, for CAP use 10**2  (to convert m/s to cm/s)
remove_response = True

# username and password for accessing embargoed data from IRIS
# Register here: http://ds.iris.edu/ds/nodes/dmc/forms/restricted-data-registration/
# Run example iex = 4 to check
user = ''
password = ''

for iex in range(14,15):
    # dummy values
    dummyval = -9999
    rlat = dummyval
    rlon = dummyval
    rtime = dummyval


    print("my iex is %d" %(iex))
    # Kyle Nenana basin earthquakes for basin amplification
    if iex == 1:
        # AEC source parameters
        otime = obspy.UTCDateTime("2015-11-20T10:53:48.168") 
        elat = 64.6210
        elon = -149.4024
        edep = 17113.4
        emag = 2.67

    if iex == 2:
        # AEC source parameters
        otime = obspy.UTCDateTime("2015-10-22T13:16:15.794") 
        elat = 64.7334
        elon = -149.0388
        edep = 18830.2
        emag = 2.74   

    if iex == 3:
        # AEC source parameters
        otime = obspy.UTCDateTime("2014-12-13T15:47:31.423") 
        elat = 64.4325
        elon = -149.3840
        edep = 12431.1
        emag = 3.25    

    if iex == 4:
        # AEC source parameters
        otime = obspy.UTCDateTime("2015-10-31T02:56:35.572") 
        elat = 64.4285
        elon = -149.6969
        edep = 23852.1
        emag = 3.47

    if iex == 5:
        # AEC source parameters
        otime = obspy.UTCDateTime("2015-11-06T01:20:12.712") 
        elat = 64.7552
        elon = -151.3103
        edep = 1502.1
        emag = 3.35

    if iex == 6:
        # AEC source parameters
        otime = obspy.UTCDateTime("2014-10-21T00:36:58.333") 
        elat = 65.1489
        elon = -149.0413
        edep = 13134.8
        emag = 4.90

    if iex == 7:
        # AEC source parameters
        otime = obspy.UTCDateTime("2014-10-23T16:30:23.968") 
        elat = 65.1644
        elon = -149.0523
        edep = 200665
        emag = 5.00

    if iex == 8:
        # AEC source parameters
        otime = obspy.UTCDateTime("2016-01-14T19:04:10.727") 
        elat = 64.6827
        elon = -149.2479
        edep = 22663.7
        emag = 3.80

    if iex == 9:
        # AEC source parameters
        otime = obspy.UTCDateTime("2015-09-28T11:50:12.949") 
        elat = 64.7148
        elon = -148.9769
        edep = 15112.7
        emag = 2.91

    if iex == 10:
        # Big Minto Event
        # AEC source parameters
        otime = obspy.UTCDateTime("2016-11-06T09:29:10.579") 
        elat = 64.1639
        elon = -150.0626
        edep = 23190.0
        emag = 4.00
        max_dist = 150

    if iex == 11:
        # Big Minto Event
        # AEC source parameters
        otime = obspy.UTCDateTime("2016-12-08T10:18:13.868") 
        elat = 64.1937
        elon = -150.0376
        edep = 24522.1
        emag = 4.30
        max_dist = 150

    if iex == 12:
        # Iniskin Event
        otime = obspy.UTCDateTime("2016-01-24T10:30:29.557") 
        elat = 59.6204
        elon = -153.3392
        edep = 125645.3
        emag = 7.10
        rlat = 64.7716
        rlon = -149.1465
        rtime = otime
        tafter_sec = 600

    if iex == 13:
        # Chile Event
        otime = obspy.UTCDateTime("2015-09-16T22:54:33.000")
        elat = -31.5695
        elon = -71.6543
        edep = 22400.0
        emag = 8.30
        rlat = 64.7716
        rlon = -149.1465
        rtime = obspy.UTCDateTime("2015-09-16T23:09:15.000")
        tafter_sec = 200
    
    if iex == 14:
        # Mariana Event
        otime = obspy.UTCDateTime("2016-07-29T21:18:26.000")
        elat = 18.5439
        elon = 145.541
        edep = 207620.0
        emag = 7.7
        rlat = 64.7716
        rlon = -149.1465
        rtime = obspy.UTCDateTime("2016-07-29T21:28:29.000")
        tafter_sec = 250

    if rlat == dummyval:
        # By default this should be the event time and location unless we want to grab stations centered at another location
        rlat = elat
        rlon = elon
        rtime = otime


# fetch and process waveforms
# IRIS
    if idb == 1:
    # import functions to access waveforms
    #import getwaveform_iris
        from obspy.clients.fdsn import Client
        from obspy.core.event import Event, Origin, Magnitude
        if not user and not password:
            client = Client("IRIS")
        else:
            client = Client("IRIS",user=user,password=password)
    # will only work for events in the 'IRIS' catalog
    # (future: for Alaska events, read the AEC catalog)
        if use_catalog==1:
            print("WARNING using event data from the IRIS catalog")
            cat = client.get_events(starttime = otime - sec_before_after_event,\
                                       endtime = otime + sec_before_after_event)
            ev = cat[0]
        else:
            print("WARNING using event data from user-defined catalog")
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

            ref_time_place = Event()
            ref_org = Origin()
            ref_org.latitude = rlat
            ref_org.longitude = rlon
            ref_org.time = rtime
            ref_org.depth = 0 # dummy value
            ref_time_place.origins.append(ref_org)
            ref_time_place.magnitudes.append(mag) # more dummies
            
    # Delete existing data directory
    eid = util_helpers.otime2eid(ev.origins[0].time)
    ddir = './'+ eid
#if os.path.exists('RAW'):
#    print("WARNING. %s already exists. Deleting ..." % ddir)
#    shutil.rmtree('RAW')
    if overwrite_ddir and os.path.exists(ddir):
        print("WARNING. %s already exists. Deleting ..." % ddir)
        shutil.rmtree(ddir)

# Extract waveforms, IRIS
    getwaveform.run_get_waveform(c = client, event = ev, idb = idb, ref_time_place = ref_time_place, 
                                 min_dist = min_dist, max_dist = max_dist, 
                                 before = tbefore_sec, after = tafter_sec, 
                                 network = network, station = station, channel = channel, 
                                 resample_freq = resample_freq, ifrotateRTZ = rotateRTZ, ifrotateUVW = rotateUVW,
                                 ifCapInp = output_cap_weight_file, 
                                 ifRemoveResponse = remove_response,
                                 ifDetrend = detrend, ifDemean = demean, ifTaper = taper,
                                 ifEvInfo = output_event_info,
                                 scale_factor = scale_factor,
                                 ipre_filt = ipre_filt, pre_filt = pre_filt, 
                                 icreateNull=icreateNull,
                                 ifFilter = ifFilter, fmin = f1, fmax = f2, filter_type = filter_type, 
                                 zerophase = zerophase, corners = corners, 
                                 iplot_response = iplot_response, ifplot_spectrogram = ifplot_spectrogram,
                                 outformat = outformat, ifsave_sacpaz = ifsave_sacpaz)
