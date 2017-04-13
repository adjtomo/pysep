#=============================================================
# run_getwaveform.py
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

# EXAMPLES (choose one)
iex = 4

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

# for CAP all waveforms need to have the same sample rate
resample_freq = 50.0         # =0 for no resampling
scale_factor = 10**2         # for CAP use 10**2  (to convert m/s to cm/s)
# event parameters
use_catalog = 1              # use an existing catalog (=1) or specify your own event parameters (see iex=9)
sec_before_after_event = 10  # time window to search for a target event in a catalog
min_dist = 0 
max_dist = 20000
# station parameters
network = '*'                # all networks
station = '*,-PURD,-NV33,-GPO'  # all stations
overwrite_ddir = 1           # 1 = delete data directory if it already exists
icreateNull = 1              # create Null traces so that rotation can work (obsby stream.rotate require 3 traces)

#------Filter--------------
# for 'bandpass' both f1 and f2 are used
# for 'lowpass' only f1 is used
# for 'highpass' only f2 is used
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
f1 = 1/40
f2 = 1/10
zerophase = True             # = False (causal), = True (acausal)
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
f0 = 0.5*f1
f3 = 2.0*f2
pre_filt=(f0, f1, f2, f3)    # applies for ipre_filt = 2 only
#pre_filt = (0.005, 0.006, 10.0, 15.0) # BH default

# username and password for accessing embargoed data from IRIS
# Register here: http://ds.iris.edu/ds/nodes/dmc/forms/restricted-data-registration/
# Run example iex = 4 to check
user = 'vsilwal@alaska.edu'
password = 'dmk3hjKl3Jnq'

# EXAMPLE TEMPLATE -- DO NOT USE
# This is a template for testing before creating an example.
# We can delete this if it creates more issues.
if iex == 0:
    otime = obspy.UTCDateTime("2009-04-07T20:12:55")
    min_dist = 0 
    max_dist = 500
    tbefore_sec = 100
    tafter_sec = 300
    network = 'II,IU'
    channel = '*'

#===============================================================
# Cook Inlet events
# See get_events_NEHRP.m
use_catalog = 0

if iex == 1:
    otime = util_helpers.eid2otime("20090502095052332")
    elat = 60.88
    elon = -150.93
    edep = 4340
    emag = 2.57

if iex == 2:
    otime = util_helpers.eid2otime("20100618081034529")
    elat = 61.09
    elon = -151.10
    edep = 7310
    emag = 2.64

if iex == 3:
    otime = util_helpers.eid2otime("20130326043006128")
    elat = 60.93
    elon = -150.86
    edep = 9570
    emag = 3.49

if iex == 4:
    otime = util_helpers.eid2otime("20081122053047880")
    elat = 60.99
    elon = -151.16
    edep = 99200
    emag = 2.92

if iex == 5:
    otime = util_helpers.eid2otime("20120308105743275")
    elat = 61.01
    elon = -150.91
    edep = 9990
    emag = 4.0

if iex == 6:
    otime = util_helpers.eid2otime("20120213174033674")
    elat = 60.93
    elon = -151.09
    edep = 10450
    emag = 2.51

if iex == 7:
    otime = util_helpers.eid2otime("20101013144537694")
    elat = 61.08
    elon = -150.94
    edep = 11350
    emag = 2.80

if iex == 8:
    otime = util_helpers.eid2otime("20130324152430120")
    elat = 60.92
    elon = -150.83
    edep = 11930
    emag = 2.70

if iex == 9:
    otime = util_helpers.eid2otime("20101228230830757")
    elat = 61.00
    elon = -150.94
    edep = 12710
    emag = 2.84

if iex == 10:
    otime = util_helpers.eid2otime("20080415084217360")
    elat = 60.97
    elon = -151.13
    edep = 12850
    emag = 2.85

if iex == 11:
    otime = util_helpers.eid2otime("20081006182438531")
    elat = 61.15
    elon = -150.76
    edep = 13330
    emag = 2.58

if iex == 12:
    otime = util_helpers.eid2otime("20090905015236732")
    elat = 60.94
    elon = -151.08
    edep = 14400
    emag = 2.90

if iex == 13:
    otime = util_helpers.eid2otime("20150315085611245")
    elat = 61.03
    elon = -150.79
    edep = 14560
    emag = 2.5

if iex == 14:
    otime = util_helpers.eid2otime("20150830212712960")
    elat = 61.00
    elon = -150.96
    edep = 14800
    emag = 2.56

if iex == 15:
    otime = util_helpers.eid2otime("20140203000307591")
    elat = 60.92
    elon = -151.13
    edep = 14900
    emag = 2.89

if iex == 16:
    otime = util_helpers.eid2otime("20150727022154395")
    elat = 60.98
    elon = -150.94
    edep = 16020
    emag = 3.44

if iex == 17:
    otime = util_helpers.eid2otime("20080408171630028")
    elat = 61.06
    elon = -150.85
    edep = 16.30
    emag = 2.59

# Only one good station
# Didn't perform CAP inversion
if iex == 18:
    otime = util_helpers.eid2otime("20120528010250849")
    elat = 60.91
    elon = -151.07
    edep = 17070
    emag = 2.59

if iex == 19:
    otime = util_helpers.eid2otime("20120802061138110")
    elat = 60.82
    elon = -151.02
    edep = 18280
    emag = 2.93

if iex == 20:
    otime = util_helpers.eid2otime("20151124212557376")
    elat = 60.94
    elon = -150.82
    edep = 19940
    emag = 2.68

if iex == 21:
    otime = util_helpers.eid2otime("20141228170032208")
    elat = 60.95
    elon = -150.87
    edep = 19990
    emag = 3.06

if iex == 22:
    otime = util_helpers.eid2otime("20141115030100885")
    elat = 60.76
    elon = -151.07
    edep = 20030
    emag = 2.52

if iex == 23:
    otime = util_helpers.eid2otime("20141211004839285")
    elat = 60.74
    elon = -151.03
    edep = 20050
    emag = 2.98

#=======================================================
# Beluga events

if iex == 24:
    otime = util_helpers.eid2otime("20080205035142446")
    elat = 61.55
    elon = -151.28
    edep = 1960
    emag = 2.62

if iex == 25:
    otime = util_helpers.eid2otime("20100328160536582")
    elat = 61.69
    elon = -151.34
    edep = 5190
    emag = 3.19

if iex == 26:
    otime = util_helpers.eid2otime("20090516015104343")
    elat = 61.66
    elon = -151.25
    edep = 6890
    emag = 2.90

if iex == 27:
    otime = util_helpers.eid2otime("20160418180212156")
    elat = 61.61
    elon = -151.22
    edep = 7490
    emag = 2.60

if iex == 28:
    otime = util_helpers.eid2otime("20120629110739385")
    elat = 61.62
    elon = -151.30
    edep = 8530
    emag = 2.69

if iex == 29:
    otime = util_helpers.eid2otime("20140124120703813")
    elat = 61.65
    elon = -151.26
    edep = 8700
    emag = 3.56

if iex == 30:
    otime = util_helpers.eid2otime("20080126042942584")
    elat = 61.56
    elon = -151.23
    edep = 11270
    emag = 3.15

if iex == 31:
    otime = util_helpers.eid2otime("20140714060410234")
    elat = 61.59
    elon = -151.29
    edep = 11500
    emag = 3.01

if iex == 32:
    otime = util_helpers.eid2otime("20120306061258556")
    elat = 61.54 
    elon = -151.25
    edep = 12420
    emag = 2.75

#========================================================
# North susitna events

if iex == 33:
    otime = util_helpers.eid2otime("20160204082713319")
    elat = 62.01
    elon = -150.72
    edep = 140
    emag = 2.53

if iex == 34:
    otime = util_helpers.eid2otime("20160118172116475")
    elat = 62.06
    elon = -150.66
    edep = 1
    emag = 2.73

if iex == 35:
    otime = util_helpers.eid2otime("20080418041458669")
    elat = 62.05
    elon = -150.50
    edep = 3930
    emag = 3.24

if iex == 36:
    otime = util_helpers.eid2otime("20120609085418516")
    elat = 61.93
    elon = -150.59
    edep = 4020
    emag = 2.51

if iex == 37:
    otime = util_helpers.eid2otime("20101112032532608")
    elat = 62.16
    elon = -150.01
    edep = 4.59
    emag = 2.73

if iex == 38:
    otime = util_helpers.eid2otime("20121209152602476")
    elat = 62.32
    elon = -149.87
    edep = 5130
    emag = 2.51

if iex == 39:
    otime = util_helpers.eid2otime("20111203093358462")
    elat = 61.97
    elon = -150.93
    edep = 5710
    emag = 4.2

if iex == 40:
    otime = util_helpers.eid2otime("20151126124903509")
    elat = 62.00
    elon = -150.69
    edep = 6420
    emag = 2.57

if iex == 41:
    otime = util_helpers.eid2otime("20151126102704452")
    elat = 62.02
    elon = -150.72
    edep = 7300
    emag = 2.67

if iex == 42:
    otime = util_helpers.eid2otime("20130824022033577")
    elat = 61.97
    elon = -150.95
    edep = 7400
    emag = 3.35

if iex == 43:
    otime = util_helpers.eid2otime("20071219215856568")
    elat = 62.23
    elon = -150.13
    edep = 8210
    emag = 3.22

if iex == 44:
    otime = util_helpers.eid2otime("20120224210718441")
    elat = 61.81
    elon = -150.98
    edep = 8810
    emag = 2.60

if iex == 45:
    otime = util_helpers.eid2otime("20081124001242867")
    elat = 62.11
    elon = -150.55
    edep = 8950
    emag = 2.72

if iex == 46:
    otime = util_helpers.eid2otime("20111107082050735")
    elat = 62.16
    elon = -150.04
    edep = 9040
    emag = 2.58

if iex == 47:
    otime = util_helpers.eid2otime("20090105021152690")
    elat = 62.04
    elon = -150.33
    edep = 9130
    emag = 2.55

if iex == 48:
    otime = util_helpers.eid2otime("20130930063202102")
    elat = 61.92
    elon = -150.90
    edep = 9190
    emag = 3.50

if iex == 48:
    otime = util_helpers.eid2otime("20160622125133166")
    elat = 62.34
    elon = -149.96
    edep = 9450
    emag = 2.68

if iex == 49:
    otime = util_helpers.eid2otime("20150119103611702")
    elat = 62.19
    elon = -150.57
    edep = 9460
    emag = 4.01

if iex == 50:
    otime = util_helpers.eid2otime("20071128015236594")
    elat = 62.04
    elon = -150.47
    edep = 9460
    emag = 2.78

if iex == 51:
    otime = util_helpers.eid2otime("20160513193104479")
    elat = 62.14
    elon = -150.40
    edep = 9910
    emag = 2.96

if iex == 52:
    otime = util_helpers.eid2otime("20101201231944361")
    elat = 62.30
    elon = -150.11
    edep = 10000
    emag = 3.21

if iex == 53:
    otime = util_helpers.eid2otime("20160118040556098")
    elat = 62.10
    elon = -150.64
    edep = 10090
    emag = 4.5

if iex == 54:
    otime = util_helpers.eid2otime("20110416060141762")
    elat = 62.31
    elon = -149.99
    edep = 10620
    emag = 3.09

if iex == 55:
    otime = util_helpers.eid2otime("20150903162903880")
    elat = 62.24
    elon = -150.52
    edep = 10750
    emag = 2.74

if iex == 56:
    otime = util_helpers.eid2otime("20110405183024423")
    elat = 62.31
    elon = -150.03
    edep = 11710
    emag = 3.40

if iex == 57:
    otime = util_helpers.eid2otime("20130120215658808")
    elat = 62.19
    elon = -150.40
    edep = 14220
    emag = 3.41

if iex == 58:
    otime = util_helpers.eid2otime("20140121142920250")
    elat = 62.09
    elon = -150.37
    edep = 15710
    emag = 3.25

if iex == 59:
    otime = util_helpers.eid2otime("20101214022237651")
    elat = 62.28
    elon = -150.27
    edep = 19210
    emag = 3.70

if iex == 60:
    otime = util_helpers.eid2otime("20140902073120072")
    elat = 62.31
    elon = -150.40
    edep = 21110
    emag = 2.87

if iex == 61:
    otime = util_helpers.eid2otime("20150518154910522")
    elat = 61.94
    elon = -150.45
    edep = 21500
    emag = 4.27



# subset of stations
min_dist = 0
max_dist = 200
tbefore_sec = 100
tafter_sec = 300
network = 'AK,AT,AV,CN,II,IU,US,XM,XV,XZ,YV,ZE'
channel = 'BH?,HH?,EH?'
resample_freq = 50        # no resampling
scale_factor = 100         # no scale factor

#=================================================================================
# END EXAMPLES
#=================================================================================

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

# Delete existing data directory
eid = util_helpers.otime2eid(ev.origins[0].time)
ddir = './'+ eid
if os.path.exists('RAW'):
    print("WARNING. %s already exists. Deleting ..." % ddir)
    shutil.rmtree('RAW')
if overwrite_ddir and os.path.exists(ddir):
    print("WARNING. %s already exists. Deleting ..." % ddir)
    shutil.rmtree(ddir)

# Extract waveforms, IRIS
getwaveform.run_get_waveform(c = client, event = ev, idb = idb, 
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
                                  iplot_response = iplot_response, ifplot_spectrogram = ifplot_spectrogram)
