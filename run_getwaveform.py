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
# + provide better options and handling for data gaps (like: "toss waveform if gaps ar 0.01*length")
# + the KEVNM header cannot store the time to the 3rd millisecond character
#   probably the best approach is to write the spill-over character into another field
#   (or reconstruct the EID from the origin time, if that is your convention)
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
iex = 0
print("Running example iex =", iex)

# DEFAULT SETTINGS (see getwaveform_iris.py)
idb = 1    # default: =1-IRIS; =2-AEC; =3-LLNL
# Pre-processing (manily for CAP)
rotateRTZ = True
rotateUVW = False   # This option works only if 'rotateRTZ = True'
output_cap_weight_file = True
detrend = True
demean = True
taper = False       # this could also be a fraction between 0 and 1 (fraction to be tapered from both sides)
output_event_info = True
outformat = 'VEL'            # Intrument response removed waveforms could be saved as 'VEL' 'DISP' 'ACC'
ifsave_sacpaz = False        # save sac pole zero (needed as input for MouseTrap module)

# for CAP all waveforms need to have the same sample rate
resample_TF = True
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


# dummy values
dummyval = -9999
rlat = dummyval
rlon = dummyval
rtime = dummyval

# username and password for accessing embargoed data from IRIS
# Register here: http://ds.iris.edu/ds/nodes/dmc/forms/restricted-data-registration/
# Run example iex = 4 to check
user = ''
password = ''

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
    station = 'KDAK,COLA'
    channel = '*'   

#=================================================================================
# CATEGORY A EXAMPLES: simple test cases (including current and previous bugs)
#=================================================================================

# ERROR EXAMPLE [obspy]
# PROBLEM: No waveforms are returned -- perhaps related to the tbefore_sec request
# ERROR MESSAGE: ValueError: The length of the input vector x must be at least padlen, which is 39.
# SOLUTION: iex = 2 fails because the input data is much to short. The input is a trace with 38 sample but sampled at 50 Hz and you want to resample to 10 Hz. The IIR filter coefficients it calculates are just too long. I'm working on getting such a filter into ObsPy which should then have nicer error messages. I doubt its possible to meaningfully filter such short data with a very sharpy filter. For now: Just add a QA step that makes sure that at least a certain number of samples enter the decimation routine.
if iex == 1:
    otime = obspy.UTCDateTime("2016-01-24T10:30:29.557")
    min_dist = 0 
    max_dist = 500
    tbefore_sec = 42          # Crashes
    tbefore_sec = 41         # works fine
    tafter_sec = 600
    network = 'AV'       # Crashes when -  tbefore_sec = 42; Works fine when - tbefore_sec = 41
    channel = 'BH?'
    rlat = 64.7716
    rlon = -149.1465
    rtime = obspy.UTCDateTime("2009-04-07T20:20:00")

# ERROR EXAMPLE [obspy]
# PROBLEM: If a particular network is requested (either explicitly or within *), no waveforms are returned.
# ERROR MESSAGE: NotImplementedError: ResponseListResponseStage not yet implemented due to missing example data. Please contact the developers with a test data set (waveforms and StationXML metadata).
# KLUDGE: list all networks explicitly, except IM
# SOLUTION: https://github.com/obspy/obspy/issues/1514
#           If you want to use it right now you'd have to use that branch 
#           - it will be a while before we release a new ObsPy version.
if iex == 2:
    otime = obspy.UTCDateTime("2009-04-07T20:12:55")
    min_dist = 300 
    max_dist = 500
    tbefore_sec = 100
    tafter_sec = 300
    #network = '*'       # crashes
    #network = 'AK,IM'   # crashes (Error because of the IM network)
    #network = '*,-IM'   # stalls indefinitely (syntax works for stations but not networks)
    network = 'AK'       # works fine 
    channel = 'BH?'

# SALMON example (restricted data from IRIS)
if iex == 4:
    otime = obspy.UTCDateTime("2016-01-24T10:30:29.557")
    min_dist = 0 
    max_dist = 300
    tbefore_sec = 100
    tafter_sec = 600
    network = 'AK,ZE'   # ZE waveforms not returned (only ZE.MPEN)
    channel = 'BH?,HH?'

# ROTATION example for components 1,2,Z
#    All rotations should be based on the azimuth of the sensor
#    (CMPAZ header in sac), which, combined with the station-source backazimuth,
#    will provide the rotation to radial and transverse components.
#    The rotation should NOT depend on channel names like E or N.
#    (Even for GSN stations the E and N do not point exactly E and N.)
if iex == 5:
    otime = obspy.UTCDateTime("2016-01-24T10:30:29.557")
    tbefore_sec = 100
    tafter_sec = 600
    #network = 'II,AK'
    station = 'KDAK,SWD'
    channel = 'BH?'

# nuclear event: LLNL (see also iex = 7)
# GOAL: To find events in the LLNL database based on a target origin time, rather than an eid.
#       (The reference time (NZYEAR, etc) should then be assigned as the actual origin time,
#       not the target origin time.)
# DEBUGGING HELPER LINE:
#   saclst NPTS o b e NZHOUR NZMIN NZSEC NZMSEC f 19910914190000000/*.z
# (This will show clearly that the reference time is NOT the origin time.)
if iex == 6:
    # to get the LLNL client, which is a private repo from Lion Krischer:
    # cd $REPOS
    # git clone https://GITHUBUSERNAME@github.com/krischer/llnl_db_client.git
    # then follow instructions for install
    idb = 3              # LLNL database
    #resample_freq = 0    # no resampling and no cutting

    # TARGET origin time (8.031 s from actual origin time)
    otime = obspy.UTCDateTime("1991-09-14T19:00:08.031Z")
    #evid = 635527        # Hoya event id in LLNL database

    min_dist = 0 
    max_dist = 1200
    tbefore_sec = 100
    tafter_sec = 600
    network = '*'        # note that the client will look for BK stations in the list
    channel = 'BH?'      # ALL channels from LLNL are returned
    #scale_factor = 10.0**2  # original
    scale_factor = 2e-1     # Hoya  
    overwrite_ddir = 0

# same as iex=6 but for the IRIS database
# GOAL: For LLNL events, we do NOT want to use the IRIS source parameters:
#       origin time, hypocenter, magnitude.
if iex == 7:
    idb = 1            # IRIS database
    #resample_freq = 0  # no resampling
    #otime = obspy.UTCDateTime("1991-09-14T19:00:00.000Z")   # Hoya actual
    otime = obspy.UTCDateTime("1991-09-14T19:00:08.031Z")   # Hoya target
    min_dist = 0 
    max_dist = 1200
    tbefore_sec = 100
    tafter_sec = 600
    # needs to be run TWICE to get BK stations and IRIS stations
    #network = 'BK'        # BK will go to NCEDC
    network = '*'         # * will give all at IRIS DMC
    channel = 'BH?,LH?' 
    overwrite_ddir = 0

# problem 1: some stations return only vertical component. our tools crash in this case.
# problem 2: short waveforms. padding by zeros adds sharp changes in the data
# problem 3: waveform contains NAN or INF. this crashes detrend
# solution 1: (ongoing)
# solution 2: disable zero-padding
# solution 3: (ongoing) print error (consider removing trace?)
if iex==8:
    idb = 1            # IRIS database
    otime = obspy.UTCDateTime("1982-08-05T14:00:00.000000Z")
    min_dist = 0 
    max_dist = 1200
    tbefore_sec = 100
    tafter_sec = 600
    network = '*'
    channel = 'BH?,LH?,EH?' 

# 1. Waveform extraction for user defined event info
# 2. Subset of stations for quickly testing data gaps and padding tests
if iex == 9:
    otime = obspy.UTCDateTime("2009-04-07T20:12:55.351")
    elat = 61.4542
    elon =-149.7428
    edep = 33033.6  # in meters
    emag = 4.6
    min_dist = 0 
    max_dist = 500
    tbefore_sec = 100
    tafter_sec = 300
    network = 'AK,AT,AV,CN,II,IU,US,XM,XV,XZ,YV'  # note: cannot use '*' because of IM
    channel = 'BH?'
    network = 'AV'
    station = 'SPBG,KABU'                      # For testing data gaps 
    use_catalog = 0                            # To get (lat,lon, etime, dep, mag) from some catalog = 1 OR use defined = 0 (see iex=9)

# Error from util_write_cap.py
#   util_write_cap.py, line 437, in add_sac_metadata
#    and pick.waveform_id.channel_code[2].upper() == 'Z' and pick.waveform_id.location_code == tr.stats.location and pick.phase_hint == 'Pn'):
#        IndexError: string index out of range
#
# This error happens when allowing data already rotated (VRT). see getwaveform_llnl.py Line 94
if iex == 10:
    idb = 3              # LLNL database
    otime = obspy.UTCDateTime("1982-08-05T14:00:00.000000Z")
    min_dist = 0 
    max_dist = 1200
    tbefore_sec = 100
    tafter_sec = 600
    network = '*'        # note that the client will look for BK stations in the list
    channel = 'BH?'      # ALL channels from LLNL are returned regardless
    scale_factor = 10.0**2  # original
    overwrite_ddir = 0
    resample_freq = 20.0 

# error applying rotation to station LL.KNB. See util_write_cap.py Line 161
# (Disable the try block and only allow Line 161 to create this problem)
# 
#   File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/signal/rotate.py", line 257, in rotate2zne
#    x, y, z = np.dot(_t, [data_1, data_2, data_3])
#    ValueError: setting an array element with a sequence.
if iex == 11:
    idb = 1 
    otime = obspy.UTCDateTime("1995-07-31T12:34:46.860000Z")
    min_dist = 0 
    max_dist = 1200
    tbefore_sec = 100
    tafter_sec = 600
    network = '*'        # note that the client will look for BK stations in the list
    channel = 'BH?'      # ALL channels from LLNL are returned regardless
    scale_factor = 10.0**2  # original
    overwrite_ddir = 0
    resample_freq = 20.0 

# PROBLEM: 
#   File "/home/alvizuri/REPOSITORIES/GEOTOOLS/python_util/util_data_syn/getwaveform_iris.py", line 71, in run_get_waveform
#    level="response")
#  File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/clients/fdsn/client.py", line 620, in get_stations
#    inventory = read_inventory(data_stream)
#  File "<decorator-gen-81>", line 2, in read_inventory
#  File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/core/util/decorator.py", line 289, in _map_example_filename
#    return func(*args, **kwargs)
#  File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/core/inventory/inventory.py", line 72, in read_inventory
#    format=format)[0]
#  File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/core/util/base.py", line 459, in _read_from_plugin
#    list_obj = read_format(filename, **kwargs)
#  File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/io/stationxml/core.py", line 154, in _read_stationxml
#    networks.append(_read_network(network, _ns))
#  File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/io/stationxml/core.py", line 200, in _read_network
#    stations.append(_read_station(station, _ns))
#  File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/io/stationxml/core.py", line 236, in _read_station
#    channels.append(_read_channel(channel, _ns))
#  File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/io/stationxml/core.py", line 356, in _read_channel
#    channel.response = _read_response(response, _ns)
#  File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/io/stationxml/core.py", line 377, in _read_response
#    response.response_stages.append(_read_response_stage(stage, _ns))
#  File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/io/stationxml/core.py", line 441, in _read_response_stage
#    raise ValueError(msg)
#ValueError: Could not find a valid Response Stage Type.
if iex == 12:
    idb = 1 
    otime = obspy.UTCDateTime("1999-01-23T03:00:33.200000Z") # frenchman flat 1
    #otime = obspy.UTCDateTime("1999-01-27T10:44:23.310000Z") # frenchman flat 2
    #otime = obspy.UTCDateTime("2000-01-30T14:46:51.310000Z") # trona mine 2
    #otime = obspy.UTCDateTime("2002-06-14T12:40:44.450000Z") # little skull
    min_dist = 0 
    max_dist = 1200
    tbefore_sec = 100
    tafter_sec = 600
    network = '*'        # note that the client will look for BK stations in the list
    channel = 'BH?'      # ALL channels from LLNL are returned regardless
    scale_factor = 10.0**2  # original
    overwrite_ddir = 0
    resample_freq = 20.0 

# Traceback (most recent call last):
#   File "run_getwaveform.py", line 380, in <module>
#     scale_factor = scale_factor, pre_filt = pre_filt)
#   File "/home/alvizuri/REPOSITORIES/GEOTOOLS/python_util/util_data_syn/getwaveform_iris.py", line 95, in run_get_waveform
#     output="VEL")
#   File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/core/stream.py", line 3030, in remove_response
#     tr.remove_response(*args, **kwargs)
#   File "<decorator-gen-39>", line 2, in remove_response
#   File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/core/trace.py", line 231, in _add_processing_info
#     result = func(*args, **kwargs)
#   File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/core/trace.py", line 2701, in remove_response
#     output=output, **kwargs)
#   File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/core/inventory/response.py", line 762, in get_evalresp_response
#     raise ObsPyException(msg)
# obspy.core.util.obspy_types.ObsPyException: Can not use evalresp on response with no response stages.
if iex == 13:
    idb = 1 
    otime = obspy.UTCDateTime("1997-06-14T19:48:19.930000Z") # indian springs
    #otime = obspy.UTCDateTime("1996-09-05T08:16:56.090000Z") # amargosa
    min_dist = 0 
    max_dist = 1200
    tbefore_sec = 100
    tafter_sec = 600
    network = '*'        # note that the client will look for BK stations in the list
    channel = 'BH?'      # ALL channels from LLNL are returned regardless
    scale_factor = 10.0**2  # original
    overwrite_ddir = 0
    resample_freq = 20.0 

# time > 2007-01-24T11:30:06.100000Z time < 2007-01-24T11:30:26.100000Z
#Traceback (most recent call last):
#  File "run_getwaveform.py", line 478, in <module>
#    ev = cat.filter(mintime_str, maxtime_str)[0]
#  File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/core/event/catalog.py", line 162, in __getitem__
#    return self.events.__getitem__(index)
#IndexError: list index out of range
if iex == 14:
    idb = 3 
    otime = obspy.UTCDateTime("2007-01-24T11:30:16.100000Z") # ralston
    min_dist = 0 
    max_dist = 1200
    tbefore_sec = 100
    tafter_sec = 600
    network = '*'        # note that the client will look for BK stations in the list
    channel = 'BH?'      # ALL channels from LLNL are returned regardless
    scale_factor = 10.0**2  # original
    overwrite_ddir = 0
    resample_freq = 20.0 

#Traceback (most recent call last):
#  File "run_getwaveform.py", line 455, in <module>
#    scale_factor = scale_factor, pre_filt = pre_filt)
#  File "/home/alvizuri/REPOSITORIES/GEOTOOLS/python_util/util_data_syn/getwaveform_iris.py", line 214, in run_get_waveform
#    rotate_and_write_stream(st2, evname_key)
#  File "/home/alvizuri/REPOSITORIES/GEOTOOLS/python_util/util_data_syn/util_write_cap.py", line 162, in rotate_and_write_stream
#    data_array = rotate.rotate2zne(d1, az1, dip1, d2, az2, dip2, d3, az3, dip3)
#  File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/signal/rotate.py", line 257, in rotate2zne
#    x, y, z = np.dot(_t, [data_1, data_2, data_3])
#ValueError: setting an array element with a sequence.
if iex == 15:
    idb = 1 
    otime = obspy.UTCDateTime("1995-07-31T12:34:46.860000Z") # timber mountain
    min_dist = 0 
    max_dist = 1200
    tbefore_sec = 100
    tafter_sec = 600
    network = '*'        # note that the client will look for BK stations in the list
    channel = 'BH?'      # ALL channels from LLNL are returned regardless
    scale_factor = 10.0**2  # original
    overwrite_ddir = 0
    resample_freq = 20.0 

# Iniskin earthquake - all channels at IU.COLA and II.KDAK
# note different channels and location codes for strong motion:
#       IU.COLA.20.HNZ.sac       http://ds.iris.edu/mda/IU/COLA
#       II.KDAK.00.ENZ.sac       http://ds.iris.edu/mda/II/KDAK
if iex == 16:
    idb = 1
    overwrite_ddir = 1       # delete data dir if it exists
    use_catalog = 0          # do not use event catalog for source parameters
    # GCMT source parameters
    # the otime is the centroid time and accounts for tshift
    otime = obspy.UTCDateTime("2016-01-24T10:30:37.400") 
    elat = 59.75
    elon = -153.27
    edep = 110700
    emag = 7.1
    # subset of stations
    min_dist = 0
    max_dist = 800
    tbefore_sec = 100
    tafter_sec = 600
    network = 'IU,II'
    #channel = '?H?,?N?'
    channel = 'HH?,BH?,BN?,HN?,EN?'
    resample_freq = 0        # no resampling
    scale_factor = 1         # no scale factor

# PROBLEM data arrays not of the same length
# Traceback (most recent call last):
#   File "run_getwaveform.py", line 736, in <module>
#     iplot_response = iplot_response)
#   File "/home/alvizuri/REPOSITORIES/GEOTOOLS/python_util/util_data_syn/getwaveform_llnl.py", line 166, in run_get_waveform
#     rotate_and_write_stream(st2, evname_key, icreateNull)
#   File "/home/alvizuri/REPOSITORIES/GEOTOOLS/python_util/util_data_syn/util_write_cap.py", line 167, in rotate_and_write_stream
#     data_array = rotate.rotate2zne(d1, az1, dip1, d2, az2, dip2, d3, az3, dip3)
#   File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/signal/rotate.py", line 247, in rotate2zne
#     raise ValueError(msg)
# ValueError: All three data arrays must be of same length.
if iex == 17:
    idb = 3            # LLNL
    otime = obspy.UTCDateTime("1989-06-27T15:30:00.02")
    elat = 37.275
    elon = -116.354
    edep = 640
    emag = 4.90
    min_dist = 0 
    max_dist = 1200
    tbefore_sec = 100
    tafter_sec = 600
    network = '*' 
    channel = 'BH?,LH?' 
    overwrite_ddir = 0

# for debugging:
# python -m pdb run_getwaveform.py
# Keep pressing `c` when prompted. It will drop you into the debugger as soon as it encounters an error. At that point you can hop up and down the stack with "u" and "d" and print the variables to see where it goes wrong. We can also get together quickly today and find the source of the issue.
#
# Occasional error
# Similar as iex=32 but for the IRIS database
#   File "/home/vipul/REPOSITORIES/GEOTOOLS/python_util/util_data_syn/util_write_cap.py", line 167, in rotate_and_write_stream
#    data_array = rotate.rotate2zne(d1, az1, dip1, d2, az2, dip2, d3, az3, dip3)
#  File "/home/vipul/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/signal/rotate.py", line 257, in rotate2zne
#    x, y, z = np.dot(_t, [data_1, data_2, data_3])
# ValueError: setting an array element with a sequence
#
# NOTE: Rerunning the same script (without any changes) solves the error sometimes!
#
# KEY:
# The problem arises from calling resample within util_write_cap.py.
# (If we do not call resample at all, or if we use resample_cut [util_write_cap.py], it works.)
#
# Not generating RTZ and NEZ waveforms after implementing Lion's suggestions
if iex == 18:
    idb = 1
    overwrite_ddir = 1       # delete data dir if it exists
    use_catalog = 0          # do not use event catalog for source parameters
    # GCMT source parameters
    # the otime is the centroid time and accounts for tshift
    otime = obspy.UTCDateTime("2015-12-02T10:05:25.798") 
    elat = 61.70
    elon = -147.26
    edep = 36590
    emag = 4.50
    # subset of stations
    min_dist = 300
    max_dist = 400
    tbefore_sec = 100
    tafter_sec = 600
    #network = 'AK,AT,AV,CN,II,IU,US,XM,XV,XZ,YV,ZE'
    network = 'IU'
    channel = 'HH?,BH?'
    #resample_freq = 0        # no resampling -- THIS WORKS
    resample_freq = 50       # THIS FAILS BUT ONLY AFTER REPEATING
    scale_factor = 1         # no scale factor

# Test case for UVW
if iex == 19:
    idb = 1
    use_catalog = 0  
    otime = obspy.UTCDateTime("2016-12-08T10:16:00")
    elat = 64.2380
    elon = -150.0581
    edep = 18507
    emag = 4.60
    station = 'F3TN'

    min_dist = 0
    max_dist = 300
    tbefore_sec = 0
    tafter_sec = 600
    network = 'XV'
    channel = 'HH?'
    scale_factor = 1
    resample_freq = 0
    detrend = False
    demean = False
    taper = False
    ipre_filt = 1

# PROBLEM run_getwaveform crashes when processing LLNL data for event BULLION
#
# WARNING: 0 traces available for rotation. Adding NULL traces -  LL.TPH..SHR*
# Traceback (most recent call last):
#   File "run_getwaveform.py", line 858, in <module>
#     iplot_response = iplot_response)
#   File "/home/alvizuri/REPOSITORIES/GEOTOOLS/python_util/util_data_syn/getwaveform.py", line 180, in run_get_waveform
#     rotate_and_write_stream(st2, evname_key, icreateNull)
#   File "/home/alvizuri/REPOSITORIES/GEOTOOLS/python_util/util_data_syn/util_write_cap.py", line 133, in rotate_and_write_stream
#     d1 = substr[0].data
#   File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/core/stream.py", line 656, in __getitem__
#     return self.traces.__getitem__(index)
# IndexError: list index out of range
if iex == 20:
    idb = 3            # LLNL
    otime = obspy.UTCDateTime("1990-06-13T16:00:00.09")
    elat = 37.262
    elon = -116.421
    edep = 674
    emag = 5.34
    min_dist = 0 
    max_dist = 1200
    tbefore_sec = 100
    tafter_sec = 600
    network = '*' 
    channel = 'BH?,LH?' 
    overwrite_ddir = 0

#=================================================================================
# CATEGORY B EXAMPLES: important events
# 1XX: southern Alaska
# 2XX: central Alaska
# 9XX: other
#=================================================================================

# SilwalTape2016 example event (Anchorage)
if iex == 100:
    otime = obspy.UTCDateTime("2009-04-07T20:12:55.351")
    min_dist = 0 
    max_dist = 500
    tbefore_sec = 100
    tafter_sec = 300
    network = 'AK,AT,AV,CN,II,IU,US,XM,XV,XZ,YV'  # note: cannot use '*' because of IM
    channel = 'BH?'
    use_catalog = 0 
    elat = 61.45420
    elon = -149.7428
    edep = 33033.60
    emag = 4.6
    #ipre_filt = 2
    #pre_filt = (0.005, 0.006, 10.0, 15.0)

# Iniskin earthquake
# NOTE: must enter username and password above to get SALMON (ZE) stations
if iex == 101:
    idb = 1
    overwrite_ddir = 1       # delete data dir if it exists
    use_catalog = 0          # do not use event catalog for source parameters
    # GCMT source parameters
    # the otime is the centroid time and accounts for tshift
    otime = obspy.UTCDateTime("2016-01-24T10:30:37.400") 
    elat = 59.75
    elon = -153.27
    edep = 110700
    emag = 7.1
    # subset of stations
    min_dist = 0
    max_dist = 800
    tbefore_sec = 100
    tafter_sec = 600
    network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US'  # IM will probably crash it
    channel = 'BH?,HH?'
    resample_freq = 0        # no resampling
    scale_factor = 1         # no scale factor

    # parameters for examining step response (causal low-pass filter on raw waveforms)
    # delete AUQ, SPCP, SPBG
    rotateRTZ = False
    ifFilter = True
    filter_type = 'lowpass'
    f1 = 1/4
    zerophase = False
    remove_response = False
    ipre_filt = 0
    demean = False
    detrend = False

# Iniskin earthquake - all strong motion
if iex == 103:
    idb = 1
    overwrite_ddir = 1       # delete data dir if it exists
    use_catalog = 0          # do not use event catalog for source parameters
    # GCMT source parameters
    # the otime is the centroid time and accounts for tshift
    otime = obspy.UTCDateTime("2016-01-24T10:30:37.400") 
    elat = 59.75
    elon = -153.27
    edep = 110700
    emag = 7.1
    # subset of stations
    min_dist = 0
    max_dist = 800
    tbefore_sec = 100
    tafter_sec = 600
    network = '*'
    channel = 'BN?,HN?,EN?'
    resample_freq = 0        # no resampling
    scale_factor = 1         # no scale factor

# Totschunda fault
if iex == 104:
    idb = 1
    overwrite_ddir = 1       # delete data dir if it exists
    use_catalog = 0          # do not use event catalog for source parameters
    # GCMT source parameters
    # the otime is the centroid time and accounts for tshift
    otime = obspy.UTCDateTime("2017-02-13T07:17:12.000") 
    elat = 62.5154
    elon = -142.7485
    edep = 7564
    emag = 5.3
    # subset of stations
    min_dist = 0
    max_dist = 500
    tbefore_sec = 50
    tafter_sec = 300
    network = 'AV'
    channel = 'BH?,HH?'
    resample_freq = 50        # no resampling
    scale_factor = 100         # no scale factor

# Klukwan earthquakes
if iex == 105:
    idb = 1
    overwrite_ddir = 1       # delete data dir if it exists
    use_catalog = 0          # do not use event catalog for source parameters
    # GCMT source parameters
    # the otime is the centroid time and accounts for tshift
    otime = obspy.UTCDateTime("2017-05-01T12:31:53.000") 
    elat = 59.8522
    elon = -136.6618
    edep = 4000
    emag = 6.2
    otime = obspy.UTCDateTime("2017-05-01T14:18:14.000") 
    elat = 59.8184
    elon = -136.7163
    edep = 1000
    emag = 6.0
    
    # subset of stations
    min_dist = 0
    max_dist = 500
    tbefore_sec = 100
    tafter_sec = 500
    network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
    channel = 'BH?,HH?'
    resample_freq = 50        # no resampling
    scale_factor = 100         # no scale factor

# Kantishna earthquakes
if iex == 106:
    idb = 1
    overwrite_ddir = 1       # delete data dir if it exists
    use_catalog = 0          # do not use event catalog for source parameters
    # GCMT source parameters
    # the otime is the centroid time and accounts for tshift
    otime = obspy.UTCDateTime("2017-04-29T11:15:48.000") 
    elat = 63.1296
    elon = -151.1517
    edep = 11000
    emag = 5.2
    otime = obspy.UTCDateTime("2017-01-31T09:38:37.000") 
    elat = 63.0817
    elon = -150.9427
    edep = 135000
    emag = 5.2
    
    # subset of stations
    min_dist = 0
    max_dist = 500
    tbefore_sec = 100
    tafter_sec = 500
    network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
    channel = 'BH?,HH?'
    resample_freq = 50        # no resampling
    scale_factor = 100         # no scale factor

# Cook Inlet earthquake
if iex == 107:
    idb = 1
    overwrite_ddir = 1       # delete data dir if it exists
    use_catalog = 0          # do not use event catalog for source parameters
    # GCMT source parameters
    # the otime is the centroid time and accounts for tshift
    otime = obspy.UTCDateTime("2017-05-07T04:25:19.000") 
    elat = 60.1945
    elon = -151.6743
    edep = 64000
    emag = 5.2
    
    # subset of stations
    min_dist = 0
    max_dist = 500
    tbefore_sec = 100
    tafter_sec = 500
    network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
    channel = 'BH?,HH?'
    resample_freq = 50        # no resampling
    scale_factor = 100         # no scale factor


# MFFZ earthquakes for investigating the step response
# LISTED IN CHRONOLOGICAL ORDER
if iex == 200:
    idb = 1
    use_catalog = 0
    #-------------------------------------------------
    otime = obspy.UTCDateTime("2014-08-31T03:06:57.111")
    elat = 65.1526
    elon = -149.0398
    edep = 16614.7
    emag = 5.20
    #-------------------------------------------------
    otime = obspy.UTCDateTime("2014-10-21T00:36:58.333")
    elat = 65.1489
    elon = -149.0413
    edep = 13134.8
    emag = 4.90
    #-------------------------------------------------
    otime = obspy.UTCDateTime("2014-10-23T16:30:23.968")
    elat = 65.1644
    elon = -149.0523
    edep = 20066.5
    emag = 5.00
    #-------------------------------------------------
    otime = obspy.UTCDateTime("2015-09-12T03:25:12.711")
    elat = 65.1207
    elon = -148.6646
    edep = 15556.8
    emag = 2.63
    #-------------------------------------------------
    otime = obspy.UTCDateTime("2015-10-31T02:56:35.572")
    elat = 64.4285
    elon = -149.6969
    edep = 23852.1
    emag = 3.47
    #-------------------------------------------------
    otime = obspy.UTCDateTime("2016-01-14T19:04:10.727")
    elat = 64.6827
    elon = -149.2479
    edep = 22663.7
    emag = 3.80
    #-------------------------------------------------
    otime = obspy.UTCDateTime("2016-11-06T9:29:10.579")
    elat = 64.1639
    elon = -150.0626
    edep = 23190.0
    emag = 4.00
    #-------------------------------------------------
    otime = obspy.UTCDateTime("2016-12-08T10:18:13.868")
    elat = 64.1937
    elon = -150.0376
    edep = 24522.1
    emag = 4.30
    #-------------------------------------------------
    otime = obspy.UTCDateTime("2017-05-08T05:09:02.000") 
    elat = 65.2643
    elon = -146.922
    edep = 9000   # AEC/Vipul
    emag = 3.60   # Vipul
    #-------------------------------------------------
    otime = obspy.UTCDateTime("2017-06-28T12:58:52.000") 
    elat = 64.7569
    elon = -148.8883
    edep = 18000  # Vipul
    emag = 3.50   # Vipul
    #-------------------------------------------------
    # VIPUL: WHAT ARE THESE? OTHER EVENTS?
    #otime = obspy.UTCDateTime("2015-03-30T12:33:19.000")
    #otime = obspy.UTCDateTime("2015-10-20T19:14:16.000")
    #otime = obspy.UTCDateTime("2011-12-21T16:28:41.000")
    #-------------------------------------------------

    min_dist = 0
    max_dist = 300
    tbefore_sec = 200
    tafter_sec = 600
    network = 'AV,CN,AT,TA,AK,XV,II,IU,US'
    #network = 'XV,AK'
    channel = 'BH?,HH?'
    scale_factor = 1
    resample_freq = 0
    # for CAP
    #scale_factor = 10.0**2
    #resample_freq = 50 

    # to investigate step response
    rotateRTZ = True
    rotateUVW = True
    ifFilter = True
    zerophase = False    # causal
    #filter_type = 'lowpass'
    filter_type = 'bandpass'
    f1 = 1/100  # fmin
    f2 = 1/10  # fmax
    corners = 4
    #remove_response = False
    remove_response = True
    ipre_filt = 1
    demean = True
    detrend = True
    output_cap_weight_file = False
    #outformat = 'DISP'
    ifsave_sacpaz = True
    taper = 0.2

# NENNUC event (from Steve)
if iex == 201:
    idb = 1
    use_catalog = 0          # manually enter AEC catalog parameters
    otime = obspy.UTCDateTime("2016-01-14T19:04:10.727")
    elat = 64.6827
    elon = -149.2479
    edep = 22663.7
    emag = 3.8
    #-------------------------------------------------
    otime = obspy.UTCDateTime("2015-09-12T03:25:12.711")
    elat = 65.1207
    elon = -148.6646
    edep = 15556.8
    emag = 2.6
    #-------------------------------------------------
    otime = obspy.UTCDateTime("2013-03-12T07:39:50.214")
    elat = 64.7161
    elon = -148.9505
    edep = 20000
    emag = 2.1

    # For CAP
    min_dist = 0
    max_dist = 200
    tbefore_sec = 500
    tafter_sec = 2000
    taper = 0.1
    network = 'AV,CN,AT,TA,AK,XV,II,IU,US'
    # network = 'AK,TA,II,IU,US'
    # station = "-RAG"  # RAG has choppy data; gives error: Exception: Can't merge traces with same ids but differing sampling rates!
    channel = 'BH?,HH?'
    pre_filt = (0.0025, 0.003, 10.0, 15.0) # BH default
    ipre_filt = 1

    # to investigate step response
    #tbefore_sec = 2000
    #tafter_sec = 2000
    #resample_freq = 0 
    #scale_factor = 1
    #ifFilter = True
    #zerophase = False
    #filter_type = 'bandpass'
    #f1 = 1/100
    #f2 = 1/20
    #remove_response = True
    #ipre_filt = 0
    #demean = True
    #detrend = True

# Kyle Nenana basin earthquakes for basin amplification
if iex == 210:
    idb = 1
    overwrite_ddir = 1       # delete data dir if it exists
    use_catalog = 0          # do not use event catalog for source parameters
    # AEC source parameters
    otime = obspy.UTCDateTime("2015-11-20T10:53:48.168") 
    elat = 64.6210
    elon = -149.4024
    edep = 17113.4
    emag = 2.67
    # subset of stations
    min_dist = 0
    max_dist = 100
    tbefore_sec = 100
    tafter_sec = 400
    network = 'AK,AT,II,IU,US,XM,XV,XZ,TA'  # no CN,AV,YV,ZE
    channel = 'BH?,HH?'
    #resample_freq = 0        # no resampling
    #scale_factor = 1         # no scale factor
    # For CAP moment tensor
    resample_freq = 50         # same as greens function 
    scale_factor = 100         # change from m/s to cm/s
    #ipre_filt = 0
    remove_response = True
    #demean = False
    #detrend = False

# gets waveforms from M 8.3 Chile event with stations centered in Minto
if iex == 211:
    idb = 1
    overwrite_ddir = 1       # delete data dir if it exists
    use_catalog = 0
    otime = obspy.UTCDateTime("2015-09-16T22:54:33.000")
    elat = -31.5695
    elon = -71.6543
    edep = 22400.0
    emag = 8.30
    rlat = 64.7716
    rlon = -149.1465
    rtime = obspy.UTCDateTime("2015-09-16T23:09:15.000")
    tbefore_sec = 100
    tafter_sec = 200
    min_dist = 0
    max_dist = 100
    network = 'AK,AT,II,IU,US,XM,XV,XZ,TA'  # no CN,AV,YV,ZE
    channel = 'BH?,HH?'
    resample_freq = 0        
    scale_factor = 1        
    remove_response = True

# gets a 
if iex == 212:
    idb = 1
    overwrite_ddir = 1       # delete data dir if it exists
    use_catalog = 0
    otime = obspy.UTCDateTime("2016-01-01T00:00:00.000")
    elat = 64.6130
    elon = -149.0992
    edep = 0
    emag = 0.00
    #rlat = 64.7716
    #rlon = -149.1465
    #rtime = obspy.UTCDateTime("2015-09-16T23:09:15.000")
    tbefore_sec = 0
    tafter_sec = 3600
    min_dist = 0
    max_dist = 100
    network = 'XV,AK,TA'  # no CN,AV,YV,ZE
    channel = 'HH?,BH?'
    resample_freq = 0        
    scale_factor = 1        
    remove_response = True
    rotateRTZ = False

# Chatnika earthquake
if iex == 213:
    idb = 1
    overwrite_ddir = 1       # delete data dir if it exists
    use_catalog = 0          # do not use event catalog for source parameters
    # GCMT source parameters
    # the otime is the centroid time and accounts for tshift
    otime = obspy.UTCDateTime("2017-05-08T05:09:02.000") 
    elat = 65.2643
    elon = -146.922
    edep = 9000
    emag = 3.8
    
    # subset of stations
    min_dist = 0
    max_dist = 500
    tbefore_sec = 100
    tafter_sec = 500
    network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
    channel = 'BH?,HH?'
    resample_freq = 50        # no resampling
    scale_factor = 100         # no scale factor

# NE Nenana earthquake
if iex == 214:
    idb = 1
    overwrite_ddir = 1       # delete data dir if it exists
    use_catalog = 0          # do not use event catalog for source parameters
    # GCMT source parameters
    # the otime is the centroid time and accounts for tshift
    otime = obspy.UTCDateTime("2017-06-28T12:58:52") 
    elat = 64.7569
    elon = -148.8883
    edep = 18000
    emag = 3.5
    
    # subset of stations
    min_dist = 0
    max_dist = 500
    tbefore_sec = 100
    tafter_sec = 500
    network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
    channel = 'BH?,HH?'
    resample_freq = 50        
    scale_factor = 100        
#------------------------------------------------

#-----------------------------------------------------------
# Event from HutchisonGhosh2016
# Crashes when using all networks [i.e. network = '*']
if iex == 900:
    idb = 1
    overwrite_ddir = 1       # delete data dir if it exists
    use_catalog = 0          # do not use event catalog for source parameters
    # GCMT source parameters
    # the otime is the centroid time and accounts for tshift
    otime = obspy.UTCDateTime("2014-12-12T00:13:17") 
    elat = 49.10
    elon = -124.05
    edep = 60000
    emag = 4.10
    # subset of stations
    min_dist = 0
    max_dist = 250
    tbefore_sec = 100
    tafter_sec = 300
    network = '*'
    channel = 'LH?,BH?'
    resample_freq = 50        # no resampling
    scale_factor = 100         # no scale factor

    # For plotting filtered waveforms
    tbefore_sec = 500
    tafter_sec = 500
    resample_freq = 0 
    scale_factor = 1
    ifFilter = True
    zerophase = True
    filter_type = 'bandpass'
    f1 = 1/50
    f2 = 1/20
    remove_response = True
    ipre_filt = 2
    demean = True
    detrend = True
    taper = True

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
        
        ref_time_place = ev
        if rlat != dummyval:
            ref_time_place.origins[0].latitude = rlat
            ref_time_place.origins[0].longitude = rlon
            ref_time_place.origins[0].time = rtime 
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

        if rlat == dummyval:
            # By default this should be the event time and location unless we want to grab stations centered at another location
            rlat = elat
            rlon = elon
            rtime = otime
        
        ref_time_place = Event()
        ref_org = Origin()
        ref_org.latitude = rlat
        ref_org.longitude = rlon
        ref_org.time = rtime
        ref_org.depth = 0 # dummy value
        ref_time_place.origins.append(ref_org)
        ref_time_place.magnitudes.append(mag) # more dummies

# LLNL
if idb == 3:
    import llnl_db_client
    #import getwaveform_llnl
    client = llnl_db_client.LLNLDBClient(
            "/store/raw/LLNL/UCRL-MI-222502/westernus.wfdisc")

    # get event time and event ID
    cat = client.get_catalog()
    mintime_str = "time > %s" % (otime - sec_before_after_event)
    maxtime_str = "time < %s" % (otime + sec_before_after_event)
    print(mintime_str + "\n" + maxtime_str)
    #ev = cat.filter(mintime_str, maxtime_str)[0]
    ev = cat.filter(mintime_str, maxtime_str)
    
    if len(ev) > 0:
        ev = ev[0]
        # Nothing happens here.  We can change later
        ref_time_place = ev
        print(len(ev))
    else:
        print("No events in the catalog for the given time period. Stop.")
        sys.exit(0)

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
                             network = network, station = station, channel = channel, ifresample = resample_TF,
                             resample_freq = resample_freq, ifrotateRTZ = rotateRTZ, ifrotateUVW = rotateUVW,
                             ifCapInp = output_cap_weight_file, 
                             ifRemoveResponse = remove_response,
                             ifDetrend = detrend, ifDemean = demean, Taper = taper,
                             ifEvInfo = output_event_info,
                             scale_factor = scale_factor,
                             ipre_filt = ipre_filt, pre_filt = pre_filt, 
                             icreateNull=icreateNull,
                             ifFilter = ifFilter, fmin = f1, fmax = f2, filter_type = filter_type, 
                             zerophase = zerophase, corners = corners, 
                             iplot_response = iplot_response, ifplot_spectrogram = ifplot_spectrogram,
                             outformat = outformat, ifsave_sacpaz = ifsave_sacpaz)
