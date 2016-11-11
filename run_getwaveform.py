# To run this script:
# > python run_getwaveform.py
import obspy
import copy
import util_helpers
import shutil   # only used for deleting data directory
import os

# EXAMPLES (choose one)
iex = 0

# DEFAULT SETTINGS (see getwaveform_iris.py)
idb = 1    # default: =1-IRIS; =2-AEC; =3-LLNL
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
overwrite_ddir = 0          # 1 = delete data directory if it already exists
icreateNull = 1
#------Filter--------------
ifFilter = True
filt_type = 'bandpass'
fmin = .02
fmax = .1
zerophase = False     # = False (causal); = True (acausal); 
corners = 4    # Is corner in Obspy same as Pole in SAC?

# username and password for accessing embargoed data from IRIS
# Register here: http://ds.iris.edu/ds/nodes/dmc/forms/restricted-data-registration/
# Run example iex = 4 to check
user = ''
password = ''

# (Do not use)
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

# SilwalTape2016 example event (Anchorage)
if iex == 1:
    otime = obspy.UTCDateTime("2009-04-07T20:12:55")
    min_dist = 0 
    max_dist = 500
    tbefore_sec = 100
    tafter_sec = 300
    network = 'AK,AT,AV,CN,II,IU,US,XM,XV,XZ,YV'  # note: cannot use '*' because of IM
    channel = 'BH?'

# ERROR EXAMPLE [obspy]
# PROBLEM: No waveforms are returned -- perhaps related to the tbefore_sec request
# ERROR MESSAGE: ValueError: The length of the input vector x must be at least padlen, which is 39.
# SOLUTION: iex = 2 fails because the input data is much to short. The input is a trace with 38 sample but sampled at 50 Hz and you want to resample to 10 Hz. The IIR filter coefficients it calculates are just too long. I'm working on getting such a filter into ObsPy which should then have nicer error messages. I doubt its possible to meaningfully filter such short data with a very sharpy filter. For now: Just add a QA step that makes sure that at least a certain number of samples enter the decimation routine.
if iex == 2:
    otime = obspy.UTCDateTime("2016-01-24T10:30:29.557")
    min_dist = 0 
    max_dist = 500
    tbefore_sec = 42          # Crashes
    tbefore_sec = 41         # works fine
    tafter_sec = 600
    network = 'AV'       # Crashes when -  tbefore_sec = 42; Works fine when - tbefore_sec = 41
    channel = 'BH?'

# ERROR EXAMPLE [obspy]
# PROBLEM: If a particular network is requested (either explicitly or within *), no waveforms are returned.
# ERROR MESSAGE: NotImplementedError: ResponseListResponseStage not yet implemented due to missing example data. Please contact the developers with a test data set (waveforms and StationXML metadata).
# KLUDGE: list all networks explicitly, except IM
# SOLUTION: https://github.com/obspy/obspy/issues/1514
#           If you want to use it right now you'd have to use that branch 
#           - it will be a while before we release a new ObsPy version.
if iex == 3:
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
    pre_filt = (0.005, 0.006, 10.0, 15.0)

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
    pre_filt = (0.005, 0.006, 10.0, 15.0)

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
    pre_filt = (0.005, 0.006, 10.0, 15.0)

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
    min_dist = 0 
    max_dist = 1200
    tbefore_sec = 100
    tafter_sec = 600
    network = '*'        # note that the client will look for BK stations in the list
    channel = 'BH?'      # ALL channels from LLNL are returned regardless
    scale_factor = 10.0**2  # original
    overwrite_ddir = 0
    resample_freq = 20.0 
    pre_filt = (0.005, 0.006, 10.0, 15.0)

# same as #12
# ValueError: Could not find a valid Response Stage Type.
if iex == 13:
    idb = 1 
    otime = obspy.UTCDateTime("1999-01-27T10:44:23.310000Z") # frenchman flat 2
    min_dist = 0 
    max_dist = 1200
    tbefore_sec = 100
    tafter_sec = 600
    network = '*'        # note that the client will look for BK stations in the list
    channel = 'BH?'      # ALL channels from LLNL are returned regardless
    scale_factor = 10.0**2  # original
    overwrite_ddir = 0
    resample_freq = 20.0 
    pre_filt = (0.005, 0.006, 10.0, 15.0)

# same as #12
# ValueError: Could not find a valid Response Stage Type.
if iex == 14:
    idb = 1 
    otime = obspy.UTCDateTime("2002-06-14T12:40:44.450000Z") # little skull
    min_dist = 0 
    max_dist = 1200
    tbefore_sec = 100
    tafter_sec = 600
    network = '*'        # note that the client will look for BK stations in the list
    channel = 'BH?'      # ALL channels from LLNL are returned regardless
    scale_factor = 10.0**2  # original
    overwrite_ddir = 0
    resample_freq = 20.0 
    pre_filt = (0.005, 0.006, 10.0, 15.0)

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
if iex == 15:
    idb = 1 
    otime = obspy.UTCDateTime("1997-06-14T19:48:19.930000Z") # indian springs
    min_dist = 0 
    max_dist = 1200
    tbefore_sec = 100
    tafter_sec = 600
    network = '*'        # note that the client will look for BK stations in the list
    channel = 'BH?'      # ALL channels from LLNL are returned regardless
    scale_factor = 10.0**2  # original
    overwrite_ddir = 0
    resample_freq = 20.0 
    pre_filt = (0.005, 0.006, 10.0, 15.0)

# time > 2007-01-24T11:30:06.100000Z time < 2007-01-24T11:30:26.100000Z
#Traceback (most recent call last):
#  File "run_getwaveform.py", line 478, in <module>
#    ev = cat.filter(mintime_str, maxtime_str)[0]
#  File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/core/event/catalog.py", line 162, in __getitem__
#    return self.events.__getitem__(index)
#IndexError: list index out of range
if iex == 16:
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
    pre_filt = (0.005, 0.006, 10.0, 15.0)

# same as #12
# ValueError: Could not find a valid Response Stage Type.
if iex == 17:
    idb = 1 
    otime = obspy.UTCDateTime("2000-01-30T14:46:51.310000Z") # trona mine 2
    min_dist = 0 
    max_dist = 1200
    tbefore_sec = 100
    tafter_sec = 600
    network = '*'        # note that the client will look for BK stations in the list
    channel = 'BH?'      # ALL channels from LLNL are returned regardless
    scale_factor = 10.0**2  # original
    overwrite_ddir = 0
    resample_freq = 20.0 
    pre_filt = (0.005, 0.006, 10.0, 15.0)

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
if iex == 18:
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
    pre_filt = (0.005, 0.006, 10.0, 15.0)

# same as # 15
# obspy.core.util.obspy_types.ObsPyException: Can not use evalresp on response with no response stages.
if iex == 19:
    idb = 1 
    otime = obspy.UTCDateTime("1996-09-05T08:16:56.090000Z") # amargosa
    min_dist = 0 
    max_dist = 1200
    tbefore_sec = 100
    tafter_sec = 600
    network = '*'        # note that the client will look for BK stations in the list
    channel = 'BH?'      # ALL channels from LLNL are returned regardless
    scale_factor = 10.0**2  # original
    overwrite_ddir = 0
    resample_freq = 20.0 
    pre_filt = (0.005, 0.006, 10.0, 15.0)

# Iniskin earthquake
# NOTE: must enter username and password above to get SALMON (ZE) stations
if iex == 20:
    idb = 1
    overwrite_ddir = 1       # delete data dir if it exists
    use_catalog = 0          # do not use event catalog for source parameters
    # GCMT source parameters
    # the otime is the centroid time and accounts for tshift
    otime = obspy.UTCDateTime("2016-01-24T10:30:37.400") 
    elat = 59.75
    elon = -153.27
    edep = 110700  # in meters
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
    pre_filt = (0.005, 0.006, 10.0, 15.0)

# Iniskin earthquake
# NOTE: must enter username and password above to get SALMON (ZE) stations
if iex == 20:
    idb = 1
    overwrite_ddir = 1       # delete data dir if it exists
    use_catalog = 0          # do not use event catalog for source parameters
    # GCMT source parameters
    # the otime is the centroid time and accounts for tshift
    otime = obspy.UTCDateTime("2016-01-24T10:30:37.400") 
    elat = 59.75
    elon = -153.27
    edep = 110700  # in meters
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
    pre_filt = (0.005, 0.006, 10.0, 15.0)

# Iniskin earthquake - all strong motion
if iex == 21:
    idb = 1
    overwrite_ddir = 1       # delete data dir if it exists
    use_catalog = 0          # do not use event catalog for source parameters
    # GCMT source parameters
    # the otime is the centroid time and accounts for tshift
    otime = obspy.UTCDateTime("2016-01-24T10:30:37.400") 
    elat = 59.75
    elon = -153.27
    edep = 110700  # in meters
    emag = 7.1
    # subset of stations
    min_dist = 0
    max_dist = 800
    tbefore_sec = 100
    tafter_sec = 600
    network = '*'
    channel = '?N?'
    resample_freq = 0        # no resampling
    scale_factor = 1         # no scale factor
    pre_filt = (0.005, 0.006, 10.0, 15.0)   # WHAT SHOULD THIS BE?

# Iniskin earthquake - all channels at IU.COLA and II.KDAK
# note: strong motion at COLA is IU.COLA.20.HNZ.sac
#       strong motion at KDAK is II.KDAK.00.ENZ.sac
# http://ds.iris.edu/mda/II/KDAK
# http://ds.iris.edu/mda/IU/COLA
if iex == 22:
    idb = 1
    overwrite_ddir = 1       # delete data dir if it exists
    use_catalog = 0          # do not use event catalog for source parameters
    # GCMT source parameters
    # the otime is the centroid time and accounts for tshift
    otime = obspy.UTCDateTime("2016-01-24T10:30:37.400") 
    elat = 59.75
    elon = -153.27
    edep = 110700  # in meters
    emag = 7.1
    # subset of stations
    min_dist = 0
    max_dist = 800
    tbefore_sec = 100
    tafter_sec = 600
    network = 'IU,II'
    channel = '?H?,?N?'
    resample_freq = 0        # no resampling
    scale_factor = 1         # no scale factor
    pre_filt = (0.005, 0.006, 10.0, 15.0)   # WHAT SHOULD THIS BE?

# MFFZ earthquake near Clear - Mw 4.1 
if iex == 30:
    idb = 1
    otime = obspy.UTCDateTime("2016-11-06T9:29:10")
    min_dist = 0
    max_dist = 300
    tbefore_sec = 100
    tafter_sec = 300
    network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US'
    channel = 'BH?,HH?'
    scale_factor = 10.0**2  # original
    overwrite_ddir = 1
    resample_freq = 50 
    pre_filt = (0.005, 0.006, 10.0, 15.0)

# fetch and process waveforms
if idb == 1:
    # import functions to access waveforms
    import getwaveform_iris
    from obspy.clients.fdsn import Client
    from obspy.core.event import Event, Origin, Magnitude
    if not user and not password:
        client = Client("IRIS")
    else:
        client = Client("IRIS",user=user,password=password)
    # will only work for events in the 'IRIS' catalog
    # (future: for Alaska events, read the AEC catalog)
    if use_catalog==1:
        cat = client.get_events(starttime = otime - sec_before_after_event,\
                                endtime = otime + sec_before_after_event)
        ev = cat[0]
    else:
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
        shutil.rmtree('RAW')
    if overwrite_ddir and os.path.exists(ddir):
        shutil.rmtree(ddir)

    # Extract waveforms, IRIS
    getwaveform_iris.run_get_waveform(c = client, event = ev, 
                                      min_dist = min_dist, max_dist = max_dist, 
                                      before = tbefore_sec, after = tafter_sec, 
                                      network = network, station = station, channel = channel, 
                                      resample_freq = resample_freq, ifrotate = rotate,
                                      ifCapInp = output_cap_weight_file, 
                                      ifRemoveResponse = remove_response,
                                      ifDetrend = detrend, ifDemean = demean, 
                                      ifEvInfo = output_event_info, 
                                      scale_factor = scale_factor, pre_filt = pre_filt, icreateNull=icreateNull,
                                      ifFilter = ifFilter, fmin = fmin, fmax = fmax, filt_type = filt_type, 
                                      zerophase = zerophase, corners = corners)

if idb == 3:
    import llnl_db_client
    import getwaveform_llnl
    client = llnl_db_client.LLNLDBClient(
            "/store/raw/LLNL/UCRL-MI-222502/westernus.wfdisc")

    # get event time and event ID
    cat = client.get_catalog()
    mintime_str = "time > %s" % (otime - sec_before_after_event)
    maxtime_str = "time < %s" % (otime + sec_before_after_event)
    print(mintime_str, maxtime_str)
    ev = cat.filter(mintime_str, maxtime_str)[0]
    print(ev)

    # Delete existing data directory 
    eid = util_helpers.otime2eid(ev.origins[0].time)
    ddir = './'+ eid
    if overwrite_ddir and os.path.exists(ddir):
        print("WARNING. Deleting data directory! (already exists)")
        shutil.rmtree(ddir)

    # Extract waveforms, LLNL
    getwaveform_llnl.run_get_waveform(llnl_db_client = client, event = ev, 
                                      min_dist = min_dist, max_dist = max_dist, 
                                      before = tbefore_sec, after = tafter_sec, 
                                      network = network, station = station, channel = channel, 
                                      resample_freq = resample_freq, ifrotate = rotate,
                                      ifCapInp = output_cap_weight_file, 
                                      ifRemoveResponse = remove_response,
                                      ifDetrend = detrend, ifDemean = demean, 
                                      ifEvInfo = output_event_info, 
                                      scale_factor = scale_factor, pre_filt = pre_filt)
