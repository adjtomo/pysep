import obspy
import read_event_obspy_file as reof
from getwaveform_new import *

def get_ev_info(ev_info,iex):
# ===============================================================
# EXAMPLE TEMPLATE -- USE THIS ONLY FOR QUICK TESTING
# This is a template for testing before creating an example.
# We can delete this if it creates more issues.
    if iex == 0:
        ev_info.otime = obspy.UTCDateTime("2009-04-07T20:12:55")
        ev_info.min_dist = 0 
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300
        ev_info.network = 'II,IU'
        ev_info.station = 'KDAK,COLA'
        ev_info.channel = '*'
        
#=================================================================================
# CATEGORY A EXAMPLES: simple test cases (including current and previous bugs)
#=================================================================================

# ERROR EXAMPLE [obspy]
# PROBLEM: No waveforms are returned -- perhaps related to the tbefore_sec request
# ERROR MESSAGE: ValueError: The length of the input vector x must be at least padlen, which is 39.
# SOLUTION: iex = 2 fails because the input data is much to short. The input is a trace with 38 sample but sampled at 50 Hz and you want to resample to 10 Hz. The IIR filter coefficients it calculates are just too long. I'm working on getting such a filter into ObsPy which should then have nicer error messages. I doubt its possible to meaningfully filter such short data with a very sharpy filter. For now: Just add a QA step that makes sure that at least a certain number of samples enter the decimation routine.
    if iex == 1:
        ev_info.otime = obspy.UTCDateTime("2016-01-24T10:30:29.557")
        ev_info.min_dist = 0 
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 42          # Crashes
        ev_info.tbefore_sec = 41         # works fine
        ev_info.tafter_sec = 600
        ev_info.network = 'AV'       # Crashes when -  tbefore_sec = 42; Works fine when - tbefore_sec = 41
        ev_info.channel = 'BH?'
        ev_info.rlat = 64.7716
        ev_info.rlon = -149.1465
        ev_info.rtime = obspy.UTCDateTime("2009-04-07T20:20:00")
        
# ERROR EXAMPLE [obspy]
# PROBLEM: If a particular network is requested (either explicitly or within *), no waveforms are returned.
# ERROR MESSAGE: NotImplementedError: ResponseListResponseStage not yet implemented due to missing example data. Please contact the developers with a test data set (waveforms and StationXML metadata).
# KLUDGE: list all networks explicitly, except IM
# SOLUTION: https://github.com/obspy/obspy/issues/1514
#           If you want to use it right now you'd have to use that branch 
#           - it will be a while before we release a new ObsPy version.
    if iex == 2:
        ev_info.otime = obspy.UTCDateTime("2009-04-07T20:12:55")
        ev_info.min_dist = 300 
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300
        # ev_info.network = '*'       # crashes
        # ev_info.network = 'AK,IM'   # crashes (Error because of the IM network)
        # ev_info.network = '*,-IM'   # stalls indefinitely (syntax works for stations but not networks)
        ev_info.network = 'AK'       # works fine 
        ev_info.channel = 'BH?'
        
# SALMON example (restricted data from IRIS)
    if iex == 4:
        ev_info.otime = obspy.UTCDateTime("2016-01-24T10:30:29.557")
        ev_info.min_dist = 0 
        ev_info.max_dist = 300
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = 'AK,ZE'   # ZE waveforms not returned (only ZE.MPEN)
        ev_info.channel = 'BH?,HH?'
        
# ROTATION example for components 1,2,Z
#    All rotations should be based on the azimuth of the sensor
#    (CMPAZ header in sac), which, combined with the station-source backazimuth,
#    will provide the rotation to radial and transverse components.
#    The rotation should NOT depend on channel names like E or N.
#    (Even for GSN stations the E and N do not point exactly E and N.)
    if iex == 5:
        ev_info.otime = obspy.UTCDateTime("2016-01-24T10:30:29.557")
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        # ev_info.network = 'II,AK'
        ev_info.station = 'KDAK,SWD'
        ev_info.channel = 'BH?'
        
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
        ev_info.idb = 3              # LLNL database
        # resample_freq = 0    # no resampling and no cutting
        
        # TARGET origin time (8.031 s from actual origin time)
        ev_info.otime = obspy.UTCDateTime("1991-09-14T19:00:08.031Z")
        # evid = 635527        # Hoya event id in LLNL database
        
        ev_info.min_dist = 0 
        ev_info.max_dist = 1200
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = '*'        # note that the client will look for BK stations in the list
        ev_info.channel = 'BH?'      # ALL channels from LLNL are returned
        # ev_info.scale_factor = 10.0**2  # original
        ev_info.scale_factor = 2e-1     # Hoya  
        ev_info.overwrite_ddir = 0

# same as iex=6 but for the IRIS database
# GOAL: For LLNL events, we do NOT want to use the IRIS source parameters:
#       origin time, hypocenter, magnitude.
    if iex == 7:
        ev_info.idb = 1            # IRIS database
        # ev_info.resample_freq = 0  # no resampling
        # ev_info.otime = obspy.UTCDateTime("1991-09-14T19:00:00.000Z")   # Hoya actual
        ev_info.otime = obspy.UTCDateTime("1991-09-14T19:00:08.031Z")   # Hoya target
        ev_info.min_dist = 0 
        ev_info.max_dist = 1200
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        # needs to be run TWICE to get BK stations and IRIS stations
        # ev_info.network = 'BK'        # BK will go to NCEDC
        ev_info.network = '*'         # * will give all at IRIS DMC
        ev_info.channel = 'BH?,LH?' 
        ev_info.overwrite_ddir = 0

# problem 1: some stations return only vertical component. our tools crash in this case.
# problem 2: short waveforms. padding by zeros adds sharp changes in the data
# problem 3: waveform contains NAN or INF. this crashes detrend
# solution 1: (ongoing)
# solution 2: disable zero-padding
# solution 3: (ongoing) print error (consider removing trace?)
    if iex==8:
        ev_info.idb = 1            # IRIS database
        ev_info.otime = obspy.UTCDateTime("1982-08-05T14:00:00.000000Z")
        ev_info.min_dist = 0 
        ev_info.max_dist = 1200
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = '*'
        ev_info.channel = 'BH?,LH?,EH?' 
        
# 1. Waveform extraction for user defined event info
# 2. Subset of stations for quickly testing data gaps and padding tests
    if iex == 9:
        ev_info.otime = obspy.UTCDateTime("2009-04-07T20:12:55.351")
        ev_info.elat = 61.4542
        ev_info.elon =-149.7428
        ev_info.edep = 33033.6  # in meters
        ev_info.emag = 4.6
        ev_info.min_dist = 0 
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300
        ev_info.network = 'AK,AT,AV,CN,II,IU,US,XM,XV,XZ,YV'  # note: cannot use '*' because of IM
        ev_info.channel = 'BH?'
        ev_info.network = 'AV'
        ev_info.station = 'SPBG,KABU'                      # For testing data gaps 
        ev_info.use_catalog = 0                            # To get (lat,lon, etime, dep, mag) from some catalog = 1 OR use defined = 0 (see iex=9)

# Error from util_write_cap.py
#   util_write_cap.py, line 437, in add_sac_metadata
#    and pick.waveform_id.channel_code[2].upper() == 'Z' and pick.waveform_id.location_code == tr.stats.location and pick.phase_hint == 'Pn'):
#        IndexError: string index out of range
#
# This error happens when allowing data already rotated (VRT). see getwaveform_llnl.py Line 94
    if iex == 10:
        ev_info.idb = 3              # LLNL database
        ev_info.otime = obspy.UTCDateTime("1982-08-05T14:00:00.000000Z")
        ev_info.min_dist = 0 
        ev_info.max_dist = 1200
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = '*'        # note that the client will look for BK stations in the list
        ev_info.channel = 'BH?'      # ALL channels from LLNL are returned regardless
        ev_info.scale_factor = 10.0**2  # original
        ev_info.overwrite_ddir = 0
        ev_info.resample_freq = 20.0 
        
# error applying rotation to station LL.KNB. See util_write_cap.py Line 161
# (Disable the try block and only allow Line 161 to create this problem)
# 
#   File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/signal/rotate.py", line 257, in rotate2zne
#    x, y, z = np.dot(_t, [data_1, data_2, data_3])
#    ValueError: setting an array element with a sequence.
    if iex == 11:
        ev_info.idb = 1 
        ev_info.otime = obspy.UTCDateTime("1995-07-31T12:34:46.860000Z")
        ev_info.min_dist = 0 
        ev_info.max_dist = 1200
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = '*'        # note that the client will look for BK stations in the list
        ev_info.channel = 'BH?'      # ALL channels from LLNL are returned regardless
        ev_info.scale_factor = 10.0**2  # original
        ev_info.overwrite_ddir = 0
        ev_info.resample_freq = 20.0 

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
        ev_info.idb = 1 
        ev_info.otime = obspy.UTCDateTime("1999-01-23T03:00:33.200000Z") # frenchman flat 1
        # ev_info.otime = obspy.UTCDateTime("1999-01-27T10:44:23.310000Z") # frenchman flat 2
        # ev_info.otime = obspy.UTCDateTime("2000-01-30T14:46:51.310000Z") # trona mine 2
        # ev_info.otime = obspy.UTCDateTime("2002-06-14T12:40:44.450000Z") # little skull
        ev_info.min_dist = 0 
        ev_info.max_dist = 1200
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = '*'        # note that the client will look for BK stations in the list
        ev_info.channel = 'BH?'      # ALL channels from LLNL are returned regardless
        ev_info.scale_factor = 10.0**2  # original
        ev_info.overwrite_ddir = 0
        ev_info.resample_freq = 20.0 

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
        ev_info.idb = 1 
        ev_info.otime = obspy.UTCDateTime("1997-06-14T19:48:19.930000Z") # indian springs
        # ev_info.otime = obspy.UTCDateTime("1996-09-05T08:16:56.090000Z") # amargosa
        ev_info.min_dist = 0 
        ev_info.max_dist = 1200
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = '*'        # note that the client will look for BK stations in the list
        ev_info.channel = 'BH?'      # ALL channels from LLNL are returned regardless
        ev_info.scale_factor = 10.0**2  # original
        ev_info.overwrite_ddir = 0
        ev_info.resample_freq = 20.0 

# time > 2007-01-24T11:30:06.100000Z time < 2007-01-24T11:30:26.100000Z
#Traceback (most recent call last):
#  File "run_getwaveform.py", line 478, in <module>
#    ev = cat.filter(mintime_str, maxtime_str)[0]
#  File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/core/event/catalog.py", line 162, in __getitem__
#    return ev_info.events.__getitem__(index)
#IndexError: list index out of range
    if iex == 14:
        ev_info.idb = 3 
        ev_info.otime = obspy.UTCDateTime("2007-01-24T11:30:16.100000Z") # ralston
        ev_info.min_dist = 0 
        ev_info.max_dist = 1200
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = '*'        # note that the client will look for BK stations in the list
        ev_info.channel = 'BH?'      # ALL channels from LLNL are returned regardless
        ev_info.scale_factor = 10.0**2  # original
        ev_info.overwrite_ddir = 0
        ev_info.resample_freq = 20.0 

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
        ev_info.idb = 1 
        ev_info.otime = obspy.UTCDateTime("1995-07-31T12:34:46.860000Z") # timber mountain
        ev_info.min_dist = 0 
        ev_info.max_dist = 1200
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = '*'        # note that the client will look for BK stations in the list
        ev_info.channel = 'BH?'      # ALL channels from LLNL are returned regardless
        ev_info.scale_factor = 10.0**2  # original
        ev_info.overwrite_ddir = 0
        ev_info.resample_freq = 20.0 

# Iniskin earthquake - all channels at IU.COLA and II.KDAK
# note different channels and location codes for strong motion:
#       IU.COLA.20.HNZ.sac       http://ds.iris.edu/mda/IU/COLA
#       II.KDAK.00.ENZ.sac       http://ds.iris.edu/mda/II/KDAK
    if iex == 16:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2016-01-24T10:30:37.400") 
        ev_info.elat = 59.75
        ev_info.elon = -153.27
        ev_info.edep = 110700
        ev_info.emag = 7.1
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 800
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = 'IU,II'
        # ev_info.channel = '?H?,?N?'
        ev_info.channel = 'HH?,BH?,BN?,HN?,EN?'
        ev_info.resample_freq = 0        # no resampling
        ev_info.scale_factor = 1         # no scale factor

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
        ev_info.idb = 3            # LLNL
        ev_info.otime = obspy.UTCDateTime("1989-06-27T15:30:00.02")
        ev_info.elat = 37.275
        ev_info.elon = -116.354
        ev_info.edep = 640
        ev_info.emag = 4.90
        ev_info.min_dist = 0 
        ev_info.max_dist = 1200
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = '*' 
        ev_info.channel = 'BH?,LH?' 
        ev_info.overwrite_ddir = 0

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
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2015-12-02T10:05:25.798") 
        ev_info.elat = 61.70
        ev_info.elon = -147.26
        ev_info.edep = 36590
        ev_info.emag = 4.50
        # subset of stations
        ev_info.min_dist = 300
        ev_info.max_dist = 400
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        # ev_info.network = 'AK,AT,AV,CN,II,IU,US,XM,XV,XZ,YV,ZE'
        ev_info.network = 'IU'
        ev_info.channel = 'HH?,BH?'
        # ev_info.resample_freq = 0        # no resampling -- THIS WORKS
        ev_info.resample_freq = 50       # THIS FAILS BUT ONLY AFTER REPEATING
        ev_info.scale_factor = 1         # no scale factor

# Test case for UVW
    if iex == 19:
        ev_info.idb = 1
        ev_info.use_catalog = 0  
        ev_info.otime = obspy.UTCDateTime("2016-12-08T10:16:00")
        ev_info.elat = 64.2380
        ev_info.elon = -150.0581
        ev_info.edep = 18507
        ev_info.emag = 4.60
        ev_info.station = 'F3TN'
            
        ev_info.min_dist = 0
        ev_info.max_dist = 300
        ev_info.tbefore_sec = 0
        ev_info.tafter_sec = 600
        ev_info.network = 'XV'
        ev_info.channel = 'HH?'
        ev_info.scale_factor = 1
        ev_info.resample_freq = 0
        ev_info.detrend = False
        ev_info.demean = False
        ev_info.taper = False
        ev_info.ipre_filt = 1

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
#     return ev_info.traces.__getitem__(index)
# IndexError: list index out of range
    if iex == 20:
        ev_info.idb = 3            # LLNL
        ev_info.otime = obspy.UTCDateTime("1990-06-13T16:00:00.09")
        ev_info.elat = 37.262
        ev_info.elon = -116.421
        ev_info.edep = 674
        ev_info.emag = 5.34
        ev_info.min_dist = 0 
        ev_info.max_dist = 1200
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = '*' 
        ev_info.channel = 'BH?,LH?' 
        ev_info.overwrite_ddir = 0

#=================================================================================
# CATEGORY B EXAMPLES: important events
# 1XX: southern Alaska
# 2XX: central Alaska
# 9XX: other
#=================================================================================
# SilwalTape2016 example event (Anchorage)
    if iex == 100:
        ev_info.otime = obspy.UTCDateTime("2009-04-07T20:12:55.351")
        ev_info.min_dist = 0 
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300
        ev_info.network = 'AK,AT,AV,CN,II,IU,US,XM,XV,XZ,YV'  # note: cannot use '*' because of IM
        ev_info.channel = 'BH?'
        ev_info.use_catalog = 0 
        ev_info.elat = 61.45420
        ev_info.elon = -149.7428
        ev_info.edep = 33033.60
        ev_info.emag = 4.6
        ev_info.resample_freq = 50


# Iniskin earthquake
# NOTE: must enter username and password above to get SALMON (ZE) stations
    if iex == 101:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2016-01-24T10:30:37.400") 
        ev_info.elat = 59.75
        ev_info.elon = -153.27
        ev_info.edep = 110700
        ev_info.emag = 7.1
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 800
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US'  # IM will probably crash it
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 0        # no resampling
        ev_info.scale_factor = 1         # no scale factor
        
        # parameters for examining step response (causal low-pass filter on raw waveforms)
        # delete AUQ, SPCP, SPBG
        ev_info.rotateRTZ = False
        ev_info.ifFilter = True
        ev_info.filter_type = 'lowpass'
        ev_info.f1 = 1/4
        ev_info.zerophase = False
        ev_info.remove_response = False
        ev_info.ipre_filt = 0
        ev_info.demean = False
        ev_info.detrend = False

# Iniskin earthquake - all strong motion
    if iex == 103:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2016-01-24T10:30:37.400") 
        ev_info.elat = 59.75
        ev_info.elon = -153.27
        ev_info.edep = 110700
        ev_info.emag = 7.1
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 800
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 600
        ev_info.network = '*'
        ev_info.channel = 'BN?,HN?,EN?'
        ev_info.resample_freq = 0        # no resampling
        ev_info.scale_factor = 1         # no scale factor

# Totschunda fault
    if iex == 104:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2017-02-13T07:17:12.000") 
        ev_info.elat = 62.5154
        ev_info.elon = -142.7485
        ev_info.edep = 7564
        ev_info.emag = 5.3
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 50
        ev_info.tafter_sec = 300
        ev_info.network = 'AV'
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor

# Klukwan earthquakes
    if iex == 105:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2017-05-01T12:31:53.000") 
        ev_info.elat = 59.8522
        ev_info.elon = -136.6618
        ev_info.edep = 4000
        ev_info.emag = 6.2
        ev_info.otime = obspy.UTCDateTime("2017-05-01T14:18:14.000") 
        ev_info.elat = 59.8184
        ev_info.elon = -136.7163
        ev_info.edep = 1000
        ev_info.emag = 6.0
        
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 500
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor

# Kantishna earthquakes
    if iex == 106:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2017-04-29T11:15:48.000") 
        ev_info.elat = 63.1296
        ev_info.elon = -151.1517
        ev_info.edep = 11000
        ev_info.emag = 5.2
        ev_info.otime = obspy.UTCDateTime("2017-01-31T09:38:37.000") 
        ev_info.elat = 63.0817
        ev_info.elon = -150.9427
        ev_info.edep = 135000
        ev_info.emag = 5.2
        
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 500
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor

# Cook Inlet earthquake
    if iex == 107:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2017-05-07T04:25:19.000") 
        ev_info.elat = 60.1945
        ev_info.elon = -151.6743
        ev_info.edep = 64000
        ev_info.emag = 5.2
    
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 500
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor
            
# MFFZ earthquakes for investigating the step response
# LISTED IN CHRONOLOGICAL ORDER
    if iex == 200:
        ev_info.idb = 1
        ev_info.use_catalog = 0
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2014-08-31T03:06:57.111")
        ev_info.elat = 65.1526
        ev_info.elon = -149.0398
        ev_info.edep = 16614.7
        ev_info.emag = 5.20
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2014-10-21T00:36:58.333")
        ev_info.elat = 65.1489
        ev_info.elon = -149.0413
        ev_info.edep = 13134.8
        ev_info.emag = 4.90
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2014-10-23T16:30:23.968")
        ev_info.elat = 65.1644
        ev_info.elon = -149.0523
        ev_info.edep = 20066.5
        ev_info.emag = 5.00
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2015-09-12T03:25:12.711")
        ev_info.elat = 65.1207
        ev_info.elon = -148.6646
        ev_info.edep = 15556.8
        ev_info.emag = 2.63
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2015-10-31T02:56:35.572")
        ev_info.elat = 64.4285
        ev_info.elon = -149.6969
        ev_info.edep = 23852.1
        ev_info.emag = 3.47
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2016-01-14T19:04:10.727")
        ev_info.elat = 64.6827
        ev_info.elon = -149.2479
        ev_info.edep = 22663.7
        ev_info.emag = 3.80
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2016-11-06T9:29:10.579")
        ev_info.elat = 64.1639
        ev_info.elon = -150.0626
        ev_info.edep = 23190.0
        ev_info.emag = 4.00
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2016-12-08T10:18:13.868")
        ev_info.elat = 64.1937
        ev_info.elon = -150.0376
        ev_info.edep = 24522.1
        ev_info.emag = 4.30
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2017-05-08T05:09:02.000") 
        ev_info.elat = 65.2643
        ev_info.elon = -146.922
        ev_info.edep = 9000   # AEC/Vipul
        ev_info.emag = 3.60   # Vipul
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2017-06-28T12:58:52.000") 
        ev_info.elat = 64.7569
        ev_info.elon = -148.8883
        ev_info.edep = 18000  # Vipul
        ev_info.emag = 3.50   # Vipul
        # -------------------------------------------------
        # VIPUL: WHAT ARE THESE? OTHER EVENTS?
        # ev_info.otime = obspy.UTCDateTime("2015-03-30T12:33:19.000")
        # ev_info.otime = obspy.UTCDateTime("2015-10-20T19:14:16.000")
        # ev_info.otime = obspy.UTCDateTime("2011-12-21T16:28:41.000")
        # -------------------------------------------------
        
        ev_info.min_dist = 0
        ev_info.max_dist = 300
        ev_info.tbefore_sec = 200
        ev_info.tafter_sec = 600
        ev_info.network = 'AV,CN,AT,TA,AK,XV,II,IU,US'
        # ev_info.network = 'XV,AK'
        ev_info.channel = 'BH?,HH?'
        ev_info.scale_factor = 1
        ev_info.resample_freq = 0
        # for CAP
        # ev_info.scale_factor = 10.0**2
        # ev_info.resample_freq = 50 
        
        # to investigate step response
        ev_info.rotateRTZ = True
        ev_info.rotateUVW = True
        ev_info.ifFilter = True
        ev_info.zerophase = False    # causal
        # filter_type = 'lowpass'
        ev_info.filter_type = 'bandpass'
        ev_info.f1 = 1/100  # fmin
        ev_info.f2 = 1/10  # fmax
        ev_info.corners = 4
        # ev_info.remove_response = False
        ev_info.remove_response = True
        ev_info.ipre_filt = 1
        ev_info.demean = True
        ev_info.detrend = True
        ev_info.output_cap_weight_file = False
        # ev_info.outformat = 'DISP'
        ev_info.ifsave_sacpaz = True
        ev_info.taper = 0.2
        
# NENNUC event (from Steve)
    if iex == 201:
        ev_info.idb = 1
        ev_info.use_catalog = 0          # manually enter AEC catalog parameters
        ev_info.otime = obspy.UTCDateTime("2016-01-14T19:04:10.727")
        ev_info.elat = 64.6827
        ev_info.elon = -149.2479
        ev_info.edep = 22663.7
        ev_info.emag = 3.8
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2015-09-12T03:25:12.711")
        ev_info.elat = 65.1207
        ev_info.elon = -148.6646
        ev_info.edep = 15556.8
        ev_info.emag = 2.6
        # -------------------------------------------------
        ev_info.otime = obspy.UTCDateTime("2013-03-12T07:39:50.214")
        ev_info.elat = 64.7161
        ev_info.elon = -148.9505
        ev_info.edep = 20000
        ev_info.emag = 2.1
            
        # For CAP
        ev_info.min_dist = 0
        ev_info.max_dist = 200
        ev_info.tbefore_sec = 500
        ev_info.tafter_sec = 2000
        ev_info.taper = 0.1
        ev_info.network = 'AV,CN,AT,TA,AK,XV,II,IU,US'
        # ev_info.network = 'AK,TA,II,IU,US'
        # ev_info.station = "-RAG"  # RAG has choppy data; gives error: Exception: Can't merge traces with same ids but differing sampling rates!
        ev_info.channel = 'BH?,HH?'
        ev_info.pre_filt = (0.0025, 0.003, 10.0, 15.0) # BH default
        ev_info.ipre_filt = 1

        # to investigate step response
        # ev_info.ev_info.tbefore_sec = 2000
        # ev_info.tafter_sec = 2000
        # ev_info.resample_freq = 0 
        # ev_info.scale_factor = 1
        # ev_info.ifFilter = True
        # ev_info.zerophase = False
        # ev_info.filter_type = 'bandpass'
        # ev_info.f1 = 1/100
        # ev_info.f2 = 1/20
        # ev_info.remove_response = True
        # ev_info.ipre_filt = 0
        # ev_info.demean = True
        # ev_info.detrend = True
        
# Kyle Nenana basin earthquakes for basin amplification
    if iex == 210:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # AEC source parameters
        ev_info.otime = obspy.UTCDateTime("2015-11-20T10:53:48.168") 
        ev_info.elat = 64.6210
        ev_info.elon = -149.4024
        ev_info.edep = 17113.4
        ev_info.emag = 2.67
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 200
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 400
        ev_info.network = 'AK,AT,II,IU,US,XM,XV,XZ,TA'  # no CN,AV,YV,ZE
        ev_info.channel = 'BH?,HH?'
        # ev_info.resample_freq = 0        # no resampling
        # ev_info.scale_factor = 1         # no scale factor
        # For CAP moment tensor
        ev_info.resample_freq = 50         # same as greens function 
        ev_info.scale_factor = 100         # change from m/s to cm/s
        # ev_info.ipre_filt = 0
        ev_info.remove_response = True
        # ev_info.demean = False
        # ev_info.detrend = False

# gets waveforms from M 8.3 Chile event with stations centered in Minto
    if iex == 211:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0
        ev_info.otime = obspy.UTCDateTime("2015-09-16T22:54:33.000")
        ev_info.elat = -31.5695
        ev_info.elon = -71.6543
        ev_info.edep = 22400.0
        ev_info.emag = 8.30
        ev_info.rlat = 64.7716
        ev_info.rlon = -149.1465
        ev_info.rtime = obspy.UTCDateTime("2015-09-16T23:09:15.000")
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 200
        ev_info.min_dist = 0
        ev_info.max_dist = 100
        ev_info.network = 'AK,AT,II,IU,US,XM,XV,XZ,TA'  # no CN,AV,YV,ZE
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 0        
        ev_info.scale_factor = 1        
        ev_info.remove_response = True

# gets a ???
    if iex == 212:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0
        ev_info.otime = obspy.UTCDateTime("2016-06-06T00:00:00.000")
        ev_info.elat = 64.6130
        ev_info.elon = -149.0992
        ev_info.edep = 0
        ev_info.emag = 0.00
        # ev_info.rlat = 64.7716
        # ev_info.rlon = -149.1465
        # ev_info.rtime = obspy.UTCDateTime("2015-09-16T23:09:15.000")
        ev_info.tbefore_sec = 0
        ev_info.tafter_sec = 3600
        ev_info.min_dist = 0
        ev_info.max_dist = 100
        ev_info.network = 'XV,AK,TA'  # no CN,AV,YV,ZE
        ev_info.channel = 'HH?,BH?'
        ev_info.resample_freq = 0        
        ev_info.scale_factor = 1        
        ev_info.remove_response = True
        ev_info.rotateRTZ = False
        # ev_info.pre_filt=(f0*0.001, f1*0.001, f2*1000, f3*1000)
        # ev_info.ipre_filt = 2
        ev_info.ipre_filt = 0

# Chatnika earthquake
    if iex == 213:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2017-05-08T05:09:02.000") 
        ev_info.elat = 65.2643
        ev_info.elon = -146.922
        ev_info.edep = 9000
        ev_info.emag = 3.8
        
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 500
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor
            
# NE Nenana earthquake
    if iex == 214:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2017-06-28T12:58:52") 
        ev_info.elat = 64.7569
        ev_info.elon = -148.8883
        ev_info.edep = 18000
        ev_info.emag = 3.5
        
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 500
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 500
        ev_info.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
        ev_info.channel = 'BH?,HH?'
        ev_info.resample_freq = 50        
        ev_info.scale_factor = 100 

# Loop over multiple event       
    if iex == 215:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        
        events_file = "/home/ksmith/REPOSITORIES/manuscripts/kyle/papers/basinamp/data/basinamp_obspy.txt"
        eids,otimes,elons,elats,edeps,emags = reof.read_events_obspy_file(events_file)

        ev_info_list = []
        for xx in range(0,3):
            ev_info_temp = ev_info.copy()
            ev_info_temp.otime = obspy.UTCDateTime(otimes[xx])
            ev_info_temp.elat = elats[xx]
            ev_info_temp.elon = elons[xx]
            ev_info_temp.edep = edeps[xx]
            ev_info_temp.emag = emags[xx]
            ev_info_temp.eid = eids[xx]
            
            # subset of stations
            ev_info_temp.min_dist = 0
            ev_info_temp.max_dist = 500
            ev_info_temp.tbefore_sec = 100
            ev_info_temp.tafter_sec = 500
            ev_info_temp.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
            ev_info_temp.channel = 'BH?,HH?'
            ev_info_temp.resample_freq = 50        
            ev_info_temp.scale_factor = 100 

            # append getwaveform objects
            ev_info_list.append(ev_info_temp)
        
        # always return ev_info
        ev_info = ev_info_list
#----------------------------------------------------------

#-----------------------------------------------------------
# Event from HutchisonGhosh2016
# Crashes when using all networks [i.e. network = '*']
    if iex == 900:
        ev_info.idb = 1
        ev_info.overwrite_ddir = 1       # delete data dir if it exists
        ev_info.use_catalog = 0          # do not use event catalog for source parameters
        # GCMT source parameters
        # the otime is the centroid time and accounts for tshift
        ev_info.otime = obspy.UTCDateTime("2014-12-12T00:13:17") 
        ev_info.elat = 49.10
        ev_info.elon = -124.05
        ev_info.edep = 60000
        ev_info.emag = 4.10
        # subset of stations
        ev_info.min_dist = 0
        ev_info.max_dist = 250
        ev_info.tbefore_sec = 100
        ev_info.tafter_sec = 300
        ev_info.network = '*'
        ev_info.channel = 'LH?,BH?'
        ev_info.resample_freq = 50        # no resampling
        ev_info.scale_factor = 100         # no scale factor
        
        # For plotting filtered waveforms
        ev_info.tbefore_sec = 500
        ev_info.tafter_sec = 500
        ev_info.resample_freq = 0 
        ev_info.scale_factor = 1
        ev_info.ifFilter = True
        ev_info.zerophase = True
        ev_info.filter_type = 'bandpass'
        ev_info.f1 = 1/50
        ev_info.f2 = 1/20
        ev_info.remove_response = True
        ev_info.ipre_filt = 2
        ev_info.demean = True
        ev_info.detrend = True
        ev_info.taper = True

    return(ev_info)
#=================================================================================
# END EXAMPLES
#=================================================================================
