import obspy

class getwaveform_input:
    def __init__(self):
        # DEFAULT PARAMETERS
        
        # For CAP
        self.resample_freq = 40           # 0 causes errors. Use resample_TF instead
        self.scale_factor = 10**2         # for CAP use 10**2  (to convert m/s to cm/s)

        # event parameters
        self.use_catalog = 1              # use an existing catalog (=1) or specify your own event parameters (see iex=9)
        self.sec_before_after_event = 10  # time window to search for a target event in a catalog
        self.min_dist = 0 
        self.max_dist = 20000
        self.tbefore_sec = 100
        self.tafter_sec = 300

        # station parameters
        self.network = '*'                # all networks
        self.station = '*,-PURD,-NV33,-GPO'  # all stations
        self.channel = '*'                # all channels     
        self.overwrite_ddir = 1           # 1 = delete data directory if it already exists
        self.icreateNull = 1              # create Null traces so that rotation can work (obsby stream.rotate require 3 traces)
        
        # Filter parameters
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
        self.filter_type = 'bandpass'
        # f1 should consider the requested length of the time series
        # f2 should consider the sampling rate for the desired channels
        self.f1 = 1/40 # fmin - highpass will keep frequencies larger than fmin
        self.f2 = 1/10 # fmax - lowpass will keep frequencies lower than fmax
        self.zerophase = True             # = False (causal), = True (acausal)
        # 4 pole filter is more sharper at the edges than 2 pole
        self.corners = 4                  # Is corner in Obspy same as Pole in SAC?
        
        # For tapering down the pre-filter
        self.f0 = 0.5*self.f1
        self.f3 = 2.0*self.f2
        self.pre_filt=(self.f0, self.f1, self.f2, self.f3)    # applies for ipre_filt = 2 only
        # pre_filt = (0.005, 0.006, 10.0, 15.0) # BH default

        # Pre-processing (manily for CAP)
        self.rotateRTZ = True
        self.rotateUVW = False   # This option works only if 'rotateRTZ = True'
        self.output_cap_weight_file = True
        self.detrend = True
        self.demean = True
        self.taper = False       # this could also be a fraction between 0 and 1 (fraction to be tapered from both sides)
        self.output_event_info = True
        self.outformat = 'VEL'            # Intrument response removed waveforms could be saved as 'VEL' 'DISP' 'ACC'
        self.ifsave_sacpaz = False        # save sac pole zero (needed as input for MouseTrap module)

        # for CAP all waveforms need to have the same sample rate
        self.resample_TF = True
        
        # Perhaps you want to set ipre_filt = 0 to prevent extra filtering
        self.ifFilter = False 
        # pre-filter for deconvolution
        # https://ds.iris.edu/files/sac-manual/commands/transfer.html
        # Pre-filter will not be applied if remove_response = False 
        self.ipre_filt = 1                # =0 No pre_filter
                                     # =1 default pre_filter (see getwaveform_iris.py)
                                     # =2 user-defined pre_filter (use this if you are using bandpass filter)

        self.remove_response = True
        self.iplot_response = False
        self.ifplot_spectrogram = False

    def get_extraction_info(self,iex):
        # ===============================================================
        # EXAMPLE TEMPLATE -- DO NOT USE
        # This is a template for testing before creating an example.
        # We can delete this if it creates more issues.
        if iex == 0:
            self.otime = obspy.UTCDateTime("2009-04-07T20:12:55")
            self.min_dist = 0 
            self.max_dist = 500
            self.tbefore_sec = 100
            self.tafter_sec = 300
            self.network = 'II,IU'
            self.station = 'KDAK,COLA'
            self.channel = '*'
            
#=================================================================================
# CATEGORY A EXAMPLES: simple test cases (including current and previous bugs)
#=================================================================================

# ERROR EXAMPLE [obspy]
# PROBLEM: No waveforms are returned -- perhaps related to the tbefore_sec request
# ERROR MESSAGE: ValueError: The length of the input vector x must be at least padlen, which is 39.
# SOLUTION: iex = 2 fails because the input data is much to short. The input is a trace with 38 sample but sampled at 50 Hz and you want to resample to 10 Hz. The IIR filter coefficients it calculates are just too long. I'm working on getting such a filter into ObsPy which should then have nicer error messages. I doubt its possible to meaningfully filter such short data with a very sharpy filter. For now: Just add a QA step that makes sure that at least a certain number of samples enter the decimation routine.
        if iex == 1:
            self.otime = obspy.UTCDateTime("2016-01-24T10:30:29.557")
            self.min_dist = 0 
            self.max_dist = 500
            self.tbefore_sec = 42          # Crashes
            self.tbefore_sec = 41         # works fine
            self.tafter_sec = 600
            self.network = 'AV'       # Crashes when -  tbefore_sec = 42; Works fine when - tbefore_sec = 41
            self.channel = 'BH?'
            self.rlat = 64.7716
            self.rlon = -149.1465
            self.rtime = obspy.UTCDateTime("2009-04-07T20:20:00")
            
            # ERROR EXAMPLE [obspy]
            # PROBLEM: If a particular network is requested (either explicitly or within *), no waveforms are returned.
            # ERROR MESSAGE: NotImplementedError: ResponseListResponseStage not yet implemented due to missing example data. Please contact the developers with a test data set (waveforms and StationXML metadata).
            # KLUDGE: list all networks explicitly, except IM
            # SOLUTION: https://github.com/obspy/obspy/issues/1514
            #           If you want to use it right now you'd have to use that branch 
            #           - it will be a while before we release a new ObsPy version.
        if iex == 2:
            self.otime = obspy.UTCDateTime("2009-04-07T20:12:55")
            self.min_dist = 300 
            self.max_dist = 500
            self.tbefore_sec = 100
            self.tafter_sec = 300
            # self.network = '*'       # crashes
            # self.network = 'AK,IM'   # crashes (Error because of the IM network)
            # self.network = '*,-IM'   # stalls indefinitely (syntax works for stations but not networks)
            self.network = 'AK'       # works fine 
            self.channel = 'BH?'
            
# SALMON example (restricted data from IRIS)
        if iex == 4:
            self.otime = obspy.UTCDateTime("2016-01-24T10:30:29.557")
            self.min_dist = 0 
            self.max_dist = 300
            self.tbefore_sec = 100
            self.tafter_sec = 600
            self.network = 'AK,ZE'   # ZE waveforms not returned (only ZE.MPEN)
            self.channel = 'BH?,HH?'

# ROTATION example for components 1,2,Z
#    All rotations should be based on the azimuth of the sensor
#    (CMPAZ header in sac), which, combined with the station-source backazimuth,
#    will provide the rotation to radial and transverse components.
#    The rotation should NOT depend on channel names like E or N.
#    (Even for GSN stations the E and N do not point exactly E and N.)
        if iex == 5:
            self.otime = obspy.UTCDateTime("2016-01-24T10:30:29.557")
            self.tbefore_sec = 100
            self.tafter_sec = 600
            # self.network = 'II,AK'
            self.station = 'KDAK,SWD'
            self.channel = 'BH?'

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
            self.idb = 3              # LLNL database
            # resample_freq = 0    # no resampling and no cutting
            
            # TARGET origin time (8.031 s from actual origin time)
            self.otime = obspy.UTCDateTime("1991-09-14T19:00:08.031Z")
            # evid = 635527        # Hoya event id in LLNL database
            
            self.min_dist = 0 
            self.max_dist = 1200
            self.tbefore_sec = 100
            self.tafter_sec = 600
            self.network = '*'        # note that the client will look for BK stations in the list
            self.channel = 'BH?'      # ALL channels from LLNL are returned
            # self.scale_factor = 10.0**2  # original
            self.scale_factor = 2e-1     # Hoya  
            self.overwrite_ddir = 0

# same as iex=6 but for the IRIS database
# GOAL: For LLNL events, we do NOT want to use the IRIS source parameters:
#       origin time, hypocenter, magnitude.
        if iex == 7:
            self.idb = 1            # IRIS database
            # self.resample_freq = 0  # no resampling
            # self.otime = obspy.UTCDateTime("1991-09-14T19:00:00.000Z")   # Hoya actual
            self.otime = obspy.UTCDateTime("1991-09-14T19:00:08.031Z")   # Hoya target
            self.min_dist = 0 
            self.max_dist = 1200
            self.tbefore_sec = 100
            self.tafter_sec = 600
            # needs to be run TWICE to get BK stations and IRIS stations
            # self.network = 'BK'        # BK will go to NCEDC
            self.network = '*'         # * will give all at IRIS DMC
            self.channel = 'BH?,LH?' 
            self.overwrite_ddir = 0

# problem 1: some stations return only vertical component. our tools crash in this case.
# problem 2: short waveforms. padding by zeros adds sharp changes in the data
# problem 3: waveform contains NAN or INF. this crashes detrend
# solution 1: (ongoing)
# solution 2: disable zero-padding
# solution 3: (ongoing) print error (consider removing trace?)
        if iex==8:
            self.idb = 1            # IRIS database
            self.otime = obspy.UTCDateTime("1982-08-05T14:00:00.000000Z")
            self.min_dist = 0 
            self.max_dist = 1200
            self.tbefore_sec = 100
            self.tafter_sec = 600
            self.network = '*'
            self.channel = 'BH?,LH?,EH?' 

# 1. Waveform extraction for user defined event info
# 2. Subset of stations for quickly testing data gaps and padding tests
        if iex == 9:
            self.otime = obspy.UTCDateTime("2009-04-07T20:12:55.351")
            self.elat = 61.4542
            self.elon =-149.7428
            self.edep = 33033.6  # in meters
            self.emag = 4.6
            self.min_dist = 0 
            self.max_dist = 500
            self.tbefore_sec = 100
            self.tafter_sec = 300
            self.network = 'AK,AT,AV,CN,II,IU,US,XM,XV,XZ,YV'  # note: cannot use '*' because of IM
            self.channel = 'BH?'
            self.network = 'AV'
            self.station = 'SPBG,KABU'                      # For testing data gaps 
            self.use_catalog = 0                            # To get (lat,lon, etime, dep, mag) from some catalog = 1 OR use defined = 0 (see iex=9)

# Error from util_write_cap.py
#   util_write_cap.py, line 437, in add_sac_metadata
#    and pick.waveform_id.channel_code[2].upper() == 'Z' and pick.waveform_id.location_code == tr.stats.location and pick.phase_hint == 'Pn'):
#        IndexError: string index out of range
#
# This error happens when allowing data already rotated (VRT). see getwaveform_llnl.py Line 94
        if iex == 10:
            self.idb = 3              # LLNL database
            self.otime = obspy.UTCDateTime("1982-08-05T14:00:00.000000Z")
            self.min_dist = 0 
            self.max_dist = 1200
            self.tbefore_sec = 100
            self.tafter_sec = 600
            self.network = '*'        # note that the client will look for BK stations in the list
            self.channel = 'BH?'      # ALL channels from LLNL are returned regardless
            self.scale_factor = 10.0**2  # original
            self.overwrite_ddir = 0
            self.resample_freq = 20.0 

# error applying rotation to station LL.KNB. See util_write_cap.py Line 161
# (Disable the try block and only allow Line 161 to create this problem)
# 
#   File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/signal/rotate.py", line 257, in rotate2zne
#    x, y, z = np.dot(_t, [data_1, data_2, data_3])
#    ValueError: setting an array element with a sequence.
        if iex == 11:
            self.idb = 1 
            self.otime = obspy.UTCDateTime("1995-07-31T12:34:46.860000Z")
            self.min_dist = 0 
            self.max_dist = 1200
            self.tbefore_sec = 100
            self.tafter_sec = 600
            self.network = '*'        # note that the client will look for BK stations in the list
            self.channel = 'BH?'      # ALL channels from LLNL are returned regardless
            self.scale_factor = 10.0**2  # original
            self.overwrite_ddir = 0
            self.resample_freq = 20.0 

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
            self.idb = 1 
            self.otime = obspy.UTCDateTime("1999-01-23T03:00:33.200000Z") # frenchman flat 1
            # self.otime = obspy.UTCDateTime("1999-01-27T10:44:23.310000Z") # frenchman flat 2
            # self.otime = obspy.UTCDateTime("2000-01-30T14:46:51.310000Z") # trona mine 2
            # self.otime = obspy.UTCDateTime("2002-06-14T12:40:44.450000Z") # little skull
            self.min_dist = 0 
            self.max_dist = 1200
            self.tbefore_sec = 100
            self.tafter_sec = 600
            self.network = '*'        # note that the client will look for BK stations in the list
            self.channel = 'BH?'      # ALL channels from LLNL are returned regardless
            self.scale_factor = 10.0**2  # original
            self.overwrite_ddir = 0
            self.resample_freq = 20.0 

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
            self.idb = 1 
            self.otime = obspy.UTCDateTime("1997-06-14T19:48:19.930000Z") # indian springs
            # self.otime = obspy.UTCDateTime("1996-09-05T08:16:56.090000Z") # amargosa
            self.min_dist = 0 
            self.max_dist = 1200
            self.tbefore_sec = 100
            self.tafter_sec = 600
            self.network = '*'        # note that the client will look for BK stations in the list
            self.channel = 'BH?'      # ALL channels from LLNL are returned regardless
            self.scale_factor = 10.0**2  # original
            self.overwrite_ddir = 0
            self.resample_freq = 20.0 

# time > 2007-01-24T11:30:06.100000Z time < 2007-01-24T11:30:26.100000Z
#Traceback (most recent call last):
#  File "run_getwaveform.py", line 478, in <module>
#    ev = cat.filter(mintime_str, maxtime_str)[0]
#  File "/home/alvizuri/miniconda2/envs/sln/lib/python3.5/site-packages/obspy/core/event/catalog.py", line 162, in __getitem__
#    return self.events.__getitem__(index)
#IndexError: list index out of range
        if iex == 14:
            self.idb = 3 
            self.otime = obspy.UTCDateTime("2007-01-24T11:30:16.100000Z") # ralston
            self.min_dist = 0 
            self.max_dist = 1200
            self.tbefore_sec = 100
            self.tafter_sec = 600
            self.network = '*'        # note that the client will look for BK stations in the list
            self.channel = 'BH?'      # ALL channels from LLNL are returned regardless
            self.scale_factor = 10.0**2  # original
            self.overwrite_ddir = 0
            self.resample_freq = 20.0 

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
            self.idb = 1 
            self.otime = obspy.UTCDateTime("1995-07-31T12:34:46.860000Z") # timber mountain
            self.min_dist = 0 
            self.max_dist = 1200
            self.tbefore_sec = 100
            self.tafter_sec = 600
            self.network = '*'        # note that the client will look for BK stations in the list
            self.channel = 'BH?'      # ALL channels from LLNL are returned regardless
            self.scale_factor = 10.0**2  # original
            self.overwrite_ddir = 0
            self.resample_freq = 20.0 

# Iniskin earthquake - all channels at IU.COLA and II.KDAK
# note different channels and location codes for strong motion:
#       IU.COLA.20.HNZ.sac       http://ds.iris.edu/mda/IU/COLA
#       II.KDAK.00.ENZ.sac       http://ds.iris.edu/mda/II/KDAK
        if iex == 16:
            self.idb = 1
            self.overwrite_ddir = 1       # delete data dir if it exists
            self.use_catalog = 0          # do not use event catalog for source parameters
            # GCMT source parameters
            # the otime is the centroid time and accounts for tshift
            self.otime = obspy.UTCDateTime("2016-01-24T10:30:37.400") 
            self.elat = 59.75
            self.elon = -153.27
            self.edep = 110700
            self.emag = 7.1
            # subset of stations
            self.min_dist = 0
            self.max_dist = 800
            self.tbefore_sec = 100
            self.tafter_sec = 600
            self.network = 'IU,II'
            # self.channel = '?H?,?N?'
            self.channel = 'HH?,BH?,BN?,HN?,EN?'
            self.resample_freq = 0        # no resampling
            self.scale_factor = 1         # no scale factor

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
            self.idb = 3            # LLNL
            self.otime = obspy.UTCDateTime("1989-06-27T15:30:00.02")
            self.elat = 37.275
            self.elon = -116.354
            self.edep = 640
            self.emag = 4.90
            self.min_dist = 0 
            self.max_dist = 1200
            self.tbefore_sec = 100
            self.tafter_sec = 600
            self.network = '*' 
            self.channel = 'BH?,LH?' 
            self.overwrite_ddir = 0

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
            self.idb = 1
            self.overwrite_ddir = 1       # delete data dir if it exists
            self.use_catalog = 0          # do not use event catalog for source parameters
            # GCMT source parameters
            # the otime is the centroid time and accounts for tshift
            self.otime = obspy.UTCDateTime("2015-12-02T10:05:25.798") 
            self.elat = 61.70
            self.elon = -147.26
            self.edep = 36590
            self.emag = 4.50
            # subset of stations
            self.min_dist = 300
            self.max_dist = 400
            self.tbefore_sec = 100
            self.tafter_sec = 600
            # self.network = 'AK,AT,AV,CN,II,IU,US,XM,XV,XZ,YV,ZE'
            self.network = 'IU'
            self.channel = 'HH?,BH?'
            # self.resample_freq = 0        # no resampling -- THIS WORKS
            self.resample_freq = 50       # THIS FAILS BUT ONLY AFTER REPEATING
            self.scale_factor = 1         # no scale factor

# Test case for UVW
        if iex == 19:
            self.idb = 1
            self.use_catalog = 0  
            self.otime = obspy.UTCDateTime("2016-12-08T10:16:00")
            self.elat = 64.2380
            self.elon = -150.0581
            self.edep = 18507
            self.emag = 4.60
            self.station = 'F3TN'
            
            self.min_dist = 0
            self.max_dist = 300
            self.tbefore_sec = 0
            self.tafter_sec = 600
            self.network = 'XV'
            self.channel = 'HH?'
            self.scale_factor = 1
            self.resample_freq = 0
            self.detrend = False
            self.demean = False
            self.taper = False
            self.ipre_filt = 1

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
            self.idb = 3            # LLNL
            self.otime = obspy.UTCDateTime("1990-06-13T16:00:00.09")
            self.elat = 37.262
            self.elon = -116.421
            self.edep = 674
            self.emag = 5.34
            self.min_dist = 0 
            self.max_dist = 1200
            self.tbefore_sec = 100
            self.tafter_sec = 600
            self.network = '*' 
            self.channel = 'BH?,LH?' 
            self.overwrite_ddir = 0

#=================================================================================
# CATEGORY B EXAMPLES: important events
# 1XX: southern Alaska
# 2XX: central Alaska
# 9XX: other
#=================================================================================
        # SilwalTape2016 example event (Anchorage)
        if iex == 100:
            self.otime = obspy.UTCDateTime("2009-04-07T20:12:55.351")
            self.min_dist = 0 
            self.max_dist = 500
            self.tbefore_sec = 100
            self.tafter_sec = 300
            self.network = 'AK,AT,AV,CN,II,IU,US,XM,XV,XZ,YV'  # note: cannot use '*' because of IM
            self.channel = 'BH?'
            self.use_catalog = 0 
            self.elat = 61.45420
            self.elon = -149.7428
            self.edep = 33033.60
            self.emag = 4.6


# Iniskin earthquake
# NOTE: must enter username and password above to get SALMON (ZE) stations
        if iex == 101:
            self.idb = 1
            self.overwrite_ddir = 1       # delete data dir if it exists
            self.use_catalog = 0          # do not use event catalog for source parameters
            # GCMT source parameters
            # the otime is the centroid time and accounts for tshift
            self.otime = obspy.UTCDateTime("2016-01-24T10:30:37.400") 
            self.elat = 59.75
            self.elon = -153.27
            self.edep = 110700
            self.emag = 7.1
            # subset of stations
            self.min_dist = 0
            self.max_dist = 800
            self.tbefore_sec = 100
            self.tafter_sec = 600
            self.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US'  # IM will probably crash it
            self.channel = 'BH?,HH?'
            self.resample_freq = 0        # no resampling
            self.scale_factor = 1         # no scale factor

            # parameters for examining step response (causal low-pass filter on raw waveforms)
            # delete AUQ, SPCP, SPBG
            self.rotateRTZ = False
            self.ifFilter = True
            self.filter_type = 'lowpass'
            self.f1 = 1/4
            self.zerophase = False
            self.remove_response = False
            self.ipre_filt = 0
            self.demean = False
            self.detrend = False

# Iniskin earthquake - all strong motion
        if iex == 103:
            self.idb = 1
            self.overwrite_ddir = 1       # delete data dir if it exists
            self.use_catalog = 0          # do not use event catalog for source parameters
            # GCMT source parameters
            # the otime is the centroid time and accounts for tshift
            self.otime = obspy.UTCDateTime("2016-01-24T10:30:37.400") 
            self.elat = 59.75
            self.elon = -153.27
            self.edep = 110700
            self.emag = 7.1
            # subset of stations
            self.min_dist = 0
            self.max_dist = 800
            self.tbefore_sec = 100
            self.tafter_sec = 600
            self.network = '*'
            self.channel = 'BN?,HN?,EN?'
            self.resample_freq = 0        # no resampling
            self.scale_factor = 1         # no scale factor

# Totschunda fault
        if iex == 104:
            self.idb = 1
            self.overwrite_ddir = 1       # delete data dir if it exists
            self.use_catalog = 0          # do not use event catalog for source parameters
            # GCMT source parameters
            # the otime is the centroid time and accounts for tshift
            self.otime = obspy.UTCDateTime("2017-02-13T07:17:12.000") 
            self.elat = 62.5154
            self.elon = -142.7485
            self.edep = 7564
            self.emag = 5.3
            # subset of stations
            self.min_dist = 0
            self.max_dist = 500
            self.tbefore_sec = 50
            self.tafter_sec = 300
            self.network = 'AV'
            self.channel = 'BH?,HH?'
            self.resample_freq = 50        # no resampling
            self.scale_factor = 100         # no scale factor

# Klukwan earthquakes
        if iex == 105:
            self.idb = 1
            self.overwrite_ddir = 1       # delete data dir if it exists
            self.use_catalog = 0          # do not use event catalog for source parameters
            # GCMT source parameters
            # the otime is the centroid time and accounts for tshift
            self.otime = obspy.UTCDateTime("2017-05-01T12:31:53.000") 
            self.elat = 59.8522
            self.elon = -136.6618
            self.edep = 4000
            self.emag = 6.2
            self.otime = obspy.UTCDateTime("2017-05-01T14:18:14.000") 
            self.elat = 59.8184
            self.elon = -136.7163
            self.edep = 1000
            self.emag = 6.0
    
            # subset of stations
            self.min_dist = 0
            self.max_dist = 500
            self.tbefore_sec = 100
            self.tafter_sec = 500
            self.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
            self.channel = 'BH?,HH?'
            self.resample_freq = 50        # no resampling
            self.scale_factor = 100         # no scale factor

# Kantishna earthquakes
        if iex == 106:
            self.idb = 1
            self.overwrite_ddir = 1       # delete data dir if it exists
            self.use_catalog = 0          # do not use event catalog for source parameters
            # GCMT source parameters
            # the otime is the centroid time and accounts for tshift
            self.otime = obspy.UTCDateTime("2017-04-29T11:15:48.000") 
            self.elat = 63.1296
            self.elon = -151.1517
            self.edep = 11000
            self.emag = 5.2
            self.otime = obspy.UTCDateTime("2017-01-31T09:38:37.000") 
            self.elat = 63.0817
            self.elon = -150.9427
            self.edep = 135000
            self.emag = 5.2
    
            # subset of stations
            self.min_dist = 0
            self.max_dist = 500
            self.tbefore_sec = 100
            self.tafter_sec = 500
            self.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
            self.channel = 'BH?,HH?'
            self.resample_freq = 50        # no resampling
            self.scale_factor = 100         # no scale factor

# Cook Inlet earthquake
        if iex == 107:
            self.idb = 1
            self.overwrite_ddir = 1       # delete data dir if it exists
            self.use_catalog = 0          # do not use event catalog for source parameters
            # GCMT source parameters
            # the otime is the centroid time and accounts for tshift
            self.otime = obspy.UTCDateTime("2017-05-07T04:25:19.000") 
            self.elat = 60.1945
            self.elon = -151.6743
            self.edep = 64000
            self.emag = 5.2
    
            # subset of stations
            self.min_dist = 0
            self.max_dist = 500
            self.tbefore_sec = 100
            self.tafter_sec = 500
            self.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
            self.channel = 'BH?,HH?'
            self.resample_freq = 50        # no resampling
            self.scale_factor = 100         # no scale factor
            
# MFFZ earthquakes for investigating the step response
# LISTED IN CHRONOLOGICAL ORDER
        if iex == 200:
            self.idb = 1
            self.use_catalog = 0
            # -------------------------------------------------
            self.otime = obspy.UTCDateTime("2014-08-31T03:06:57.111")
            self.elat = 65.1526
            self.elon = -149.0398
            self.edep = 16614.7
            self.emag = 5.20
            # -------------------------------------------------
            self.otime = obspy.UTCDateTime("2014-10-21T00:36:58.333")
            self.elat = 65.1489
            self.elon = -149.0413
            self.edep = 13134.8
            self.emag = 4.90
            # -------------------------------------------------
            self.otime = obspy.UTCDateTime("2014-10-23T16:30:23.968")
            self.elat = 65.1644
            self.elon = -149.0523
            self.edep = 20066.5
            self.emag = 5.00
            # -------------------------------------------------
            self.otime = obspy.UTCDateTime("2015-09-12T03:25:12.711")
            self.elat = 65.1207
            self.elon = -148.6646
            self.edep = 15556.8
            self.emag = 2.63
            # -------------------------------------------------
            self.otime = obspy.UTCDateTime("2015-10-31T02:56:35.572")
            self.elat = 64.4285
            self.elon = -149.6969
            self.edep = 23852.1
            self.emag = 3.47
            # -------------------------------------------------
            self.otime = obspy.UTCDateTime("2016-01-14T19:04:10.727")
            self.elat = 64.6827
            self.elon = -149.2479
            self.edep = 22663.7
            self.emag = 3.80
            # -------------------------------------------------
            self.otime = obspy.UTCDateTime("2016-11-06T9:29:10.579")
            self.elat = 64.1639
            self.elon = -150.0626
            self.edep = 23190.0
            self.emag = 4.00
            # -------------------------------------------------
            self.otime = obspy.UTCDateTime("2016-12-08T10:18:13.868")
            self.elat = 64.1937
            self.elon = -150.0376
            self.edep = 24522.1
            self.emag = 4.30
            # -------------------------------------------------
            self.otime = obspy.UTCDateTime("2017-05-08T05:09:02.000") 
            self.elat = 65.2643
            self.elon = -146.922
            self.edep = 9000   # AEC/Vipul
            self.emag = 3.60   # Vipul
            # -------------------------------------------------
            self.otime = obspy.UTCDateTime("2017-06-28T12:58:52.000") 
            self.elat = 64.7569
            self.elon = -148.8883
            self.edep = 18000  # Vipul
            self.emag = 3.50   # Vipul
            # -------------------------------------------------
            # VIPUL: WHAT ARE THESE? OTHER EVENTS?
            # self.otime = obspy.UTCDateTime("2015-03-30T12:33:19.000")
            # self.otime = obspy.UTCDateTime("2015-10-20T19:14:16.000")
            # self.otime = obspy.UTCDateTime("2011-12-21T16:28:41.000")
    # -------------------------------------------------
            
            self.min_dist = 0
            self.max_dist = 300
            self.tbefore_sec = 200
            self.tafter_sec = 600
            self.network = 'AV,CN,AT,TA,AK,XV,II,IU,US'
            # self.network = 'XV,AK'
            self.channel = 'BH?,HH?'
            self.scale_factor = 1
            self.resample_freq = 0
            # for CAP
            # self.scale_factor = 10.0**2
            # self.resample_freq = 50 
            
            # to investigate step response
            self.rotateRTZ = True
            self.rotateUVW = True
            self.ifFilter = True
            self.zerophase = False    # causal
            # filter_type = 'lowpass'
            self.filter_type = 'bandpass'
            self.f1 = 1/100  # fmin
            self.f2 = 1/10  # fmax
            self.corners = 4
            # self.remove_response = False
            self.remove_response = True
            self.ipre_filt = 1
            self.demean = True
            self.detrend = True
            self.output_cap_weight_file = False
            # self.outformat = 'DISP'
            self.ifsave_sacpaz = True
            self.taper = 0.2

# NENNUC event (from Steve)
        if iex == 201:
            self.idb = 1
            self.use_catalog = 0          # manually enter AEC catalog parameters
            self.otime = obspy.UTCDateTime("2016-01-14T19:04:10.727")
            self.elat = 64.6827
            self.elon = -149.2479
            self.edep = 22663.7
            self.emag = 3.8
            # -------------------------------------------------
            self.otime = obspy.UTCDateTime("2015-09-12T03:25:12.711")
            self.elat = 65.1207
            self.elon = -148.6646
            self.edep = 15556.8
            self.emag = 2.6
            # -------------------------------------------------
            self.otime = obspy.UTCDateTime("2013-03-12T07:39:50.214")
            self.elat = 64.7161
            self.elon = -148.9505
            self.edep = 20000
            self.emag = 2.1
            
            # For CAP
            self.min_dist = 0
            self.max_dist = 200
            self.tbefore_sec = 500
            self.tafter_sec = 2000
            self.taper = 0.1
            self.network = 'AV,CN,AT,TA,AK,XV,II,IU,US'
            # self.network = 'AK,TA,II,IU,US'
            # self.station = "-RAG"  # RAG has choppy data; gives error: Exception: Can't merge traces with same ids but differing sampling rates!
            self.channel = 'BH?,HH?'
            self.pre_filt = (0.0025, 0.003, 10.0, 15.0) # BH default
            self.ipre_filt = 1

            # to investigate step response
            # self.self.tbefore_sec = 2000
            # self.tafter_sec = 2000
            # self.resample_freq = 0 
            # self.scale_factor = 1
            # self.ifFilter = True
            # self.zerophase = False
            # self.filter_type = 'bandpass'
            # self.f1 = 1/100
            # self.f2 = 1/20
            # self.remove_response = True
            # self.ipre_filt = 0
            # self.demean = True
            # self.detrend = True

# Kyle Nenana basin earthquakes for basin amplification
        if iex == 210:
            self.idb = 1
            self.overwrite_ddir = 1       # delete data dir if it exists
            self.use_catalog = 0          # do not use event catalog for source parameters
            # AEC source parameters
            self.otime = obspy.UTCDateTime("2015-11-20T10:53:48.168") 
            self.elat = 64.6210
            self.elon = -149.4024
            self.edep = 17113.4
            self.emag = 2.67
            # subset of stations
            self.min_dist = 0
            self.max_dist = 200
            self.tbefore_sec = 100
            self.tafter_sec = 400
            self.network = 'AK,AT,II,IU,US,XM,XV,XZ,TA'  # no CN,AV,YV,ZE
            self.channel = 'BH?,HH?'
            # self.resample_freq = 0        # no resampling
            # self.scale_factor = 1         # no scale factor
            # For CAP moment tensor
            self.resample_freq = 50         # same as greens function 
            self.scale_factor = 100         # change from m/s to cm/s
            # self.ipre_filt = 0
            self.remove_response = True
            # self.demean = False
            # self.detrend = False

# gets waveforms from M 8.3 Chile event with stations centered in Minto
        if iex == 211:
            self.idb = 1
            self.overwrite_ddir = 1       # delete data dir if it exists
            self.use_catalog = 0
            self.otime = obspy.UTCDateTime("2015-09-16T22:54:33.000")
            self.elat = -31.5695
            self.elon = -71.6543
            self.edep = 22400.0
            self.emag = 8.30
            self.rlat = 64.7716
            self.rlon = -149.1465
            self.rtime = obspy.UTCDateTime("2015-09-16T23:09:15.000")
            self.tbefore_sec = 100
            self.tafter_sec = 200
            self.min_dist = 0
            self.max_dist = 100
            self.network = 'AK,AT,II,IU,US,XM,XV,XZ,TA'  # no CN,AV,YV,ZE
            self.channel = 'BH?,HH?'
            self.resample_freq = 0        
            self.scale_factor = 1        
            self.remove_response = True

# gets a 
        if iex == 212:
            self.idb = 1
            self.overwrite_ddir = 1       # delete data dir if it exists
            self.use_catalog = 0
            self.otime = obspy.UTCDateTime("2016-06-06T00:00:00.000")
            self.elat = 64.6130
            self.elon = -149.0992
            self.edep = 0
            self.emag = 0.00
            # self.rlat = 64.7716
            # self.rlon = -149.1465
            # self.rtime = obspy.UTCDateTime("2015-09-16T23:09:15.000")
            self.tbefore_sec = 0
            self.tafter_sec = 3600
            self.min_dist = 0
            self.max_dist = 100
            self.network = 'XV,AK,TA'  # no CN,AV,YV,ZE
            self.channel = 'HH?,BH?'
            self.resample_freq = 0        
            self.scale_factor = 1        
            self.remove_response = True
            self.rotateRTZ = False
            # self.pre_filt=(f0*0.001, f1*0.001, f2*1000, f3*1000)
            # self.ipre_filt = 2
            self.ipre_filt = 0

# Chatnika earthquake
        if iex == 213:
            self.idb = 1
            self.overwrite_ddir = 1       # delete data dir if it exists
            self.use_catalog = 0          # do not use event catalog for source parameters
            # GCMT source parameters
            # the otime is the centroid time and accounts for tshift
            self.otime = obspy.UTCDateTime("2017-05-08T05:09:02.000") 
            self.elat = 65.2643
            self.elon = -146.922
            self.edep = 9000
            self.emag = 3.8
    
            # subset of stations
            self.min_dist = 0
            self.max_dist = 500
            self.tbefore_sec = 100
            self.tafter_sec = 500
            self.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
            self.channel = 'BH?,HH?'
            self.resample_freq = 50        # no resampling
            self.scale_factor = 100         # no scale factor

# NE Nenana earthquake
        if iex == 214:
            self.idb = 1
            self.overwrite_ddir = 1       # delete data dir if it exists
            self.use_catalog = 0          # do not use event catalog for source parameters
            # GCMT source parameters
            # the otime is the centroid time and accounts for tshift
            self.otime = obspy.UTCDateTime("2017-06-28T12:58:52") 
            self.elat = 64.7569
            self.elon = -148.8883
            self.edep = 18000
            self.emag = 3.5
            
            # subset of stations
            self.min_dist = 0
            self.max_dist = 500
            self.tbefore_sec = 100
            self.tafter_sec = 500
            self.network = 'AV,CN,ZE,AT,TA,AK,XV,II,IU,US' 
            self.channel = 'BH?,HH?'
            self.resample_freq = 50        
            self.scale_factor = 100        
#----------------------------------------------------------

#-----------------------------------------------------------
# Event from HutchisonGhosh2016
# Crashes when using all networks [i.e. network = '*']
        if iex == 900:
            self.idb = 1
            self.overwrite_ddir = 1       # delete data dir if it exists
            self.use_catalog = 0          # do not use event catalog for source parameters
            # GCMT source parameters
            # the otime is the centroid time and accounts for tshift
            self.otime = obspy.UTCDateTime("2014-12-12T00:13:17") 
            self.elat = 49.10
            self.elon = -124.05
            self.edep = 60000
            self.emag = 4.10
            # subset of stations
            self.min_dist = 0
            self.max_dist = 250
            self.tbefore_sec = 100
            self.tafter_sec = 300
            self.network = '*'
            self.channel = 'LH?,BH?'
            self.resample_freq = 50        # no resampling
            self.scale_factor = 100         # no scale factor

            # For plotting filtered waveforms
            self.tbefore_sec = 500
            self.tafter_sec = 500
            self.resample_freq = 0 
            self.scale_factor = 1
            self.ifFilter = True
            self.zerophase = True
            self.filter_type = 'bandpass'
            self.f1 = 1/50
            self.f2 = 1/20
            self.remove_response = True
            self.ipre_filt = 2
            self.demean = True
            self.detrend = True
            self.taper = True
            
#=================================================================================
# END EXAMPLES
#=================================================================================
