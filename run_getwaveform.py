# To run this script:
# > python run_getwaveform.py
import obspy
import copy

# EXAMPLES (choose one)
iex = 1

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
resample_freq = 20.0                      # =0 for no resampling
scale_factor = 0                          # 10**2 to convert m/s to cm/s

# username and password for accessing embargoed data from IRIS
# Register here: http://ds.iris.edu/ds/nodes/dmc/forms/restricted-data-registration/
# Run example iex = 4 to check
user = ""
password = ""

# (Do not use)
# This is a template for testing before creating an example.
# We can delete this if it creates more issues.
if iex == 0:
    otime = obspy.UTCDateTime("2009-04-07T20:12:55")
    min_dist = 0 
    max_dist = 500
    tbefore_sec = 100
    tafter_sec = 300
    network = 'XP'
    channel = '*Z'

# SilwalTape2016 example event (Anchorage)
if iex == 1:
    otime = obspy.UTCDateTime("2009-04-07T20:12:55")
    min_dist = 0 
    max_dist = 500
    tbefore_sec = 100
    tafter_sec = 300
    network = 'AK,AT,AV,CN,II,IU,XM,XV,XZ,YV'  # note: cannot use '*' because of IM
    channel = 'BH*'

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
    channel = 'BH*'

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
    network = 'AK'       # works fine 
    channel = 'BH*'

# SALMON example (restricted data from IRIS)
if iex == 4:
    otime = obspy.UTCDateTime("2016-01-24T10:30:29.557")
    min_dist = 0 
    max_dist = 700
    tbefore_sec = 100
    tafter_sec = 600
    network = 'AK,ZE'   # ZE waveforms not returned (only ZE.MPEN)
    channel = 'BH*,HH*'

# ERROR EXAMPLE
# PROBLEM 1: output file names should be BHR and BHT (not BH1 and BH2)
# PROBLEM 2: output files are NOT being rotated
# PROBLEM 3: output file names should be EID.NN.SSS.LL.CCC.sac
# MAIN PROBLEM: All rotations should be based on the azimuth of the sensor
#    (CMPAZ header in sac), which, combined with the station-source backazimuth,
#    will provide the rotation to radial and transverse components.
#    The rotation should NOT depend on channel names like E or N.
#    (Even for GSN stations the E and N do not point exactly E and N.)
if iex == 5:
    otime = obspy.UTCDateTime("2016-01-24T10:30:29.557")
    min_dist = 0 
    max_dist = 500
    tbefore_sec = 100
    tafter_sec = 600
    network = 'II'
    channel = 'BH*'

# nuclear event: LLNL (see also iex = 7)
# GOAL: To find events in the LLNL database based on a target origin time,
#       rather than an eid. Perhaps a second parameter is needed to say
#       "find closest event within 30 seconds of target time".
#       (The reference time (NZYEAR, etc) should then be assined as the actual origin time,
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
    #otime = obspy.UTCDateTime("1991-09-14T19:00:08.031Z")
    evid = 635527        # Hoya event id in LLNL database

    min_dist = 0 
    max_dist = 1200
    tbefore_sec = 100
    tafter_sec = 600
    network = '*'        # note that the client will look for BK stations in the list
    channel = 'BH*'      # ALL channels from LLNL are returned
    #scale_factor = 10.0**2  # original
    scale_factor = 2e-1     # Hoya  

# same as iex=6 but for the IRIS database
# GOAL: For LLNL events, we do NOT want to use the IRIS source parameters:
#       origin time, hypocenter, magnitude.
#       (Unsure what the threshold is for IRIS catalog.)
#       We need a flag to decide whether to find/use source parameters in the IRIS
#       catalog, or to assign them ourselves.
if iex == 7:
    idb = 1            # IRIS database
    #resample_freq = 0  # no resampling
    #otime = obspy.UTCDateTime("1991-09-14T19:00:00.000Z")   # Hoya actual
    otime = obspy.UTCDateTime("1991-09-14T19:00:08.031Z")   # Hoya target
    min_dist = 0 
    max_dist = 1200
    tbefore_sec = 100
    tafter_sec = 600
    network = '*'
    channel = 'BH*,LH*' 

# fetch and process waveforms
if idb == 1:
    # import functions to access waveforms
    import getwaveform_iris
    from obspy.clients.fdsn import Client
    if not user and not password:
        client = Client("IRIS")
    else:
        client = Client("IRIS",user=user,password=password)
    # will only work for events in the 'IRIS' catalog
    # (future: for Alaska events, read the AEC catalog)
    cat = client.get_events(starttime = otime-10, endtime = otime+10)

    # Extract waveforms, IRIS
    getwaveform_iris.run_get_waveform(c = client, event = cat[0], min_dist = min_dist, max_dist = max_dist, 
            before = tbefore_sec, after = tafter_sec, network = network, channel = channel, 
            resample_freq = resample_freq, ifrotate = rotate,
            ifCapInp = output_cap_weight_file, ifRemoveResponse = remove_response,
            ifDetrend = detrend, ifDemean = demean, ifEvInfo = output_event_info, 
            scale_factor = scale_factor, pre_filt = pre_filt)

if idb == 3:
    import llnl_db_client
    import getwaveform_llnl
    client = llnl_db_client.LLNLDBClient(
            "/store/raw/LLNL/UCRL-MI-222502/westernus.wfdisc")

    # Extract waveforms, LLNL
    getwaveform_llnl.run_get_waveform(llnl_db_client = client, event = evid, 
            network = network, channel = channel, 
            min_dist = min_dist, max_dist = max_dist, 
            before = tbefore_sec, after = tafter_sec, 
            resample_freq = resample_freq,
            ifrotate=True, ifCapInp=True, ifRemoveResponse=True,
            ifDetrend=True, ifDemean=True, ifEvInfo=True,
            scale_factor = scale_factor,
            pre_filt=(0.005, 0.006, 10.0, 15.0))
