# To run this script:
# > python run_getwaveform.py
import obspy
import copy
import util_helpers
import shutil   # only used for deleting data directory
import os

# EXAMPLES (choose one)
iex = 10

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
station = '*'                # all stations
overwrite_ddir = 0          # 1 = delete data directory if it already exists

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
    channel = 'BH*'

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
    channel = 'BH*'      # ALL channels from LLNL are returned
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
    channel = 'BH*,LH*' 
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
    channel = 'BH*,LH*,EH*' 

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
    network = 'AK,AT,AV,CN,II,IU,XM,XV,XZ,YV'  # note: cannot use '*' because of IM
    channel = 'BH*'
    network = 'AV'
    station = 'SPBG,KABU'                      # For testing data gaps 
    use_catalog = 0                            # To get (lat,lon, etime, dep, mag) from some catalog = 1 OR use defined = 0 (see iex=9)

# Error from util_write_cap.py
#   unsupported operand type(s) for +: 'float' and 'UTCDateTime'
# QUESTION Why the error? both are UTCDateTime
#  evtime = UTCDateTime("1992-06-29T10:14:22.280000Z"), 
#  tr.stats.edntime = UTCDateTime("1992-06-29T10:17:49.103995Z")
if iex == 10:
    idb = 3              # LLNL database
    otime = obspy.UTCDateTime("1992-06-29T10:14:22.480000")
    min_dist = 0 
    max_dist = 1200
    tbefore_sec = 100
    tafter_sec = 600
    network = '*'        # note that the client will look for BK stations in the list
    channel = 'BH*'      # ALL channels from LLNL are returned regardless
    scale_factor = 10.0**2  # original
    overwrite_ddir = 0
    resample_freq = 20.0 
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
            scale_factor = scale_factor, pre_filt = pre_filt)

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
