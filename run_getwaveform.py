# To run this script:
# > python run_getwaveform.py
import obspy
import copy

# EXAMPLES (choose one)
iex = 1

# default settings
idb = 1    # default: =1-IRIS; =2-AEC; =3-LLNL
# username and password for accessing embargoed data from IRIS
# Register here: http://ds.iris.edu/ds/nodes/dmc/forms/restricted-data-registration/
# Run example iex = 4 to check
user = ""
password = ""
# Pre-processing (manily for CAP)
rotate = True
output_cap_weight_file = True
remove_response = True
detrend = True
demean = True
output_event_info = True
pre_filt = (0.005, 0.006, 5.0, 10.0)      # Why is this needed?
resample_freq = 20.0                      # =0 for no resampling
scale_factor = 10.0**2                    # =10.0**2 (convert m/s to cm/s)

# SilwalTape2016 example event (Anchorage)
if iex == 1:
    otime = obspy.UTCDateTime("2009-04-07T20:12:55")
    min_dist = 0 
    max_dist = 500
    before = 100
    after = 300
    network = 'AK,AT,AV,CN,II,IU,XM,XV,XZ,YV'  # note: cannot use '*' because of IM
    channel = 'BH*'

# ERROR EXAMPLE [obspy]
# PROBLEM: No waveforms are returned -- perhaps related to the tbefore request
# ERROR MESSAGE: ValueError: The length of the input vector x must be at least padlen, which is 39.
# KLUDGE: change the before time to 41 or do not use AV network
if iex == 2:
    otime = obspy.UTCDateTime("2016-01-24T10:30:29.557")
    min_dist = 0 
    max_dist = 500
    before = 42          # Crashes
    before = 41         # works fine
    after = 600
    network = 'AV'       # Crashes when -  before = 42; Works fine when - before = 41
    channel = 'BH*'

# ERROR EXAMPLE [obspy]
# PROBLEM: If a particular network is requested (either explicitly or within *), no waveforms are returned.
# ERROR MESSAGE: NotImplementedError: ResponseListResponseStage not yet implemented due to missing example data. Please contact the developers with a test data set (waveforms and StationXML metadata).
# KLUDGE: list all networks explicitly, except IM
if iex == 3:
    otime = obspy.UTCDateTime("2009-04-07T20:12:55")
    min_dist = 300 
    max_dist = 500
    before = 100
    after = 300
    #network = '*'       # crashes
    #network = 'AK,IM'   # crashes (Error because of the IM network)
    network = 'AK'       # works fine 
    channel = 'BH*'

# ERROR EXAMPLE [IRIS + obspy]
# PROBLEM: cannot get embargoed waveforms (ZE waveforms = SALMON)
# --> need to consult with IRIS and/or obspy
if iex == 4:
    otime = obspy.UTCDateTime("2016-01-24T10:30:29.557")
    min_dist = 0 
    max_dist = 700
    before = 100
    after = 600
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
    before = 100
    after = 600
    network = 'II'
    channel = 'BH*'

# ERROR EXAMPLE LLNL #1
# see also iex = 7
# PROBLEM: SAC headers 'b' and 'e' are set based on user-secified 'before' and 'after',
#          whereas they need to be set based on whatever data is available.
# KLUDGE: use large (incorrect) time shifts in CAP
if iex == 6:
    # to get the LLNL client, which is a private repo from Lion Krischer:
    # cd $REPOS
    # git clone https://GITHUBUSERNAME@github.com/krischer/llnl_db_client.git
    # then follow instructions for install
    idb = 3 # NOTE this is a request for LLNL. See also iex=7 (IRIS)
    otime = obspy.UTCDateTime("1991-09-14T19:00:00.000Z")   # Hoya
    evid = 635527   # Hoya
    min_dist = 0 
    max_dist = 1200
    before = 20
    after = 600
    network = '*'  # note that the client will look for BK stations in the list
    channel = 'BH*'
    #scale_factor = 10.0**2  # original
    scale_factor = 2e-1     # Hoya  

# ERROR EXAMPLE FOR HOYA EVENT
# see also iex = 6
# PROBLEM: there is an 80 second offset between IRIS and LLNL data. For example
# see station ANMO.
if iex == 7:
    idb = 1 # NOTE this is a request for IRIS. See also iex=6 (LLNL)
    otime = obspy.UTCDateTime("1991-09-14T19:00:00.000Z")   # Hoya
    evid = 635527   # Hoya
    min_dist = 0 
    max_dist = 1200
    before = 20
    after = 600
    network = '*'  # note that the client will look for BK stations in the list
    channel = 'BH*'
    #scale_factor = 10.0**2  # original
    scale_factor = 2e-1     # Hoya  

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
            before = before, after = after, network = network, channel = channel, 
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
            before = before, after = after, 
            resample_freq = resample_freq,
            ifrotate=True, ifCapInp=True, ifRemoveResponse=True,
            ifDetrend=True, ifDemean=True, ifEvInfo=True,
            scale_factor = scale_factor,
            pre_filt=(0.005, 0.006, 10.0, 15.0))
