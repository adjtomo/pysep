# To run this script:
# > python run_getwaveform.py
import obspy
import copy

# EXAMPLES:
iex = 1    # Choose example to run (default = 1)
idb = 1    # default: =1-IRIS; =2-AEC; =3-LLNL

# default settings
rotate = True
output_cap_weight_file = True
remove_response = True
detrend = True
demean = True
output_event_info = True
pre_filt = (0.005, 0.006, 5.0, 10.0)    # Why is this needed?
resample_freq = 20.0                      # =0 for no resampling
scale_factor = 10.0**2                    # =10.0**2 (convert m/s to cm/s)

# WORKING EXAMPLE (default)
# SilwalTape2016 example event (Anchorage)
if iex == 1:
    otime = obspy.UTCDateTime("2009-04-07T20:12:55")
    min_dist = 0 
    max_dist = 500
    before = 100
    after = 300
    network = 'AK,AT,YV,PS,AV,IU,II,XZ,XM'
    channel = 'BH*'

# ERROR EXAMPLE [obspy]
# GOAL: to get all waveforms for this event
# PROBLEM: No waveforms are returned -- perhaps related to the tbefore request
# ERROR MESSAGE: ValueError: The length of the input vector x must be at least padlen, which is 39.
if iex == 2:
    otime = obspy.UTCDateTime("2016-01-24T10:30:29.557")
    min_dist = 0 
    max_dist = 500
    before = 42          # Crashes
    #before = 41         # works fine
    after = 600
    #network = 'AK,AT,CN,II,IU,TA,XM,XV,XZ,ZE'
    network = 'AV'       # Crashes when -  before = 42; Works fine when - before = 41
    channel = 'BH*'

# ERROR EXAMPLE [obspy]
# GOAL: to be able to set network = '*'
# PROBLEM: is a particular network is requested (either explicitly or within *), no waveforms are returned
# ERROR MESSAGE: NotImplementedError: ResponseListResponseStage not yet implemented due to missing example data. Please contact the developers with a test data set (waveforms and StationXML metadata).
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
# GOAL: to be able to get waveforms that are embargoed
# PROBLEM: cannot get ZE waveforms (SALMON)
if iex == 4:
    otime = obspy.UTCDateTime("2016-01-24T10:30:29.557")
    min_dist = 0 
    max_dist = 700
    before = 100
    after = 600
    network = 'AK,ZE'   # no ZE waveforms returned
    channel = 'BH*,HH*'

# ERROR EXAMPLE [UAF]
# GOAL: to get multiple waveforms for same station, different location
# PROBLEM: file name does not have location code
if iex == 5:
    otime = obspy.UTCDateTime("2016-01-24T10:30:29.557")
    min_dist = 0 
    max_dist = 500
    before = 100
    after = 600
    network = 'II'
    channel = 'BH*'

# fetch and process waveforms
if idb == 1:
    # import functions to access waveforms
    import getwaveform_iris
    from obspy.clients.fdsn import Client
    client = Client("IRIS")
    # will only work for events in the 'IRIS' catalog
    # (future: for Alaska events, read the AEC catalog)
    cat = client.get_events(starttime = otime-10, endtime = otime+10)

    # Extract waveforms, IRIS
    getwaveform_iris.run_get_waveform(event = cat[0], min_dist = min_dist, max_dist = max_dist, 
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
