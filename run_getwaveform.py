# To run this script:
# > python run_getwaveform.py
import obspy
import copy

# EXAMPLES:
iex = 2    # Choose example to run 
idb = 1    # default: =1-IRIS; =2-AEC; =3-LLNL

# default settings
rotate = True
output_cap_weight_file = True
remove_response = True
detrend = True
demean = True
output_event_info = True
pre_filt = (0.005, 0.006, 5.0, 10.0)    # Why is this needed?
resamp_freq = 20.0                      # =0 for no resampling
scale_factor = 10.0**2                    # =10.0**2 (convert m/s to cm/s)

# Iniskin earthquake
# This produces NO waveforms and exits with this error:
# ValueError: The length of the input vector x must be at least padlen, which is 39.
if iex == 1:
    otime = obspy.UTCDateTime("2016-01-24T10:30:29.557")
    min_dist = 0 
    max_dist = 500
    #before = 42   # Crashes
    before = 41    # works fine
    after = 600
    # Vipul: figure out which network it is
    network = 'AK,AT,AV,CN,II,IU,TA,XM,XV,XZ,ZE'
    channel = 'BH*'

# SilwalTape2016 example event (Anchorage)
if iex == 2:
    otime = obspy.UTCDateTime("2009-04-07T20:12:55")
    min_dist = 0 
    max_dist = 500
    before = 100
    after = 300
    network = 'AK,AT,YV,PS,AV,IU,II,XZ,XM'
    channel = 'BH*'

# Example for showing the script crashes when requesting for IM network stations
# NotImplementedError: ResponseListResponseStage not yet implemented due to missing example data. Please contact the developers with a test data set (waveforms and StationXML metadata).
if iex == 3:
    otime = obspy.UTCDateTime("2009-04-07T20:12:55")
    min_dist = 300 
    max_dist = 500
    before = 100
    after = 300
    #network = 'AK,IM' # Crashes (Error because of the IM network)
    network = 'AK'   # works fine 
    channel = 'BH*'

# Iniskin earthquake -- trying to get waveforms from ZE (SALMON), which is embargoed
if iex == 4:
    otime = obspy.UTCDateTime("2016-01-24T10:30:29.557")
    min_dist = 0 
    max_dist = 500
    before = 100
    after = 600
    network = 'AK,ZE'   # no ZE waveforms returned
    channel = '*H*'

# Iniskin earthquake -- trying to get KDAK waveforms (different location codes)
# note that we should have B = -100 but that is not the case
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
    c = Client("IRIS")
    # will only work for events in the 'IRIS' catalog
    # (future: for Alaska events, read the AEC catalog)
    cat = c.get_events(starttime = otime-10, endtime = otime+10)

# Extract waveforms 
getwaveform_iris.run_get_waveform(ev = cat[0], min_dist = min_dist, max_dist = max_dist, 
                                  before = before, after = after, network = network, channel = channel, 
                                  samp_freq = resamp_freq, ifrotate = rotate,
                                  ifCapInp = output_cap_weight_file, ifRemoveResponse = remove_response,
                                  ifDetrend = detrend, ifDemean = demean, ifEvInfo = output_event_info, 
                                  scale_factor = scale_factor, pre_filt = pre_filt)
