import obspy
from getwaveform import *

def get_ev_info(ev_info,iex):

    # center at F3TN
    ev_info.rlat = 64.7716  
    ev_info.rlon = -149.1465 
    ev_info.scale_factor = 1
    ev_info.min_dist = 0
    ev_info.max_dist = 150
    ev_info.tbefore_sec = 100
    ev_info.tafter_sec = 200 
    ev_info.network = 'AK,AT,II,IU,US,XM,XV,XZ,TA'
    ev_info.channel = 'BH?,HH?'
    ev_info.use_catalog = 0
    ev_info.station = 'I23K'

    if iex == 0:
        # AEC parameters
        ev_info.otime = obspy.UTCDateTime("2015-11-20T10:53:48.168") 
        ev_info.elat = 64.6210
        ev_info.elon = -149.4024
        ev_info.edep = 17113.4
        ev_info.emag = 2.67
        ev_info.rtime = ev_info.otime

    if iex == 1:
        # AEC parameters
        ev_info.otime = obspy.UTCDateTime("2015-10-22T13:16:15.794") 
        ev_info.elat = 64.7334
        ev_info.elon = -149.0388
        ev_info.edep = 18830.2
        ev_info.emag = 2.74  
        ev_info.rtime = ev_info.otime

    if iex == 2:
        # AEC parameters
        ev_info.otime = obspy.UTCDateTime("2014-12-13T15:47:31.423") 
        ev_info.elat = 64.4325
        ev_info.elon = -149.3840
        ev_info.edep = 12431.1
        ev_info.emag = 3.25  
        ev_info.rtime = ev_info.otime
        
    if iex == 3:
        # AEC source parameters
        ev_info.otime = obspy.UTCDateTime("2015-10-31T02:56:35.572") 
        ev_info.elat = 64.4285
        ev_info.elon = -149.6969
        ev_info.edep = 23852.1
        ev_info.emag = 3.47
        ev_info.rtime = ev_info.otime

    if iex == 4:
        # AEC source parameters
        ev_info.otime = obspy.UTCDateTime("2015-11-06T01:20:12.712") 
        ev_info.elat = 64.7552
        ev_info.elon = -151.3103
        ev_info.edep = 1502.1
        ev_info.emag = 3.35
        ev_info.rtime = ev_info.otime

    if iex == 5:
        # AEC source parameters
        ev_info.otime = obspy.UTCDateTime("2014-10-21T00:36:58.333") 
        ev_info.elat = 65.1489
        ev_info.elon = -149.0413
        ev_info.edep = 13134.8
        ev_info.emag = 4.90
        ev_info.rtime = ev_info.otime

    if iex == 6:
        # AEC source parameters
        ev_info.otime = obspy.UTCDateTime("2014-10-23T16:30:23.968") 
        ev_info.elat = 65.1644
        ev_info.elon = -149.0523
        ev_info.edep = 200665
        ev_info.emag = 5.00
        ev_info.rtime = ev_info.otime

    if iex == 7:
        # AEC source parameters
        ev_info.otime = obspy.UTCDateTime("2016-01-14T19:04:10.727") 
        ev_info.elat = 64.6827
        ev_info.elon = -149.2479
        ev_info.edep = 22663.7
        ev_info.emag = 3.80
        ev_info.rtime = ev_info.otime

    if iex == 8:
        # AEC source parameters
        ev_info.otime = obspy.UTCDateTime("2015-09-28T11:50:12.949") 
        ev_info.elat = 64.7148
        ev_info.elon = -148.9769
        ev_info.edep = 15112.7
        ev_info.emag = 2.91
        ev_info.rtime = ev_info.otime

    if iex == 9:
        # Big Minto Event
        # AEC source parameters
        ev_info.otime = obspy.UTCDateTime("2016-11-06T09:29:10.579") 
        ev_info.elat = 64.1639
        ev_info.elon = -150.0626
        ev_info.edep = 23190.0
        ev_info.emag = 4.00
        ev_info.rtime = ev_info.otime

    if iex == 10:
        # Big Minto Event
        # AEC source parameters
        ev_info.otime = obspy.UTCDateTime("2016-12-08T10:18:13.868") 
        ev_info.elat = 64.1937
        ev_info.elon = -150.0376
        ev_info.edep = 24522.1
        ev_info.emag = 4.30
        ev_info.rtime = ev_info.otime

    if iex == 11:
        # Iniskin Event
        ev_info.otime = obspy.UTCDateTime("2016-01-24T10:30:29.557") 
        ev_info.elat = 59.6204
        ev_info.elon = -153.3392
        ev_info.edep = 125645.3
        ev_info.emag = 7.10
        ev_info.rtime = ev_info.otime
        ev_info.tafter_sec = 600

    if iex == 12:
        # Chile Event
        ev_info.otime = obspy.UTCDateTime("2015-09-16T22:54:33.000")
        ev_info.elat = -31.5695
        ev_info.elon = -71.6543
        ev_info.edep = 22400.0
        ev_info.emag = 8.30
        ev_info.rtime = obspy.UTCDateTime("2015-09-16T23:09:15.000")
        ev_info.tafter_sec = 200

    if iex == 13:
        # Mariana Event observed by FLATS
        ev_info.otime = obspy.UTCDateTime("2016-07-29T21:18:26.000")
        ev_info.elat = 18.5439
        ev_info.elon = 145.541
        ev_info.edep = 207620.0
        ev_info.emag = 7.7
        ev_info.rtime = obspy.UTCDateTime("2016-07-29T21:28:19.000")
        
    if iex == 14:
        ev_info.otime = obspy.UTCDateTime("2017-01-31T09:38:37.576")
        #otime = obspy.UTCDateTime("2017-01-31T09:38:37.000")
        ev_info.elat = 63.0817
        ev_info.elon = -150.9427
        ev_info.edep = 132900
        ev_info.emag = 5.2
        ev_info.rtime = ev_info.otime

    if iex == 15:
        ev_info.otime = obspy.UTCDateTime("2017-04-29T11:15:48.000")
        ev_info.elat = 63.1296
        ev_info.elon = -151.1517
        ev_info.edep = 10800
        ev_info.emag = 5.2
        ev_info.rtime = ev_info.otime

    if iex == 16:
        ev_info.otime = obspy.UTCDateTime("2017-06-28T12:58:52")
        ev_info.elat = 64.7569
        ev_info.elon = -148.8883
        ev_info.edep = 18000
        ev_info.emag = 3.5
        ev_info.rtime = ev_info.otime

    return(ev_info)
