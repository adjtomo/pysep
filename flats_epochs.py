#!/usr/bin/python

####!/opt/antelope/python2.7.6/bin/python
# this code gets the onage, offsets, and install times.  Noise data was used to determine a rough estimate of time.  Installation times found with the webpage: http://ds.iris.edu/mda/XV/F3TN/--/HHE

from get_soh import get_iris_soh
from get_soh import get_aec_soh
import matplotlib.pyplot as plt
import numpy
import matplotlib.dates as dates
import datetime
import obspy

import obspy
from obspy.core import *
from obspy.clients.fdsn import Client


print(sys.version)


# Set parameters
akstcorr = -9
net = "XV"
loc = "*"
chan = "HHZ"
stationlist = ["F2TN","FPAP","FNN2","F4TN"]
time1 = [UTCDateTime(2017,5,5,12,0,0),UTCDateTime(2016,5,8,18,0,0),UTCDateTime(2016,8,1,20,0,0),UTCDateTime(2016,8,2,6,0,0)]
time2 = [UTCDateTime(2017,5,5,13,0,0),UTCDateTime(2016,5,8,19,0,0),UTCDateTime(2016,8,2,21,0,0),UTCDateTime(2016,8,2,7,0,0)]

count = 0

for sta in stationlist:
    # get data
    client = Client("IRIS")
    try:
        st = client.get_waveforms(net, sta, loc,chan, time1[count], time2[count], attach_response = False)
    except:
        elog.die("No data available for station %s, cha %s. Exiting..." % (sta, cha))
    count = count+1    
    #print(st)
    print("outage for ",sta,".",chan,":",st[0].stats.endtime+akstcorr*3600,"AKST")
# make plots
    #fig = plt.figure(figsize=(8.5,11))
    #plt.subplot(111)
    #plt.title(sta + ' ' + wft[0].strftime('%Y/%m/%dT%H:%M:%S.%f') + ' - ' + wft[len(wft)-1].strftime('%Y/%m/%dT%H:%M:%S.%f') + "\n HHZ ")
    #plt.plot_date(wft,wf,'b-')

    #ymax = max(abs(i) for i in wf)

    #plt.ylim(-ymax,ymax)

    #fig.autofmt_xdate()
    #plt.show()
    #print(wft[len(wft)-1].strftime('%Y/%m/%dT%H:%M:%S.%f'))
    
#stationlist = ["F2TN","FPAP","FNN2","F4TN"]
time1 = [UTCDateTime(2017,6,6,18,0,0),UTCDateTime(2016,6,6,21,0,0),UTCDateTime(2016,9,13,19,0,0),UTCDateTime(2016,9,16,21,0,0)]
time2 = [UTCDateTime(2017,6,6,19,0,0),UTCDateTime(2016,6,6,22,0,0),UTCDateTime(2016,9,13,20,0,0),UTCDateTime(2016,9,16,22,0,0)]

count = 0

for sta in stationlist:
    # get data
    client = Client("IRIS")
    try:
        st = client.get_waveforms(net, sta, loc, chan, time1[count], time2[count], attach_response = False)
    except:
        elog.die("No data available for station %s, cha %s. Exiting..." % (sta, cha))
    count = count+1    
    #print(st)
    print("onage for ",sta,".",chan,":",st[0].stats.starttime+akstcorr*3600,"AKST")

stationlist = ["F1TN","F2TN","F3TN","F4TN","F5MN","F6TP","F7TV","F8KN","FAPT","FNN1","FNN2","FPAP","FTGH"]
time1 = [UTCDateTime("2015-08-31T17:00:00"),UTCDateTime("2015-08-27T21:00:00"),UTCDateTime("2014-09-29T15:00:00"),UTCDateTime("2015-08-31T18:00:00"),UTCDateTime("2015-09-02T00:00:00"),UTCDateTime("2015-09-04T00:00:00"),UTCDateTime("2015-09-03T07:00:00"),UTCDateTime("2015-09-03T00:00:00"),UTCDateTime("2015-08-12T18:00:00"),UTCDateTime("2015-08-14T02:00:00"),UTCDateTime("2015-08-13T22:00:00"),UTCDateTime("2014-09-03T23:00:00"),UTCDateTime("2015-08-12T22:00:00")]

time2 = [UTCDateTime("2015-09-02T18:00:00"),UTCDateTime("2015-08-28T22:00:00"),UTCDateTime("2014-09-30T17:00:00"),UTCDateTime("2015-08-31T19:00:00"),UTCDateTime("2015-09-02T01:00:00"),UTCDateTime("2015-09-04T01:00:00"),UTCDateTime("2015-09-03T08:00:00"),UTCDateTime("2015-09-03T01:00:00"),UTCDateTime("2015-08-12T19:00:00"),UTCDateTime("2015-08-14T03:00:00"),UTCDateTime("2015-08-13T23:00:00"),UTCDateTime("2014-09-04T01:00:00"),UTCDateTime("2015-08-12T23:00:00")]

count = 0

chan2 = "HHZ"
for sta in stationlist:
    # get data
    client = Client("IRIS")
    try:
        st = client.get_waveforms(net, sta, loc, chan2, time1[count], time2[count], attach_response = False)
    except:
        elog.die("No data available for station %s, cha %s. Exiting..." % (sta, cha))
    count = count+1    
    #print(st)
    print("starttime for ",sta,".",chan2,":",st[0].stats.starttime+akstcorr*3600,"AKST")
