#!/opt/antelope/python2.7.6/bin/python

#import os
#import sys

#sys.path.append(os.environ['ANTELOPE'] + "/data/python")
#print(sys.version)

import obspy
from obspy.core import *
from obspy.clients.fdsn import Client
import datetime
from statistics import mode

#import antelope.elog as elog

def get_iris_soh(net, sta, cha, time1,time2):
    client = Client("IRIS")
    try:
        st = client.get_waveforms(net, sta, "*", cha, time1, time2, attach_response = False)
    except:
        elog.die("No data available for station %s, cha %s. Exiting..." % (sta, cha))

    # calculate statistics
    samprates = []
    dts = []
    for ii in list(range(len(st))):
        samprates.append(st[ii].stats.sampling_rate)
        dts.append(st[ii].stats.delta)

    mostpopsamp = mode(samprates)
    mostpopdt = mode(dts)

    # make mergeable
    for ii in list(range(len(st))):
        if st[ii].stats.sampling_rate != mostpopsamp:
            pctoff = (mostpopsamp-st[ii].stats.sampling_rate)/mostpopsamp*100
            if abs(pctoff) > 1:
                print("sample rates off by ", pctoff, "%")
            st[ii].stats.sampling_rate = mostpopsamp
        if st[ii].stats.delta != mostpopdt:
            st[ii].stats.delta = mostpopdt

    st.merge(method=0)
    tr = st[0]
    sttime = tr.stats.starttime
    fntime = tr.stats.endtime
    numpoints = tr.stats.npts
    date_list = [sttime + x*tr.stats.delta for x in range(0, numpoints)]
    datetime_list = [date_list[y].datetime for y in range(0, numpoints)]
    tt = datetime_list
    return tr,tt


#try:

#except:
#    sys.exit("Problem importing elog library.  Check your Antelope installation.  Exiting...")

# Datascope libraries
#import antelope.datascope as datascope

# Stock libraries
#import antelope.stock as stock


def get_aec_soh(sta, cha, yea, mon, day,hou, min, sec):
    
    # Read from state of health directory
    db_nme = '/aerun/op/db/diagnostic/diagnostic'
    # OR
    #db_nme = '/aec/db/diagnostics/diagnostics'
    #try:
    db = datascope.dbopen(db_nme, 'r')
    #except:
    #    elog.die("Problem loading the database [%s] for processing!  Exiting..." % db_nme)

    db = db.lookup(table='wfdisc')

    # The origin time
    t_or = "%s/%s/%s %s:%s:%s.00" % (mon, day, yea, hou, min, sec)

    # Convert to epoch
    t_or_st = stock.str2epoch(t_or)

    # Window around the origin time.  Time added is in units of seconds
    t1 = t_or_st + 0
    t2 = t_or_st + 60*60 

    # Load the trace
    tr = db.trloadchan(t1, t2, sta, cha)

    # Get the required data from the trace
    tr.record = 0
    (nsamp, samprate) = tr.getv('nsamp', 'samprate')
    d = array(tr.trdata()) # The actual data

    # Sanity check
    if (nsamp != len(d)):
        elog.die("For some reason the data vector has a length that does not match the number of samples!  Exiting...")

    # Convert trace object to ObsPy
    st = Stream() # Empty stream
    tr0 = Trace() # Empty ObsPy trace object 

    # Populate the trace object
    tr0.stats.network = net
    tr0.stats.station = sta
    tr0.stats.channel = cha
    tr0.stats.sampling_rate = samprate
    tr0.stats.npts = nsamp
    tr0.stats.starttime = t1
    tr0.stats.starttime = t2
    tr0.data = d

    # Append the stream with the new trace object
    st.append(tr0)

    return st
