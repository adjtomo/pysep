#!/usr/bin/python

####!/opt/antelope/python2.7.6/bin/python
# units for state of health can be found here:
# https://www.passcal.nmt.edu/content/data-archiving/documentation/passive-source

#import os
#import sys

#import signal

#signal.signal(signal.SIGINT, signal.SIG_DFL)

#sys.path.append(os.environ['ANTELOPE'] + "/data/python")
#import antelope


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
#
#from obspy.core import read, Stream, UTCDateTime
#from obspy_ext.antelope.utils import add_antelope_path
#from obspy_ext.antelope.dbobjects import Dbrecord, DbrecordList 


print(sys.version)


# Set parameters
net = "XV"
stationlist = ["F1TN", "F2TN", "F3TN","F4TN","F5MN","F6TP","F7TV","F8KN","FTGH","FNN1","FNN2","FAPT","FPAP"]

#sta = "F2TN"


time1 = UTCDateTime(2017,1,1,0,0,0)
time2 = UTCDateTime(2017,4,6,0,0,0)

for sta in stationlist:
    # get data
    vm0,vm0t = get_iris_soh(net,sta,"VM0",time1,time2)

    vm1,vm1t = get_iris_soh(net,sta,"VM1",time1,time2)

    vm2,vm2t = get_iris_soh(net,sta,"VM2",time1,time2)

    vki,vkit = get_iris_soh(net,sta,"VKI",time1,time2)

    veptemp,vept = get_iris_soh(net,sta,"VEP",time1,time2)
    vep = []
    vep = [a*0.15 for a in veptemp]

    vec,vect = get_iris_soh(net,sta,"VEC",time1,time2)

# make plots
    fig = plt.figure(figsize=(8.5,11))
    plt.subplot(411)
    plt.title(sta + ' ' + vkit[0].strftime('%Y/%m/%d') + ' - ' + vkit[len(vkit)-1].strftime('%Y/%m/%d') + "\n VKI temperature ")
    plt.plot_date(vkit,vki,'b-')
    plt.ylim(-10,50)
    plt.ylabel("Celsius")

    plt.subplot(412)
    plt.title("VEP main system voltage")
    plt.plot_date(vept,vep,'b-')
    plt.ylim(8,17)
    plt.ylabel("Volts")

    plt.subplot(413)
    plt.title("VEC main system current")
    plt.plot_date(vect,vec,'b-')
    plt.ylim(50,100)
    plt.ylabel("mA")

    plt.subplot(414)
    plt.title("Mass Positions")
    plt.plot_date(vm0t,vm0,'g-',label="VM0")
    plt.plot_date(vm1t,vm1,'r-',label="VM1")
    plt.plot_date(vm2t,vm2,'b-',label="VM2")
    plt.ylim(-35,35)
    plt.legend(loc=2,fontsize = 10)

#pp = max(abs(i) for i in vm0)
#pp2 = max(abs(j) for j in vm1)
#pp3 = max(abs(k) for k in vm2)
#ymax = max(pp,pp2,pp3)

#plt.ylim(-ymax,ymax)

    fig.autofmt_xdate()
    fig.savefig("/home/ksmith/currentSOH/" + sta + "soh.pdf", dpi=200)

    #plt.show()



'''
cha = "VM0"
vmo_aec = get_aec_soh(sta,cha,2016,8,1,0,0,0)
'''

'''
# get groundmotion waveforms
da2 = 17
ho2 = 11
mi2 = 13
oadd = 60*60
w1,t1 = get_iris_soh(net,sta,"HHZ",ye,mo,da2,ho2,mi2,se,oadd)

w2,t2 = get_iris_soh(net,sta,"HHN",ye,mo,da2,ho2,mi2,se,oadd)

w3,t3 = get_iris_soh(net,sta,"HHE",ye,mo,da2,ho2,mi2,se,oadd)

fig = plt.figure(3)
plt.subplot(311)
plt.title(sta + " HHZ")
plt.plot_date(t1,w1,'b-')
plt.subplot(312)
plt.title("HHN")
plt.plot_date(t2,w2,'b-')
plt.subplot(313)
plt.title("HHE")
plt.plot_date(t3,w3,'b-')
fig.autofmt_xdate()
fig.savefig(sta + "wavesignal.png")

plt.show()
'''
