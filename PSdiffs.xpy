#!/opt/antelope/python2.7.6/bin/python

'''
This tool allows us see how the P-S differential time changes in source depth for a fixed distance

To compile:
$ make

To run:
$ ./PSdiffs 20

-- Kyle Smith 

Created on August 8, 2016
'''

import re
from sys import exit
import argparse
import math as math
import numpy as np
import matplotlib.pyplot as plt

import obspy
from obspy.taup import TauPyModel

# a way for this function to add arguments
parser = argparse.ArgumentParser(description='Epicentral distance of an event in degrees')
parser.add_argument('edist', action='store', type=float, help='fixed epicentral distance in degrees')
args = parser.parse_args()
epidist = args.edist

# choose the model you want
model = TauPyModel(model="ak135")

# We gotta be realistic here
if(epidist < 0 or epidist> 180):
    print("Epicenter distance must be between 0 and 180!")
    exit(0)

sdepths = np.linspace(0,200,51)
#sdepths = [100]
myptimes = []
mystimes = []
keepdepths = []
mydiffs = []

# get arrival times and find the PS differential
for x in range(0,len(sdepths)):
    try:
       parrivals = model.get_travel_times(source_depth_in_km=sdepths[x],distance_in_degree=epidist,phase_list=["P"])
       sarrivals = model.get_travel_times(source_depth_in_km=sdepths[x],distance_in_degree=epidist,phase_list=["S"])
       #print(sdepths[x])
       #print(parrivals[0].name)
       #print(sarrivals[0].time)
       #print(sarrivals[0].time-parrivals[0].time)
       myptimes.append(parrivals[0].time)
       mystimes.append(sarrivals[0].time)
       mydiffs.append(sarrivals[0].time-parrivals[0].time)
       keepdepths.append(sdepths[x])
    except:
        print("function does not work with depth", sdepths[x], "km" )

# Plot the results
#plt.plot(keepdepths,myptimes,'ro')
plt.plot(keepdepths,mydiffs,'ro')
plt.ylabel('P-S Time, s')
plt.xlabel('Source Depth, km')
plt.title("P-S differential time for epicentral distance %s degrees"  % epidist)
plt.show()


