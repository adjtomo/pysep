#!/opt/antelope/python2.7.6/bin/python

'''
This tool allows us see how the P wave arrival changes with distance for a fixed depth

To compile:
$ make

To run:
$ ./Parrival 100

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

# a way for us to add in arguments
parser = argparse.ArgumentParser(description='A source depth of an event in km')
parser.add_argument('dep', action='store', type=float, help='fixed source depth in km')
args = parser.parse_args()
sdepth = args.dep

# choose a velocity model
model = TauPyModel(model="iasp91")

# Trying to be realistic
if(sdepth < 0 or sdepth > 6371):
    print("source depth must be between 0 and 6371!")
    exit(0)



dists = np.linspace(0,180,91)
mytimes = []
keepdists = []

# Calculate arrival times
for x in range(0,len(dists)):
    try:
        arrivals = model.get_travel_times(source_depth_in_km=sdepth,distance_in_degree=dists[x],phase_list=["P"])
    # if there is no output then we will get an error
        mytimes.append(arrivals[0].time)
        keepdists.append(dists[x])
    except:
        print("function does not work with distance", dists[x], " degrees" )

# Plot the results
plt.plot(keepdists,mytimes,'ro')
plt.ylabel('P Travel Time, s')
plt.xlabel('Epicentral Distance, deg')
plt.title("Source depth %s km"  % sdepth)
plt.show()
