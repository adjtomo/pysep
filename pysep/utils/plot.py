#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pysep plotting utilities
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from obspy.core.event import Catalog
from obspy.taup import TauPyModel


def plot_source_receiver_map(inv, event, fid="./station_map.png"):
    """
    Simple utility to plot station and event on a single map
    """
    # Temp change the size of all text objects to get text labels small
    # and make the event size a bit smaller
    with plt.rc_context({"font.size": 5, "lines.markersize": 4}):
        inv.plot(projection="local", resolution="i", label=True, show=False,
                 size=10, color="g", method="cartopy")
        # !!! fig=plt.gca() kwarg is a workaround for bug in
        # !!! ObsPy==1.3.0 (obspy.imaging.maps._plot_cartopy_into_axes)
        Catalog(events=[event]).plot(method="cartopy", outfile=fid,
                                     fig=plt.gca(), show=False)


def plot_phases():
    """
    Tool to plot phases

    TODO fix this up

    Kyle Smith
    January 29,2018
    :return:
    """
    phases=["S"]
    phases2=["s"]
    sourcedepth = 100
    taup_model = "ak135"
    dists = list(np.arange(0,180,1))
    model = TauPyModel(model=taup_model)
    Phase1arrivals = []
    Phase2arrivals = []
    Phase1_IA = []
    Phase2_IA = []
    for dist_deg in dists[:]:
        temparr = model.get_travel_times(source_depth_in_km=sourcedepth,distance_in_degree=dist_deg,phase_list=[phases[0]])
        temparr2 = model.get_travel_times(source_depth_in_km=sourcedepth,distance_in_degree=dist_deg,phase_list=[phases2[0]])

        # there may be something wrong with obpsy and plotting since the error is that the arrivals does not have the attribute for plot_rays(), see Vipul's email for another way to plot.
        try:
            somearray = model.get_ray_paths(source_depth_in_km=sourcedepth,distance_in_degree=dist_deg,phase_list=[phases[0]])

        except:
            print('no somearray!')

        if len(temparr)==0:
            Phase1arrivals.append(math.nan)
            Phase1_IA.append(math.nan)
        else:
            Phase1arrivals.append(temparr[0].time)
            Phase1_IA.append(temparr[0].incident_angle)
        if len(temparr2)==0:
            Phase2arrivals.append(math.nan)
            Phase2_IA.append(math.nan)
        else:
            Phase2arrivals.append(temparr2[0].time)
            Phase2_IA.append(temparr2[0].incident_angle)
    dkm = dists
    #dkm = [a*111 for a in dists]
    f1 = plt.figure(1)
    L2, = plt.plot(dkm,Phase2arrivals,'ro',label=phases2[0])
    L1, = plt.plot(dkm,Phase1arrivals,'b*',label=phases[0])
    plt.xlim(0, 180)
    plt.legend(handles=[L1,L2])
    plt.ylabel("Time after EQ origin time (s)")
    plt.xlabel("Source-Receiver Distance (deg)")
    plt.title("Source depth:" + str(sourcedepth) + " km")
    f1.show()

    f2 = plt.figure(2)
    L2, = plt.plot(dkm,Phase2_IA,'ro',label=phases2[0])
    L1, = plt.plot(dkm,Phase1_IA,'b*',label=phases[0])
    plt.xlim(0, 180)
    plt.legend(handles=[L1,L2])
    plt.ylabel("Incident Angle(deg)")
    plt.xlabel("Source-Receiver Distance(deg)")
    plt.title("Source depth:" + str(sourcedepth) + " km")
    f2.show()

    input()
