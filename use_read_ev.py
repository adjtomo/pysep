#! /usr/bin/env python
# -*- coding: utf-8 -*-

# example script for looping over events in a text file

import read_event_obspy_file as reof
import obspy
import os

events_file = "/home/ksmith/REPOSITORIES/manuscripts/kyle/papers/basinamp/data/basinamp_obspy.txt"

eids,otimes,elons,elats,edeps,emags = reof.read_events_obspy_file(events_file)

for xx in range(0,3):
    otime = obspy.UTCDateTime(otimes[xx])
    elat = elats[xx]
    elon = elons[xx]
    edep = edeps[xx]
    emag = emags[xx]
    eid = eids[xx]
   # print(emag)
    os.system('ls | head -12 > ' + eid + '_last_2git_commits.txt')

    

