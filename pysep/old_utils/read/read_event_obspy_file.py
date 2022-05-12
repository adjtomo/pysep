#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to grab source parameters from files generated from  write_seis_obspy.m
"""

def read_events_obspy_file(filename):
    f = open(filename, "r")
    lines = f.readlines()
    f.close()
    count = 0
    eids = [];otimes = [];elons = [];elats = [];edeps = [];emags = []

    for line in lines:
        line_elements = line.split()
        eids.append(line_elements[1])
        otimes.append(line_elements[2])
        elons.append(line_elements[3])
        elats.append(line_elements[4])
        edeps.append(line_elements[5])
        emags.append(line_elements[6])
        print(line_elements[1],line_elements[2],line_elements[3],
              line_elements[4],line_elements[5],line_elements[6])
    return eids,otimes,elons,elats,edeps,emags
