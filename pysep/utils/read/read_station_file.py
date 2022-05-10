#! /usr/bin/env python
# -*- coding: utf-8 -*-

def read_station_file(filename):
    f = open(filename, "r")
    lines = f.readlines()[1:]
    f.close()
    netwk = [];stn_name = [];
    # Before 2019/7/17 it was "True" but now it will be False 
    # since we do not want to redownload WFs everytime stations 
    # selection changes, just use matlab, etc to omit stations
    use_ones_stations = False; 
    for line in lines:
        line_elements = line.split()
        if (line_elements[8] == '1') and use_ones_stations:
            netwk.append(line_elements[0])
            stn_name.append(line_elements[1])
        elif not use_ones_stations:
            netwk.append(line_elements[0])
            stn_name.append(line_elements[1])
    return netwk,stn_name
