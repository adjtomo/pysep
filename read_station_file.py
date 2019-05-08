#! /usr/bin/env python
# -*- coding: utf-8 -*-

def read_station_file(filename):
    f = open(filename, "r")
    lines = f.readlines()
    f.close()
    netwk = [];stn_name = [];

    for line in lines:
        line_elements = line.split()
        if line_elements[7] == '1':
            netwk.append(line_elements[0])
            stn_name.append(line_elements[1])
    return netwk,stn_name
