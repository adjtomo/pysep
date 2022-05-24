"""
Read utilities for Pysep
"""
import re
import os
import yaml
import numpy as np

from obspy.geodetics import gps2dist_azimuth


def read_yaml(fid):
    """
    Read a YAML file and return a dictionary
    """
    # work around PyYAML bugs
    yaml.SafeLoader.add_implicit_resolver(
        u'tag:yaml.org,2002:float',
        re.compile(u'''^(?:
         [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
        |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
        |\\.[0-9_]+(?:[eE][-+][0-9]+)?
        |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
        |[-+]?\\.(?:inf|Inf|INF)
        |\\.(?:nan|NaN|NAN))$''', re.X),
        list(u'-+0123456789.'))

    with open(fid, 'r') as f:
        config = yaml.safe_load(f)

    # Replace 'None' and 'inf' values to match expectations
    for key, val in config.items():
        if val == "None":
            config[key] = None
        if val == "inf":
            config[key] = np.inf

    return config


def write_stations_file(inv, event, fid="./stations_list.txt"):
    """
    Write a list of station codes, distances, etc.
    """
    event_latitude = event.preferred_origin().latitude
    event_longitude = event.preferred_origin().longitude

    with open(fid, "w") as f:
        for net in inv:
            for sta in net:
                dist_m, az, baz = gps2dist_azimuth(lat1=event_latitude,
                                                   lon1=event_longitude,
                                                   lat2=sta.latitude,
                                                   lon2=sta.longitude
                                                   )
                dist_km = dist_m * 1E-3
                f.write(f"{sta.code:<6} {net.code:<2} "
                        f"{sta.latitude:9.4f} {sta.longitude:9.4f} "
                        f"{dist_km:8.3f} {az:6.2f}\n")
