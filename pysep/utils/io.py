"""
Read utilities for Pysep
"""
import re
import yaml
import numpy as np

from obspy import UTCDateTime
from obspy.geodetics import gps2dist_azimuth


def read_event_file(fid):
    """
    Read an event input file which is just a text file that should contain
    information on different events and their hypocenters

    :type fid: str
    :param fid: event input file
    :rtype: list of dict
    :return: parsed in event information
    """
    list_out = []
    with open(fid, "r") as f:
        lines = f.readlines()
    for line in lines:
        # Commented out lines are skipped
        if line.strip().startswith("#"):
            continue
        origin_time, longitude, latitude, depth_km, mag = line.strip().split()
        dict_out = {"origin_time": UTCDateTime(origin_time),
                    "event_longitude": float(longitude),
                    "event_latitude": float(latitude),
                    "event_depth_km": float(depth_km),
                    "event_magnitude": float(mag)
                    }
        list_out.append(dict_out)

    return list_out


def read_yaml(fid):
    """
    Read a YAML file and return a dictionary

    :type fid: str
    :param fid: YAML file to read from
    :rtype: dict
    :return: YAML keys and variables in a dictionary
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

    :type event: obspy.core.event.Event
    :param event: optional user-provided event object which will force a
        skip over QuakeML/event searching
    :type inv: obspy.core.inventory.Inventory
    :param inv: optional user-provided inventory object which will force a
        skip over StationXML/inventory searching
    :type fid: str
    :param fid: name of the file to write to. defaults to ./stations_list.txt
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
