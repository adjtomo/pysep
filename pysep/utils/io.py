"""
Read utilities for Pysep
"""
import os
import re
import yaml
import numpy as np

from obspy import UTCDateTime, Stream, Trace, read_events, Inventory
from obspy.core.inventory.network import Network
from obspy.core.inventory.station import Station
from obspy.geodetics import gps2dist_azimuth

from pysep.utils.cap_sac import append_sac_headers


def read_synthetics(fid, cmtsolution, stations, location="", precision=4):
    """
    Specfem3D outputs seismograms to ASCII (.sem?) files. Converts SPECFEM
    .sem? files into Stream objects with the correct header
    information.

    .. note::
        if origintime is None, the default start time is 1970-01-01T00:00:00,
        this is fine for quick-looking the data but it is recommended that an
        actual starttime is given.

    .. note ::
        Tested up to Specfem3D Cartesian git version 6895e2f7

    :type fid: str
    :param fid: path of the given ascii file
    :type cmtsolution: str
    :param cmtsolution: CMTSOLUTION file defining the event which generated
        the synthetics. Used to grab event information
    :type stations: str
    :param stations: STATIONS file defining the station locations for the
        SPECFEM generated synthetics
    :type location: str
    :param location: location value for a given station/component
    :rtype st: obspy.Stream.stream
    :return st: stream containing header and data info taken from ascii file
    """
    # This was tested up to SPECFEM3D Cartesian git version 6895e2f7
    try:
        times = np.loadtxt(fname=fid, usecols=0)
        data = np.loadtxt(fname=fid, usecols=1)

    # At some point in 2018, the Specfem developers changed how the ascii files
    # were formatted from two columns to comma separated values, and repeat
    # values represented as 2*value_float where value_float represents the data
    # value as a float
    except ValueError:
        times, data = [], []
        with open(fid, 'r') as f:
            lines = f.readlines()
        for line in lines:
            try:
                time_, data_ = line.strip().split(',')
            except ValueError:
                if "*" in line:
                    time_ = data_ = line.split('*')[-1]
                else:
                    raise ValueError
            times.append(float(time_))
            data.append(float(data_))

        times = np.array(times)
        data = np.array(data)

    # We assume that dt is constant after 'precision' decimal points
    delta = round(times[1] - times[0], precision)

    # Get metadata information from CMTSOLUTION and STATIONS files
    event = read_events(cmtsolution)[0]
    inv = read_stations(stations)
    origintime = event.preferred_origin().time

    # Honor that Specfem doesn't start exactly on 0
    origintime += times[0]

    # Write out the header information
    try:
        # SPECFEM2D/SPECFEM3D_Cartesian style name format, e.g., NZ.BFZ.BXE.semd
        net, sta, cha, fmt = os.path.basename(fid).split(".")
    except ValueError:
        # SPECFEM3D_Globe style name format, e.g., TA.O20K.BXR.sem.ascii
        net, sta, cha, fmt, suffix = os.path.basename(fid).split(".")

    stats = {"network": net, "station": sta, "location": location,
             "channel": cha, "starttime": origintime, "npts": len(data),
             "delta": delta, "mseed": {"dataquality": 'D'},
             "time_offset": times[0], "format": fmt
             }
    st = Stream([Trace(data=data, header=stats)])
    st = append_sac_headers(st, event, inv)

    return st


def read_stations(path_to_stations):
    """
    Convert a Specfem3D STATIONS file into an ObsPy Inventory object.

    Specfem3D STATION files contain no channel or location information, so
    the inventory can only go down to the station level.

    .. note::
        This is copied verbatim from Pyatoa.utils.read.read_stations()

    .. note::
        This assumes a row structure for the station file is
        STA, NET, LAT [deg], LON [deg], ELEVATION [m], BURIAL [m]

    :type path_to_stations: str
    :param path_to_stations: the path to the STATIONS file that is associated
        with the Specfem3D DATA directory
    :rtype: obspy.core.inventory.Inventory
    :return: a station-level Inventory object
    """
    stations = np.loadtxt(path_to_stations, dtype="str")
    if stations.size == 0:
        return Inventory()

    # Get all the unique network names, try-except to catch when there is only
    # one station in the file
    try:
        networks = {_: [] for _ in np.unique(stations[:, 1])}
    except IndexError:
        networks = {stations[1]: []}
        stations = [stations]

    for sta in stations:
        # Parse the station information
        station_ = sta[0]
        network_ = sta[1]
        latitude_ = float(sta[2])
        longitude_ = float(sta[3])
        elevation_ = float(sta[4])
        burial_ = float(sta[5])  # burial isnt an option in ObsPy, not used

        # Create the station object, temp store in a network
        station = Station(code=station_, latitude=latitude_,
                          longitude=longitude_, elevation=elevation_,
                          creation_date=UTCDateTime()
                          )
        networks[network_].append(station)

    # Create the network objects
    list_of_networks = []
    for network, stations in networks.items():
        list_of_networks.append(Network(code=network, stations=stations))

    return Inventory(networks=list_of_networks, source="PYATOA")



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
