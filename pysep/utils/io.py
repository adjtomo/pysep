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
from obspy.core.event import Event, Origin, Magnitude

from pysep import logger
from pysep.utils.mt import moment_magnitude, seismic_moment, Source
from pysep.utils.fmt import format_event_tag_legacy, channel_code
from pysep.utils.cap_sac import append_sac_headers


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


def read_sem(fid, origintime=None, source=None, stations=None, location="",
             precision=4):
    """
    Specfem3D outputs seismograms to ASCII (.sem? or .sem.ascii) files.
    Converts SPECFEM synthetics into ObsPy Stream objects with the correct header
    information.

    .. note ::
        Tested up to Specfem3D Cartesian git version 6895e2f7

    :type fid: str
    :param fid: path of the given ascii file
    :type origintime: obspy.UTCDateTime
    :param origintime: UTCDatetime object for the origintime of the event
    :type source: str
    :param source: optional SPECFEM source file (e.g., CMTSOLUTION, SOURCE)
        defining the event which generated the synthetics. Used to grab event
        information and append as SAC headers to the ObsPy Stream
    :type stations: str
    :param stations: optional STATIONS file defining the station locations for
        the SPECFEM generated synthetics, used to generate SAC headers
    :type location: str
    :param location: location value for a given station/component
    :type precision: int
    :param precision: dt precision determined by differencing two
        adjancent time steps in the underlying ascii text file.
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
    if origintime is None and source is None:
        logger.warning("no `origintime` or `event` given, setting dummy "
                       "starttime: 1970-01-01T00:00:00")
        origintime = UTCDateTime("1970-01-01T00:00:00")
    elif source:
        # !!! Add logic here to try different read options
        event = read_events(source, format="CMTSOLUTION")[0]
        origintime = event.preferred_origin().time
        logger.info(f"reading origintime from event: {origintime}")

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
             "delta": delta, "mseed": {"dataquality": 'D'}, "format": fmt
             }
    st = Stream([Trace(data=data, header=stats)])

    if source and stations:
        try:
            inv = read_stations(stations)
            st = append_sac_headers(st, event, inv)
        # Broad catch here as this is an optional step that might not always
        # work or be possible
        except Exception as e:
            logger.warning(f"could not append SAC header to trace because {e}")

    return st


def read_sem_cartesian(fid, source, stations, location="", precision=4):
    """
    Specfem2D and Specfem3D may have domains defined in a Cartesian coordinate
    system. Because of this, the read_sem() function will fail because
    the intermediate ObsPy objects and functions expect geographic coordinates.
    This function bypasses these checks with some barebones objects which
    mimic their behavior. Only used for RecSec to plot record sections.

    TODO can we combine this with cap_sac.append_sac_headers()? Currently the
        code block at the bottom (manually appending header) is redundant 
        and should use cap_sac?

    .. note::
        RecSec requires SAC header values `kevnm`, `dist`, `az`, `baz`,
        `stlo`, `stla`, `evlo`, `evla`

    :type fid: str
    :param fid: path of the given ascii file
    :type source: str
    :param source: SOURCE or CMTSOLUTION file defining the event which
        generated the synthetics. Used to grab event information.
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
    try:
        event = read_specfem3d_cmtsolution_cartesian(source)
    # Specfem2D and 3D source/cmtsolution files have different formats
    except ValueError:
        event = read_specfem2d_source(source)

    # Generate a dictionary object to store station information
    station_list = np.loadtxt(stations, dtype="str", ndmin=2)
    stations = {}
    for sta in station_list:
        # NN.SSS = {latitude, longitude}
        stations[f"{sta[1]}.{sta[0]}"] = {"stla": float(sta[2]),
                                          "stlo": float(sta[3])
                                          }

    starttime = event.preferred_origin().time

    # Honor that Specfem doesn't start exactly on 0
    starttime += times[0]

    # Write out the header information
    try:
        # SPECFEM2D/SPECFEM3D_Cartesian style name format, e.g., NZ.BFZ.BXE.semd
        net, sta, cha, fmt = os.path.basename(fid).split(".")
    except ValueError:
        # SPECFEM3D_Globe style name format, e.g., TA.O20K.BXR.sem.ascii
        net, sta, cha, fmt, suffix = os.path.basename(fid).split(".")

    stats = {"network": net, "station": sta, "location": location,
             "channel": cha, "starttime": starttime, "npts": len(data),
             "delta": delta, "mseed": {"dataquality": 'D'}, "format": fmt
             }
    st = Stream([Trace(data=data, header=stats)])

    # Manually append SAC header values here
    for tr in st:
        net_sta = f"{tr.stats.network}.{tr.stats.station}"
        stla = stations[net_sta]["stla"]
        stlo = stations[net_sta]["stlo"]
        evla = event.preferred_origin().latitude
        evlo = event.preferred_origin().longitude

        # Calculate Cartesian distance and azimuth/backazimuth
        dist_m = np.sqrt(((stlo - evlo) ** 2) + ((stla - evla) ** 2))
        azimuth = np.degrees(np.arctan2((stlo - evlo), (stla - evla))) % 360
        backazimuth = (azimuth - 180) % 360
        otime = event.preferred_origin().time
        # Only values required by RecSec
        sac_header = {
            "stla": stations[net_sta]["stla"],
            "stlo": stations[net_sta]["stlo"],
            "evla": event.preferred_origin().latitude,
            "evlo": event.preferred_origin().longitude,
            "dist": dist_m * 1E-3,
            "az": azimuth,
            "baz": backazimuth,
            "kevnm": format_event_tag_legacy(event),  # only take date code
            "nzyear": otime.year,
            "nzjday": otime.julday,
            "nzhour": otime.hour,
            "nzmin": otime.minute,
            "nzsec": otime.second,
            "nzmsec": otime.microsecond,
            }
        tr.stats.sac = sac_header

    return st


def read_specfem3d_cmtsolution_cartesian(path_to_cmtsolution):
    """
    Create a barebones ObsPy Event object from a SPECFEM3D CMTSOLUTION file with
    coordinates defined in cartesian. Required because ObsPy read_events() will
    throw a ValueError when CMTSOLUTIONS do not have geographically defined
    coordinates.
    """
    def _get_resource_id(name, res_type, tag=None):
        """
        Helper function to create consistent resource ids. From ObsPy
        """
        res_id = f"smi:local/source/{name:s}/{res_type:s}"
        if tag is not None:
            res_id += "#" + tag
        return res_id

    with open(path_to_cmtsolution, "r") as f:
        lines = f.readlines()

    # First line contains meta information about event
    _, year, month, day, hour, minute, sec, lat, lon, depth, mb, ms, _ = \
        lines[0].strip().split()

    origin_time = UTCDateTime(f"{year}-{month}-{day}T{hour}:{minute}:{sec}")

    # Remaining lines contain information on event, some useful some not
    source_dict = {}
    for line in lines[1:]:
        # Skip comments and newlines
        if line.startswith("#") or line == "\n":
            continue
        key, val = line.strip().split(":")
        # Strip trailing comments from values
        source_dict[key.strip()] = val.strip()

    origin = Origin(
        resource_id=_get_resource_id(source_dict["event name"],
                                     "origin", tag="source"),
        time=origin_time, longitude=source_dict["longorUTM"],
        latitude=source_dict["latorUTM"],
        depth=abs(float(source_dict["depth"]) * 1E3)  # units: m
    )

    magnitude = Magnitude(
        resource_id=_get_resource_id(source_dict["event name"], "magnitude"),
        mag=ms, magnitude_type="Ms", origin_id=origin.resource_id.id
    )

    event = Event(resource_id=_get_resource_id(name=source_dict["event name"],
                                               res_type="event"))
    event.origins.append(origin)
    event.magnitudes.append(magnitude)

    event.preferred_origin_id = origin.resource_id.id
    event.preferred_magnitude_id = magnitude.resource_id.id

    return event


def read_specfem2d_source(path_to_source, origin_time=None):
    """
    Create a barebones ObsPy Event object from a SPECFEM2D Source file, which
    only contains information required by Pyatoa.

    Only requires access to: event.preferred_origin(),
    event.preferred_magnitude() and event.preferred_moment_tensor().
    Moment tensor is wrapped in try-except so we only need origin and magnitude.

    Modified from:
    https://docs.obspy.org/master/_modules/obspy/io/cmtsolution/
                                            core.html#_internal_read_cmtsolution

    .. note::
        Source files do not provide origin times so we just provide an
        arbitrary value but allow user to set time
    """

    def _get_resource_id(name, res_type, tag=None):
        """Helper function to create consistent resource ids. From ObsPy"""
        res_id = f"smi:local/source/{name:s}/{res_type:s}"
        if tag is not None:
            res_id += "#" + tag
        return res_id

    # First set dummy origin time
    if origin_time is None:
        origin_time = "1970-01-01T00:00:00"
        logger.warning("no origin time set for SPECFEM2D source, setting "
                       f"dummy value: {origin_time}")

    with open(path_to_source, "r") as f:
        lines = f.readlines()

    # First line expected to be e.g.,: '## Source 1'
    source_name = lines[0].strip().split()[-1]
    source_dict = {}
    for line in lines:
        # Skip comments and newlines
        if line.startswith("#") or line == "\n":
            continue
        key, val, *_ = line.strip().split("=")
        # Strip trailing comments from values
        val = val.split("#")[0].strip()
        source_dict[key.strip()] = val.strip()

    origin = Origin(
        resource_id=_get_resource_id(source_name, "origin", tag="source"),
        time=origin_time, longitude=source_dict["xs"],
        latitude=source_dict["xs"], depth=source_dict["zs"]
    )

    # Attempt to calculate the moment magnitude from the moment tensor
    # components. This might fail because the values don't make sense or are
    # not provided
    try:
        moment = seismic_moment(mt=[float(source_dict["Mxx"]),
                                    float(source_dict["Mzz"]),
                                    float(source_dict["Mxz"])
                                    ])
        moment_mag = moment_magnitude(moment=moment)
    except Exception as e:
        moment_mag = None

    magnitude = Magnitude(
        resource_id=_get_resource_id(source_name, "magnitude"),
        mag=moment_mag, magnitude_type="Mw", origin_id=origin.resource_id.id
    )

    event = Event(resource_id=_get_resource_id(name=source_name,
                                               res_type="event"))
    event.origins.append(origin)
    event.magnitudes.append(magnitude)

    event.preferred_origin_id = origin.resource_id.id
    event.preferred_magnitude_id = magnitude.resource_id.id

    return event


def read_forcesolution(path_to_source, default_time="2000-01-01T00:00:00"):
    """
    Create a barebones Pyatoa Source object from a FORCESOLUTION Source file,
    which mimics the behavior of the more complex ObsPy Event object
    """
    with open(path_to_source, "r") as f:
        lines = f.readlines()

    # Place values from file into dict
    source_dict = {}
    for line in lines[:]:
        if ":" not in line:
            continue
        key, val = line.strip().split(":")
        source_dict[key] = val

    # First line has the source name
    source_dict["source_id"] = lines[0].strip().split()[-1]

    event = Source(
        resource_id=f"pysep:source/{source_dict['source_id']}",
        origin_time=default_time, latitude=source_dict["latorUTM"],
        longitude=source_dict["longorUTM"], depth=source_dict["depth"]
    )

    return event


def read_stations(path_to_stations):
    """
    Convert a SPECFEM STATIONS file into an ObsPy Inventory object.

    Specfem3D STATION files contain no channel or location information, so
    the inventory can only go down to the station level.

    .. note::
        This assumes a row structure for the station file is
        STA, NET, LAT [deg], LON [deg], ELEVATION [m], BURIAL [m]

    :type path_to_stations: str
    :param path_to_stations: the path to the STATIONS file that is associated
        with the Specfem3D DATA directory
    :rtype: obspy.core.inventory.Inventory
    :return: a station-level Inventory object
    :raises ValueError: if latitude and longitude values are not in geographic
        coordinates (i.e., in cartesian coordinates). Thrown by the init of the
        Station class.
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

    return Inventory(networks=list_of_networks, source="PySEP")


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


def write_sem(st, unit, path="./", time_offset=0):
    """
    Write an ObsPy Stream in the two-column ASCII format that Specfem uses

    :type st: obspy.core.stream.Stream
    :param st: stream containing synthetic data to be written
    :type unit: str
    :param unit: the units of the synthetic data, used for file extension, must
        be 'd', 'v', 'a' for displacement, velocity, acceleration, resp.
    :type path: str
    :param path: path to save data to, defaults to cwd
    :type time_offset: float
    :param time_offset: optional argument to offset the time array. Sign matters
        e.g. time_offset=-20 means t0=-20
    """
    assert(unit.lower() in ["d", "v", "a"]), "'unit' must be 'd', 'v' or 'a'"
    for tr in st:
        s = tr.stats
        fid = f"{s.network}.{s.station}.{channel_code(s.delta)}X{s.channel[-1]}"
        fid = os.path.join(path, f"{fid}.sem{unit.lower()}")
        data = np.vstack((tr.times() + time_offset, tr.data)).T
        np.savetxt(fid, data, ["%13.7f", "%17.7f"])


def write_pysep_station_file(inv, event, fid="./station_list.txt", 
                             order_station_list_by=None):
    """
    Write a list of station codes, distances, etc. useful for understanding
    characterization of all collected stations

    :type event: obspy.core.event.Event
    :param event: optional user-provided event object which will force a
        skip over QuakeML/event searching
    :type inv: obspy.core.inventory.Inventory
    :param inv: optional user-provided inventory object which will force a
        skip over StationXML/inventory searching
    :type fid: str
    :param fid: name of the file to write to. defaults to ./station_list.txt
    :type order_station_list_by: str
    :param order_station_list_by: how to order the stations written to file.
        Available are: network, station, latitude, longitude, elevation, burial.
        If not given, order is determined by the internal sorting of the
        input Inventory
    """
    # Key indices correspond to stations list
    keys = ["station", "network", "latitude", "longitude", "distance",
            "azimuth"]
    if order_station_list_by and order_station_list_by not in keys:
        logger.warning(f"`order_station_by` must be in {keys}, "
                       f"setting default")
        order_station_by = None

    event_latitude = event.preferred_origin().latitude
    event_longitude = event.preferred_origin().longitude
    
    stations = [] 
    for net in inv:
        for sta in net:
            dist_m, az, baz = gps2dist_azimuth(lat1=event_latitude,
                                               lon1=event_longitude,
                                               lat2=sta.latitude,
                                               lon2=sta.longitude
                                               )
            dist_km = dist_m * 1E-3
            stations.append([sta.code, net.code, sta.latitude, sta.longitude,
                             dist_km, az])

    # Set the order of the station file based on the order of keys
    if order_station_list_by:
        idx = keys.index(order_station_list_by)
        stations.sort(key=lambda x: x[idx])

    with open(fid, "w") as f:
        for s in stations:
            # Order follows that listed in 'keys'
            f.write(f"{s[0]:<6} {s[1]:<2} {s[2]:9.4f} {s[3]:9.4f} {s[4]:8.3f} "
                    f"{s[5]:6.2f}\n")


def write_stations_file(inv, fid="./STATIONS", order_by=None,
                        use_elevation=False, burials=None):
    """
    Write a SPECFEM STATIONS file given an ObsPy inventory object

    .. note::
        If topography is implemented in your mesh, elevation values should be
        set to 0 which means 'surface' in SPECFEM.

    .. note::
        burial information is not contained in the ObsPy inventory so it is
        always set to a constant value, which can be given as input. default 0

    :type inv: obspy.core.inventory.Inventory
    :param inv: Inventory object with station locations to write
    :type fid: str
    :param fid: path and file id to save the STATIONS file to.
    :type order_by: str
    :param order_by: how to order the stations written to file. Available
        are: network, station, latitude, longitude, elevation, burial.
        If not given, order is determined by the internal sorting of the
        input Inventory
    :type use_elevation: bool
    :param use_elevation: if True, sets the actual elevation value from the
        inventory. if False, sets elevation to 0.
    :type burials: list of float
    :param burials: If the User has burial information they want to be used as
        the last column. Length of `burials` must match the number of stations
        in the inventory when traversing by network and station
    """
    # Simply generate lists by traversing through the inventory
    i = 0
    stations = []
    keys = ["network", "station", "latitude",
            "longitude", "elevation", "burial"]
    if order_by:
        assert(order_by in keys), f"`order_by` must be in {keys}"

    with open(fid, "w") as f:
        for net in inv:
            for sta in net:
                if use_elevation:
                    elevation = sta.elevation
                else:
                    elevation = 0.
                if burials:
                    burial = burials[i]
                else:
                    burial = 0.
                stations.append([sta.code, net.code, sta.latitude,
                                 sta.longitude, elevation, burial])
                i += 1

    # Set the order of the station file based on the order of keys
    if order_by:
        idx = keys.index(order_by)
        stations.sort(key=lambda x: x[idx])

    with open(fid, "w") as f:
        for s in stations:
            f.write(f"{s[0]:>6}{s[1]:>6}{s[2]:12.4f}{s[3]:12.4f}"
                    f"{s[4]:7.1f}{s[5]:7.1f}\n"
                    )


def write_cat_to_event_list(cat, fid_out="event_input.txt"):
    """
    Writes out an event Catalog (ObsPy object) information to an ASCII file that
    PySEP can use for collecting data. The format of the output file is

    ORIGIN_TIME LONGITUDE LATITUDE DEPTH[KM] MAGNITUDE

    :type cat: obspy.core.catalog.Catalog
    :param cat: Catalog of events to write out
    """
    with open(fid_out, "w") as f:
        for event in cat:
            origintime = str(event.preferred_origin().time)
            longitude = event.preferred_origin().longitude
            latitude = event.preferred_origin().latitude
            depth = event.preferred_origin().depth * 1E-3
            mag = event.preferred_magnitude().mag

            f.write(f"{origintime:<31}{longitude:8.2f}{latitude:8.2f}"
                    f"{depth:8.2f}{mag:6.1f}\n")
    

