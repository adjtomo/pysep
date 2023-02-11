"""
Pysep-specific formatting functions
"""
from obspy import Catalog, UTCDateTime
from obspy.clients.iris import Client


def channel_code(dt):
    """
    Specfem outputs seismograms with channel band codes set by IRIS. Instrument
    codes are always X for synthetics, but band code will vary with the sampling
    rate of the data, return the correct code given a sampling rate.
    Taken from Appenix B of the Specfem3D cartesian manual (June 15, 2018)

    :type dt: float
    :param dt: sampling rate of the data in seconds
    :rtype: str
    :return: band code as specified by Iris
    :raises KeyError: when dt is specified incorrectly
    """
    if dt >= 1:
        return "L"  # long period
    elif 0.1 < dt < 1:
        return "M"  # mid period
    elif 0.0125 < dt <= 0.1:
        return "B"  # broad band
    elif 0.001 <= dt <= 0.0125:
        return "H"  # high broad band
    elif 0.004 <= dt < 0.001:
        return "C"
    elif 0.001 <= dt < 0.0002:
        return "F"
    else:
        raise KeyError("Channel code does not exist for this value of 'dt'")


def get_codes(st=None, choice=None, suffix=None, up_to=True):
    """
    Get station codes from the internal stream attribute, where station
    codes are formatted NN.SSS.LL.CCc where N=network, S=station,
    L=location, C=channel, and c=component

    :type st: obspy.core.stream.Stream
    :param st: Stream to get codes from by running: tr.get_id()
    :type choice: str
    :param choice: choice of the part of the code returned, available:
        * 'network': return unique network codes (e.g., NN)
        * 'station': return unique network + station codes (e.g., NN.SSS)
        * 'location': return up to location (e.g., NN.SSS.LL)
        * 'channel': return up to channel, no component (e.g., NN.SSS.LL.CC)
        * else: return full station code (e.g., NN.SSS.LL.CCc)
    :type suffix: str
    :param suffix: append a string `suffix` to the end of each code. Used to
        place wildcards at the end of a code, e.g., `suffix`=='?' with
        `choice`=='channel' will give codes like NN.SSS.LL.CC?
    :type up_to: bool
    :param up_to: if True gives code up to the `choice`, e.g.,
        `choice`=='station' would return NN.SSS. If False, returns only `choice`
        code, so e.g., SSS for `choice`=='station'
    :rtype: list
    :return: unique station codes filtered by choice
    """
    full_codes = [tr.get_id() for tr in st]
    if choice == "network":
        codes = [code.split(".")[0] for code in full_codes]
    elif choice == "station":
        if up_to:
            codes = [".".join(code.split(".")[:2]) for code in full_codes]
        else:
            codes = [code.split(".")[1] for code in full_codes]
    elif choice == "location":
        if up_to:
            codes = [".".join(code.split(".")[:3]) for code in full_codes]
        else:
            codes = [code.split(".")[2] for code in full_codes]
    elif choice == "channel":
        if up_to:
            codes = [code[:-1] for code in full_codes]
        else:
            codes = [code[-1] for code in full_codes]
    else:
        codes = full_codes

    if suffix is not None:
        codes = [f"{code}{suffix}" for code in codes]

    return sorted(list(set(codes)))


def get_data_availability(cat, inv):
    """
    Determine data availability based on whether stations are 'on' for a
    given event origin time. Does not check waveforms, only station
    metadata, so not foolproof.

    :type cat: obspy.core.catalog.Catalog
    :param cat: Catalog of events to consider. Events must include origin
        information `latitude` and `longitude`
    :type inv: obspy.core.inventory.Inventory
    :param inv: Inventory of stations to consider
    :rtype: dict
    :return: keys are event resource ids and values are IDs for stations
        that were on during the event origin time
    """
    # Collect install and removal (if applicaple) for all stations
    station_times = {}
    for net in inv:
        for sta in net:
            if sta.end_date is None:
                sta.end_date = UTCDateTime()  # set to right now
            station_times[f"{net.code}.{sta.code}"] = (sta.start_date,
                                                       sta.end_date)

    # Check event origin time against station uptime
    data_avail = {}
    for event in cat:
        data_avail[event.resource_id.id] = []
        for sta, time in station_times.items():
            start_date, end_date = time
            # Check that event origin time falls between start and end date
            if start_date <= event.preferred_origin().time <= end_date:
                data_avail[event.resource_id.id].append(sta)

    return data_avail


def format_event_tag(event):
    """
    Generate a unique event tag based on event origin time and location using
    Flinn Engdahl regions

    :type event: obspy.core.event.Event
    :param event: event to generate tag from
    :rtype: str
    :return: event_name + FE region
    """
    event_lat = event.preferred_origin().latitude
    event_lon = event.preferred_origin().longitude

    client = Client()
    region = client.flinnengdahl(lat=event_lat, lon=event_lon, rtype="region")
    # e.g. NORTH ISLAND, NEW ZEALAND > NORTH_ISLAND_NEW_ZEALAND
    region = region.replace(" ", "_").replace(",", "").upper()

    # Presumed to be event name specified by origintime, e..g, 20090407T201255
    event_time = event.preferred_origin().time.strftime("%Y-%m-%dT%H%M%S")

    return f"{event_time}_{region}"


def format_event_tag_legacy(event):
    """
    Generate a unique event tag based on the event time. This was how the
    previous version of PySEP named directories and files. Replaces the old
    `otime2eid` from `util_helpers`

    :type event: obspy.core.event.Event
    :param event: event to generate tag from
    :rtype: str
    :return: event_name specified by event time
    """
    return event.preferred_origin().time.strftime("%Y%m%d%H%M%S%f")[:-3]


def index_cat(cat, idxs):
    """
    ObsPy Catalog does not allow indexing by a list of values
    (e.g., cat[0, 1, 3]) so this convenience function takes care of that by
    forming a new catalog of events chosen by indices

    :type idxs: list of int
    :param idxs: list of indices to index catalog by
    :type cat: obspy.core.catalog.Catalog
    :param cat: Catalog to index. If not given defaults to internal Cat
    """
    cat_out = Catalog()
    for idx in idxs:
        cat_out.append(cat[idx])
    return cat_out