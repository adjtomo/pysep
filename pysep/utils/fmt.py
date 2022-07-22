"""
Pysep-specific formatting functions
"""
from obspy.clients.iris import Client


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
    region = region.replace(" ", "_").replace(",", "").upper()

    # Presumed to be event name specified by origintime
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
