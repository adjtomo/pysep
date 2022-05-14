"""
Utils to honor file formats from SAC and CAP
"""
from obspy.core.stream import Stream
from obspy.geodetics import gps2dist_azimuth, kilometer2degrees

from pysep import logger


def append_sac_headers(st, event, inv):
    """
    Wrapper for trace header appending to get a loop and some logic in

    :type st: obspy.core.stream.Stream
    :param st: Stream to append SAC header to
    :type event: obspy.core.event.Event
    :param event: Event with metadata for SAC header
    :type inv: obspy.core.inventory.Inventory
    :param event: StationXML with metadata for SAC header
    :rtype: obspy.core.stream.Stream
    :return: Stream with SAC headers, those that could not be appended to
        have been removed from the stream
    """
    st_out = Stream()
    for tr in st[:]:
        try:
            st_out.append(_append_sac_headers_trace(tr, event, inv))
        except Exception as e:
            logger.warning(f"{tr.get_id()} cannot write SAC headers, removing")
    return st_out


def _append_sac_headers_trace(tr, event, inv):
    """
    Append SAC headers to ObsPy streams given event and station metadata

    Rewritten from: `util_write_cap.add_sac_metadata()`

    TODO Add back in information removed from original function
        * Add P arrival time to stats.sac["a"]
        * Add TauP phase arrivals, writing to t5, kt5, user1, kuser1
        * Add sensor type somewhere, previously stored in KT? (used for pics)

    :type tr: obspy.core.trace.Trace
    :param tr: Trace to append SAC header to
    :type event: obspy.core.event.Event
    :param event: Event with metadata for SAC header
    :type inv: obspy.core.inventory.Inventory
    :param event: StationXML with metadata for SAC header
    :rtype: obspy.core.trace.Trace
    :return: Trace with appended SAC header
    """
    net_code, sta_code, loc_code, cha_code = tr.get_id().split(".")
    # An inventory that is narrowed down to a per-channel basis
    inv_unique = inv.select(network=net_code, station=sta_code,
                            channel=cha_code)
    # All of these are subsets of an inventory
    net = inv_unique[0]
    sta = net[0]
    cha = sta[0]

    dist_m, az, baz = gps2dist_azimuth(
        lat1=event.preferred_origin().latitude,
        lon1=event.preferred_origin().longitude,
        lat2=sta.latitude, lon2=sta.longitude
    )
    dist_km = dist_m * 1E-3  # units: m -> km
    dist_deg = kilometer2degrees(dist_km)  # spherical earth approximation

    sac_header = {
        "evla": event.preferred_origin().latitude,
        "evlo": event.preferred_origin().longitude,
        "evdp": event.preferred_origin().depth / 1E3,  # depth in km
        "mag": event.preferred_magnitude(),
        "stla": sta.latitude,
        "stlo": sta.longitude,
        "stel": sta.elevation / 1E3,  # elevation in km
        # 'kevnm' renders to, e.g.,: 2012-04-04T142142
        "kevnm": event.preferred_origin().time.strftime("%Y-%m-%dT%H%M%S"),
        "dist": dist_km,
        "az": az,  # degrees
        "baz": baz,  # degrees
        "gcarc": dist_deg,  # degrees
        "cmpinc": cha.dip,  # channel dip/inclination in degrees
        "cmpaz": cha.azimuth,  # channel azimuth in degrees
        "lpspol": 0,  # 1 if left-hand polarity (usually no in passive seis)
        "lcalda": 1,  # 1 if DIST, AZ, BAZ, GCARC to be calc'd from metadata
    }

    tr.stats.sac = sac_header
    tr.stats.back_azimuth = baz

    return tr
