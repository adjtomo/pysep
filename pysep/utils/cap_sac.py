"""
Utils to honor file formats from SAC and CAP

For SAC header names and descriptions, see:
    http://www.adc1.iris.edu/files/sac-manual/manual/file_format.html
"""
import os
import numpy as np
from obspy.core.stream import Stream
from obspy.geodetics import gps2dist_azimuth, kilometer2degrees

from pysep import logger
from pysep.utils.fmt import format_event_tag


def write_cap_weights_files(st, event, path_out="./", order_by="dist"):
    """
    Write CAP (Cut-and-Paste) moment tensor inversion code weight files,
    assuming that SAC headers are already present.

    Replaces `write_cap_weights`

    The weight file has columns corresponding to the following:

        0: EVENT_STATION_ID
        1: DIST_KM
        2: BODY_Z
        3: BODY_R
        4: SURF_Z
        5: SURF_R
        6: SURF_T
        7: P_ARRIVAL
        8: LEGACY (unused)
        9: S_ARRIVAL
        10: LEGACY (unused)
        11: STATIC CORRECTION RAYLEIGH

    TODO re-add Ptime setting with event.picks
    :return:
    """
    assert(order_by in ["dist", "az", "code"]), (f"CAP weights sorting must be "
                                                 f"by 'dist', 'az', or 'code'")

    # Standard look to the white space in the weights' files
    weight_fmt = ("{code:>32}{dist:8.2f} {body_z:>3}{body_r:>2}   "
                  "{surf_z:>2}{surf_r:>2}{surf_t:>2}   "
                  "{p_arr:6.2f}{leg:>2}   "
                  "{s_arr:6.2f}{leg:>2}   "
                  "{corr:>2}\n")

    # Define pre-set keys for differently weighted weights files
    weight_files = {"weights.dat": {"body_z": 1, "body_r": 1,
                                    "surf_z": 1, "surf_r": 1, "surf_t": 1},
                    "weights_body.dat": {"body_z": 1, "body_r": 1,
                                         "surf_z": 0, "surf_r": 0, "surf_t": 0},
                    "weights_surf.dat": {"body_z": 0, "body_r": 0,
                                         "surf_z": 1, "surf_r": 1, "surf_t": 1},
                    }

    # e.g. NN.SSS.LL.CC?; wildcard component
    codes = list(set([f"{tr.get_id()[:-1]}?" for tr in st]))
    code_list = []
    p_arrival = 0.  # default value
    for code in codes:
        net, sta, loc, cha = code.split(".")
        for pick in event.picks:
            # TODO place p_arrival calculation here
            pass
        # Grab any component as they all have the same dist and azimuth
        tr = st.select(network=net, station=sta, location=loc, channel=cha)[0]
        # [code, dist_km, az, p_arrival]
        code_list.append([f"{tr.stats.sac['kevnm']}.{code}",
                          tr.stats.sac["dist"], tr.stats.sac["az"],
                          p_arrival])

    # Order codes based on distance, name or azimuth
    idx = ["code", "dist", "az"].index(order_by)
    code_list = np.array(code_list)
    ordered_codes = code_list[code_list[:, idx].argsort()]

    logger.info("writing CAP weight files")
    for basename, weights in weight_files.items():
        fid = os.path.join(path_out, basename)
        logger.debug(f"writing CAP weight file for {len(ordered_codes)} "
                     f"station(s), ordered by '{order_by}': '{fid}'")
        with open(fid, "w") as f:
            for vals in ordered_codes:
                code, dist, az, p_arr = vals
                # note: code drops the wild card that we appended earlier
                f.write(weight_fmt.format(
                    code=code[:-1], dist=float(dist), body_z=weights["body_z"],
                    body_r=weights["body_r"], surf_z=weights["surf_z"],
                    surf_r=weights["surf_r"], surf_t=weights["surf_t"],
                    p_arr=float(p_arr), leg=0, s_arr=0., corr=0)
                )


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
    for tr in st.copy()[:]:
        try:
            st_out.append(_append_sac_headers_trace(tr, event, inv))
        except Exception as e:
            logger.warning(f"{tr.get_id()} can't write SAC headers: {e}")
    return st_out


def _append_sac_headers_trace(tr, event, inv):
    """
    Append SAC headers to ObsPy streams given event and station metadata

    Rewritten from: `util_write_cap.add_sac_metadata()`

    TODO Add back in information removed from original function
        * Add P arrival time to stats.sac["a"]
        * Add TauP phase arrivals, writing to t5, kt5, user1, kuser1
        * Add sensor type somewhere, previously stored in KT? (used for picks)

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

    # Get some information from TauP regarding incident angle, takeoff angle

    sac_header = {
        "evla": event.preferred_origin().latitude,
        "evlo": event.preferred_origin().longitude,
        "evdp": event.preferred_origin().depth / 1E3,  # depth in km
        "mag": event.preferred_magnitude().mag,
        "stla": sta.latitude,
        "stlo": sta.longitude,
        "stel": sta.elevation / 1E3,  # elevation in km
        "kevnm": format_event_tag(event),
        "dist": dist_km,
        "az": az,  # degrees
        "baz": baz,  # degrees
        "gcarc": dist_deg,  # degrees
        "cmpinc": cha.dip,  # channel dip/inclination in degrees  TODO this is incident angle
        "cmpaz": cha.azimuth,  # channel azimuth in degrees
        "lpspol": 0,  # 1 if left-hand polarity (usually no in passive seis)
        "lcalda": 1,  # 1 if DIST, AZ, BAZ, GCARC to be calc'd from metadata
    }

    tr.stats.sac = sac_header
    tr.stats.back_azimuth = baz

    return tr


def format_sac_headers_post_rotation(st):
    """
    SAC headers do not update when rotating so we need to apply manual
    changes to the azimuth, inclination and naming values

    TODO is this necessary? Who is using the SAC headers and what info
        do they need? Or can we just re-run SAC header appending?
    :param st:
    :return:
    """
    azimuth_dict = {"E": 90., "N": 0., "Z": 0.}
    inclin_dict = {"E": 0., "N": 0., "Z": -90.}

    st_out = st.copy()
    for tr in st_out:
        comp = tr.stats.component
        tr.stats.sac["kcmpnm"] = tr.stats.channel  # TODO component or channel?
        if comp in azimuth_dict.keys():
            tr.stats.sac["cmpaz"] = azimuth_dict[comp]
            tr.stats.sac["cmpinc"] = inclin_dict[comp]
        elif comp == "R":
            tr.stats.sac["cmpaz"] = tr.stats.sac["az"]
        elif comp == "T":
            tr.stats.sac["cmpaz"] = (tr.stats.sac["az"] + 90) % 360

    return st_out

