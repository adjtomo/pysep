"""
Utils to honor file formats from SAC and CAP

For SAC header names and descriptions, see:
    http://www.adc1.iris.edu/files/sac-manual/manual/file_format.html
"""
import os
import numpy as np
from obspy import UTCDateTime
from obspy.core.stream import Stream
from obspy.geodetics import gps2dist_azimuth, kilometer2degrees

from pysep import logger
from pysep.utils.fmt import format_event_tag_legacy
from pysep.utils.fetch import get_taup_arrivals

# SAC HEADER CONSTANTS DEFINING NON-INTUITIVE QUANTITIES
SACDICT = {
    "first_arrival_time": "a",
    "p_arrival_time": "t5",
    "p_incident_angle": "user1",
    "p_takeoff_angle": "user3",
    "s_arrival_time": "t6",
    "s_incident_angle": "user2",
    "s_takeoff_angle": "user4",
}


def write_cap_weights_files(st, path_out="./", order_by="dist"):
    """
    Write CAP (Cut-and-Paste) moment tensor inversion code weight files,
    assuming that SAC headers are already present.

    TODO re-add Ptime setting with event.picks

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

    :type st: obspy.core.stream.Stream
    :param st: input stream to use to write CAP weight files, expected to
        have SAC header
    :type path_out: str
    :param path_out: path to write the weight file, filenames are set
        by default inside the function
    :type order_by: str
    :param order_by: how to order the list of stations that gets written out
        available options are:
        * dist: order by smallest to largest source-receiver distance (default)
        * az: order by smallest to largest azimuth (deg)
        * code: order alphabetically by station name
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
    Append SAC headers to ObsPy streams given event and station metadata.
    Also add 'back_azimuth' to Stream stats which can be used for rotation.

    Rewritten from: `util_write_cap.add_sac_metadata()`

    TODO Add back in information removed from original function
        * Add sensor type somewhere, previously stored in KT? (used for picks)

    .. note::
        We explicitely set 'iztype, 'b' and 'e' in the SAC header to tell ObsPy
        that the trace start is NOT the origin time. Otherwise all the relative
        timing (e.g., picks) will be wrong.

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
    dist_m, az, baz = gps2dist_azimuth(
        lat1=event.preferred_origin().latitude,
        lon1=event.preferred_origin().longitude,
        lat2=sta.latitude, lon2=sta.longitude
    )
    dist_km = dist_m * 1E-3  # units: m -> km
    dist_deg = kilometer2degrees(dist_km)  # spherical earth approximation

    otime = event.preferred_origin().time

    sac_header = {
        "iztype": 9,  # Ref time equivalence, IB (9): Begin time
        "b": tr.stats.starttime - event.preferred_origin().time,  # begin time
        "e": tr.stats.npts * tr.stats.delta,  # end time
        "evla": event.preferred_origin().latitude,
        "evlo": event.preferred_origin().longitude,
        "stla": sta.latitude,
        "stlo": sta.longitude,
        "stel": sta.elevation / 1E3,  # elevation in km
        "kevnm": format_event_tag_legacy(event),  # only take date code
        "nzyear": otime.year,
        "nzjday": otime.julday,
        "nzhour": otime.hour,
        "nzmin": otime.minute,
        "nzsec": otime.second,
        "nzmsec": otime.microsecond,
        "dist": dist_km,
        "az": az,  # degrees
        "baz": baz,  # degrees
        "gcarc": dist_deg,  # degrees
        "lpspol": 0,  # 1 if left-hand polarity (usually no in passive seis)
        "lcalda": 1,  # 1 if DIST, AZ, BAZ, GCARC to be calc'd from metadata
    }
    # Some Inventory objects will not go all the way to channel, only to station
    try:
        cha = sta[0]
        sac_header["cmpinc"] = cha.dip  # channel dip/inclination in degrees
        sac_header["cmpaz"] = cha.azimuth  # channel azimuth in degrees
    except IndexError:
        pass

    # Allow for no magnitude and no depth information
    try:
        evdp = event.preferred_origin().depth / 1E3  # units: km
        sac_header["evdp"] = evdp
    except Exception:  # NOQA
        pass

    try:
        mag = event.preferred_magnitude().mag
        sac_header["mag"] = mag
    except Exception:  # NOQA
        pass

    # Append SAC header and include back azimuth for rotation
    tr.stats.sac = sac_header
    tr.stats.back_azimuth = baz

    return tr


def format_sac_header_w_taup_traveltimes(st, model="ak135",
                                         phase_list=("ttall",)):
    """
    Add TauP travel times to the SAC headers using information in the SAC header
    Also get some information from TauP regarding incident angle, takeoff angle
    Hardcoded to only look at P and S arrivals (both upgoing and downgoing)

    TODO Probably find better ways to store arrival time and incident angles

    .. note::
        This function expects that the Stream has been formatted with SAC header

    .. note::
        SAC header writing could probably be in a loop, but I think it's more
        readable to see P and S values getting written separately.

    :type st: obspy.core.stream.Stream
    :param st: Stream object with SAC headers which will be written to with
        new SAC header attributser
    :type model: str
    :param model: name of the TauP model to use for arrival times etc.
        defaults to 'ak135'
    :type phase_list: list of str
    :param phase_list: phase names to get ray information from TauP with.
        Defaults to 'ttall', which is ObsPy's default for getting all phase
        arrivals. Must match Phases expected by TauP (see ObsPy TauP
        documentation for acceptable phases).
    """
    st_out = st.copy()
    phase_dict = get_taup_arrivals(st=st, model=model, phase_list=phase_list)

    # Arrivals may return multiple entires for each phase, pick earliest
    for tr in st_out[:]:
        net_sta = ".".join(tr.get_id().split(".")[:2])  # NN.SSS

        # Missing SAC header values may cause certain or all stations to not 
        # be present in the `phase_dict`
        if net_sta not in phase_dict:
            logger.warning(f"{tr.get_id()} not present in TauP arrivals, cant "
                           f"append arrival time SAC headers")
            continue
        arrivals = phase_dict[net_sta]
        # If the TauP arrival calculation fails, `arrivals` will be empty
        if not arrivals:
            logger.warning(f"{tr.get_id()} has no TauP phase arrivals, cant "
                           f"append arrival time SAC headers")
            continue

        # Find earliest arriving phase (likely same as P phase but maybe not)
        idx_times = [(i, a.time) for i, a in enumerate(arrivals)]
        idx, _ = min(idx_times, key=lambda x: x[1])  # find index of min time

        phase = arrivals[idx]  # Earliest Arrival object
        tr.stats.sac[SACDICT["first_arrival_time"]] = phase.time  # 'a'
        tr.stats.sac[f"k{SACDICT['first_arrival_time']}"] = f"{phase.name}" #ka

        # Find earliest arriving P phase (must start with letter 'P')
        idx_times = [(i, a.time) for i, a in enumerate(arrivals) if
                     a.name.upper().startswith("P")]
        idx, _ = min(idx_times, key=lambda x: x[1])  # find index of min time

        p = arrivals[idx]  # Earliest P-wave Arrival object
        tr.stats.sac[SACDICT["p_arrival_time"]] = p.time
        tr.stats.sac[f"k{SACDICT['p_arrival_time']}"] = f"{p.name}"

        # P phase incident angle (ia) and takeoff angle (ta)
        tr.stats.sac[SACDICT["p_incident_angle"]] = p.incident_angle
        tr.stats.sac[f"k{SACDICT['p_incident_angle']}"] = f"p_IncAng"
        tr.stats.sac[SACDICT["p_takeoff_angle"]] = arrivals[idx].takeoff_angle
        tr.stats.sac[f"k{SACDICT['p_takeoff_angle']}"] = f"p_TkfAng"

        # Find earliest arriving S phase (must start with letter 'S')
        idx_times = [(i, a.time) for i, a in enumerate(arrivals) if
                     a.name.upper().startswith("S")]
        idx, _ = min(idx_times, key=lambda x: x[1])
        s = arrivals[idx]  # Earliest S-wave Arrival object

        tr.stats.sac[SACDICT["s_arrival_time"]] = s.time
        tr.stats.sac[f"k{SACDICT['s_arrival_time']}"] = f"{s.name}"

        tr.stats.sac[SACDICT["s_incident_angle"]] = s.incident_angle
        tr.stats.sac[f"k{SACDICT['s_incident_angle']}"] = f"s_IncAng"
        tr.stats.sac[SACDICT["s_takeoff_angle"]] = s.takeoff_angle
        tr.stats.sac[f"k{SACDICT['s_takeoff_angle']}"] = f"s_TkfAng"

    return st_out


def format_sac_headers_post_rotation(st):
    """
    SAC headers do not update when rotating so we need to apply manual
    changes to the azimuth, inclination and naming values

    TODO is this necessary? Who is using the SAC headers and what info
        do they need? Or can we just re-run SAC header appending?

    :type st: obspy.core.stream.Stream
    :param st: Stream to append SAC headers for
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


def origin_time_from_sac_header(sac_header):
    """
    Build a UTCDateTime origin time from values in the SAC header appended to
    an ObsPy trace.

    :type sac_header: obspy.core.util.attribdict.AttribDict
    :param sac_header: SAC header built by `append_sac_header()`
    :rtype: UTCDateTime
    :return: event origin time built from SAC header
    """
    year = sac_header["nzyear"]
    jday = sac_header["nzjday"]
    hour = sac_header["nzhour"]
    min_ = sac_header["nzmin"]
    sec_ = sac_header["nzsec"]
    msec = sac_header["nzmsec"]
    time_string = f"{year}-{jday:0>3}T{hour:0>2}:{min_:0>2}:{sec_:0>2}.{msec}"

    return UTCDateTime(time_string)
