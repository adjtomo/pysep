"""
Utilities for processing waveforms
"""
import numpy as np
from scipy import signal

from obspy import Stream
from pysep import logger


def get_codes(st=None, choice=None, suffix=None):
    """
    Get station codes from the internal stream attribute, where station
    codes are formatted NN.SSS.LL.CCc where N=network, S=station,
    L=location, C=channel, and c=component

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
    :rtype: list
    :return: unique station codes filtered by choice
    """
    full_codes = [tr.get_id() for tr in st]
    if choice == "network":
        codes = [code.split(".")[0] for code in full_codes]
    elif choice == "station":
        codes = [".".join(code.split(".")[:1]) for code in full_codes]
    elif choice == "location":
        codes = [".".join(code.split(".")[:2]) for code in full_codes]
    elif choice == "channel":
        codes = [code[:-1] for code in full_codes]
    else:
        codes = full_codes

    if suffix is not None:
        codes = [f"{code}{suffix}" for code in codes]

    return list(set(codes))


def format_streams_for_rotation(st):
    """
    ObsPy requires specific channel naming to get stream rotation working.
    Comb through the stream and check that this is correct before passing
    stream to rotation
    :param st:
    :return:
    """
    logger.info(f"formatting stream of length {len(st)} for component rotation")
    codes = get_codes(st, choice="channel", suffix="?")

    st_out = Stream()
    for code in codes:
        net, sta, loc, cha = code.split(".")
        st = st.select(network=net, station=sta, location=loc,
                       channel=cha).copy()

        missing_vertical, missing_horizontal = _check_component_availability(st)

        if missing_vertical and missing_horizontal:
            logger.info(f"{code} is missing vertical and at least 1 horizontal,"
                        f"skipping over")
            continue
        elif missing_vertical and not missing_horizontal:
            logger.info(f"{code} is missing vertical, replacing vertical "
                        f"component with zeros")
            st = _append_null_trace(st, component="Z", cmpaz=0, cmpinc=-90.)
        elif not missing_vertical and missing_horizontal:
            logger.info(f"{code} is missing at least 1 horizontal component, "
                        f"removing existing horizontal data and writing zeros "
                        f"to North and East components")
            for comp in ["E", "N", "1", "2"]:
                try:
                    st.remove(st.select(component=comp)[0])
                    logger.debug(f"{code} is removing existing comp: {comp}")
                except IndexError:
                    continue
            st = _append_null_trace(st, component="E", cmpaz=90., cpmpinc=0.)
            st = _append_null_trace(st, component="N", cmpaz=0., cpmpinc=0.)

        st_out += st

    logger.info(f"stream formatting returned {len(st_out)} traces")


    return st_out


def _append_null_trace(st, component, cmpaz=None, cmpinc=None):
    """
    Create a copy of a trace with zeroed out data, used for rotation
    """
    st_out = st.copy()
    trace_null = st_out[0].copy()
    trace_null.data *= 0
    trace_null.stats.component = component
    if cmpaz is not None:
        trace_null.stats.sac["cmpaz"] = cmpaz
    if cmpinc is not None:
        trace_null.stats.sac["cmpinc"] = cmpinc
    st_out.append(trace_null)

    return st_out


def _check_component_availability(st):
    """
    Quick logic check to see if stream is missing horizontals or vertical
    """
    comps = sorted([tr.stats.component.upper() for tr in st])
    missing_vertical, missing_horizontal = False, False

    if comps:
        if "Z" not in comps:
            missing_vertical = True

        # Any one of these criteria will trigger a 'missing_horizontal' flag
        for pairs in [("N", "E"), ("R", "T"), ("1", "2")]:
            first, second = pairs
            if first in comps and second not in comps:
                missing_horizontal = True
            elif second in comps and first not in comps:
                missing_horizontal = True
    # Deal with the case where we have an empty Stream
    else:
        missing_vertical, missing_horizontal = True, True

    return missing_vertical, missing_horizontal


def rotate_to_UVW(st):
    """
    UVW orthogonal frame rotation

    In Symmetric Triaxial Seismometers, the sensing elements are also arranged
    to be mutually orthogonal, but instead of one axis being vertical, all
    three are inclined upwards from the horizontal at precisely the same angle,
    as if they were aligned with the edges of a cube balanced on a corner.

    TODO test this, not sure how it's supposed to work

    See reference:
    http://link.springer.com/referenceworkentry/10.1007/978-3-642-36197-5_194-1

    Rotation matrix reference:
    http://link.springer.com/referenceworkentry/10.1007/
                                                  978-3-642-36197-5_194-1#page-1

    :type st: obspy.core.stream.Stream
    :param st: three-component stream
    :return:
    """
    st_out = st.copy()
    assert(len(st) == 3), f"UVW rotation requires 3 component stream input"

    xyz = np.array([st_out[0].data], [st_out[1].data], [st_out[2].data])
    xyz = np.asmatrix(xyz)

    rotation_matrix =  np.matrix([[2, 0, np.sqrt(2)],
                                  [-1, np.sqrt(3), np.sqrt(2)],
                                  [-1, -np.sqrt(3), np.sqrt(2)]]
                                 )
    rotation_matrix *= 1 / np.sqrt(6)

    # Theta is the angle between the axis and the horizontal plane
    theta = (np.arcsin(1 / np.sqrt(3))) * (180.0 / np.pi)

    uvw = rotation_matrix * xyz
    for i, comp in enumerate(["U", "V", "W"]):
        st_out[i].data = np.reshape(np.asarray(uvw[i], uvw.shape[1]))
        st_out[i].stats.component = comp
        st_out[i].stats.sac["kcmpnm"] = st_out[i].stats.channel
        st_out[i].stats.sac["cmpaz"] *= (120 * i)  # azimuth
        st_out[i].stats.sac["cmpinc"] -= theta  # inclination
        if comp == "W":
            st_out[i].stats.sac["cmpinc"] += 90.

    return st_out


def zerophase_chebychev_lowpass_filter(tr, freqmax):
    """
    Custom Chebychev type two zerophase lowpass filter useful for decimation
    filtering. This filter is stable up to a reduction in frequency with a
    factor of 10. If more reduction is desired, simply decimate in steps.
    Partly based on a filter in ObsPy.

    Filter parameters:
        rp: maximum ripple of passband
        rs: attenuation of stopband
        ws: stop band frequency
        wp: pass band frequency

    .. note::
        Should be replaced once ObsPy has a proper decimation filter.

    .. note::
        This function affects the trace in place!

    :type trace: obspy.core.trace.Trace
    :param trace: The trace to be filtered.
    :type freqmax: float
    :param freqmax: The desired lowpass frequency.
    """
    rp, rs, order = 1, 96, 1e99
    ws = freqmax / (tr.stats.sampling_rate * 0.5)
    wp = ws
    order, wn = None, None

    while True:
        if order <= 12:
            break
        wp *= 0.99
        order, wn = signal.cheb2ord(wp=wp, ws=ws, gpass=rp, gstop=rs,
                                    analog=False)

    b, a = signal.cheby2(N=order, rs=rs, Wn=wn, btype="low", analog=0,
                         output="ba")

    # Apply twice to get rid of the phase distortion.
    tr.data = signal.filtfilt(b, a, tr.data)
