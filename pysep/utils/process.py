"""
Utilities for processing waveforms
"""
import numpy as np
from scipy import signal

from obspy import Stream
from obspy.geodetics import gps2dist_azimuth
from pysep import logger
from pysep.utils.fmt import get_codes


def format_streams_for_rotation(st):
    """
    ObsPy requires specific channel naming to get stream rotation working.
    Comb through the stream and check that this is correct before passing
    stream to rotation

    :type st: obspy.core.stream.Stream
    :param st: Stream to format for rotation
    :rtype st: obspy.core.stream.Stream
    :return st: Stream that has been formatted for rotation
    """
    st_original = st.copy()
    logger.info(f"formatting stream of length {len(st)} for component rotation")
    codes = get_codes(st, choice="channel", suffix="?", up_to=True)

    st_out = Stream()
    for code in codes:
        net, sta, loc, cha = code.split(".")
        st = st_original.select(network=net, station=sta, location=loc,
                                channel=cha).copy()
        missing_vertical, missing_horizontal = _check_component_availability(st)

        if missing_vertical and missing_horizontal:
            logger.debug(f"{code} is missing vertical and at least 1 "
                         f"horizontal, skipping over")
            continue
        elif missing_vertical and not missing_horizontal:
            logger.debug(f"{code} is missing vertical, replacing vertical "
                         f"component with zeros")
            st = _append_null_trace(st, component="Z", cmpaz=0, cmpinc=-90.)
        elif not missing_vertical and missing_horizontal:
            logger.debug(f"{code} is missing at least 1 horizontal component, "
                         f"removing existing horizontal data and writing zeros "
                         f"to North and East components")
            for comp in ["E", "N", "1", "2"]:
                for tr in st.select(component=comp)[:]:
                    st.remove(tr)
                    logger.debug(f"{code} replacing existing component "
                                 f"'{comp}' with zeros")

            st = _append_null_trace(st, component="E", cmpaz=90., cmpinc=0.)
            st = _append_null_trace(st, component="N", cmpaz=0., cmpinc=0.)

        st_out += st

    logger.info(f"stream formatting returned {len(st_out)} traces")

    return st_out


def append_back_azimuth_to_stats(st, inv, event):
    """
    ObsPy rotation requires back azimuth information which is accessed through
    trace.stats.back_azimuth. This function appends back azimuth values to the
    stream by cacluating using metadata from the inventory and event objects.

    .. note::
        This was written and then I realized that `_append_sac_headers_trace`
        already does this, so it is unused, but left incase it's useful later.

    :type st: obspy.core.stream.Stream
    :param st: Stream with missing components
    :type inv: obspy.core.inventory.Inventory
    :param inv: optional user-provided inventory object which will force a
        skip over StationXML/inventory searching
    :type event: obspy.core.event.Event
    :param event: optional user-provided event object which will force a
        skip over QuakeML/event searching
    :rtype: obspy.core.stream.Stream
    :return: stream with back azimuth values appended
    """
    st_out = st.copy()
    ev_lat = event.preferred_origin().latitude
    ev_lon = event.preferred_origin().longitude
    for net in inv:
        for sta in net:
            _, _, baz = gps2dist_azimuth(lat1=ev_lat, lon1=ev_lon,
                                         lat2=sta.latitude, lon2=sta.longitude
                                         )
            # Apply back azimuth changes to stream in place
            for tr in st_out.select(network=net.code, station=sta.code):
                tr.stats.back_azimuth = baz
    return st_out


def _append_null_trace(st, component, cmpaz=None, cmpinc=None):
    """
    Create a copy of a trace with zeroed out data, used for rotation

    :type st: obspy.core.stream.Stream
    :param st: Stream with missing components
    :type component: str
    :param component: component of data to turn into a null trace
    :type cmpaz: float
    :param cmpaz: component azimuth for SAC header
    :type cmpinc: float
    :param cmpinc: component inclination for SAC header
    :rtype st: obspy.core.stream.Stream
    :return st: Stream that has had null traces appended
    """
    st_out = st.copy()
    trace_null = st_out[0].copy()
    trace_null.data = np.zeros(len(trace_null.data))
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

    :type st: obspy.core.stream.Stream
    :param st: Steam to check
    :rtype: tuple of bool
    :return: (is the vertical component missing?,
              is one or both horizontal components missing?)
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


def rotate_to_uvw(st):
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
    :param st: three-component stream to rotate the UVW
    :rtype: obspy.core.stream.Stream
    :return: Stream that has been rotated to UVW
    """
    st_out = st.copy()
    assert(len(st) == 3), f"UVW rotation requires 3 component stream input"

    xyz = np.array([[st_out[0].data],
                    [st_out[1].data],
                    [st_out[2].data]])
    xyz = np.asmatrix(xyz)

    rotation_matrix = np.matrix([[2, 0, np.sqrt(2)],
                                 [-1, np.sqrt(3), np.sqrt(2)],
                                 [-1, -np.sqrt(3), np.sqrt(2)]])
    rotation_matrix *= 1 / np.sqrt(6)

    # Theta is the angle between the axis and the horizontal plane
    theta = (np.arcsin(1 / np.sqrt(3))) * (180.0 / np.pi)

    uvw = rotation_matrix * xyz
    for i, comp in enumerate(["U", "V", "W"]):
        st_out[i].data = np.reshape(np.asarray(uvw[i]), uvw.shape[1])
        st_out[i].stats.component = comp
        st_out[i].stats.sac["kcmpnm"] = st_out[i].stats.channel
        st_out[i].stats.sac["cmpaz"] *= (120 * i)  # azimuth
        st_out[i].stats.sac["cmpinc"] -= theta  # inclination
        if comp == "W":
            st_out[i].stats.sac["cmpinc"] += 90.

    return st_out


def merge_and_trim_start_end_times(st):
    """
    Trim the maximum start and minumum end times for each station to get
    all waveforms to start and end at the same time. Also merges traces with
    the same ID

    Replaces old `trim_maxstart_minend`

    .. note::
        A note from the old code said we that interpolation and resampling
        are NOT intended functionalities, so we replace the old 'interpolate'
        function with a trim.

    :type st: obspy.core.stream.Stream
    :param st: stream to merge and trim start and end times for
    :rtype: obspy.core.stream.Stream
    :return: Stream that has been trimmed
    """
    st_edit = st.copy()

    logger.info("trimming start and end times on a per-station basis")
    st_out = Stream()
    # e.g., NN.SSS.LL.CC?  i.e., only getting per-station codes
    codes = get_codes(st=st_edit, choice="channel", suffix="?")

    for code in codes:
        # Subset the stream based on the N number of components
        net, sta, loc, cha = code.split(".")
        st_edit_select = st_edit.select(network=net, station=sta,
                                        location=loc, channel=cha)

        st_edit_select.merge()  # combining like trace IDs
        for tr in st_edit_select:
            if np.ma.is_masked(tr.data):
                logger.warning(f"{tr.get_id()} has data gaps, removing")
                st_edit_select.remove(tr)

        if st_edit_select:
            # Trim based on the most latest start and earliest end time, which
            # will expose stations with insufficient data
            max_start = max([tr.stats.starttime for tr in st_edit_select])
            min_end = min([tr.stats.endtime for tr in st_edit_select])

            st_edit_select.trim(starttime=max_start, endtime=min_end)
            st_out += st_edit_select

    return st_out


def resample_data(st, resample_freq, method="interpolate"):
    """
    Apply a lowpass filter and taper data before resampling it using
    one of a choice of resampling strategies

    TODO did not refactor `resample_cut` which doesn't seem to do anything but
        is included in old code

    :type st: obspy.core.stream.Stream
    :param st: Stream to resample data for
    :type resample_freq: float
    :param resample_freq: frequency to resample for
    :type method: str
    :param method: choice of ObsPy resampling strategy:
        - interpolate: interpolate data, requires anti-aliasing lowpass filter
            before interpolating
        - resample: reample data with fourier methods
    :rtype: obspy.core.stream.Stream
    :return: Stream that has been resampled and maybe filtered
    """
    if method not in ["interpolate", "resample"]:
        logger.warning(f"unexpected resample method {method}, defaulting to "
                       f"'interpolate'")
        method = "interpolate"

    logger.info(f"resampling data to sampling rate: {resample_freq}Hz with "
                f"method '{method}'")
    st_out = st.copy()
    for i, tr in enumerate(st_out[:]):
        target_nyquist = 0.5 * resample_freq
        current_nyquist = 0.5 * tr.stats.sampling_rate
        if target_nyquist < current_nyquist:
            try:
                logger.debug(f"{tr.get_id()} applying lowpass filter to prep "
                             f"for new nyquist freq: "
                             f"{current_nyquist} -> {target_nyquist}")
                # This function affects the trace in-place
                zerophase_chebychev_lowpass_filter(
                    tr=tr, freqmax=target_nyquist
                )
            except Exception as e:
                logger.warning(f"exception in lowpass filtering "
                               f"{tr.get_id()}, remove: {e}")
                logger.debug(e)
                st_out.remove(tr)
        else:
            tr.detrend("linear")
        tr.taper(max_percentage=0.01, type="hann")
        # Enforce that this data array is 'c'ontiguous
        tr.data = np.require(tr.data, requirements=["C"])
        if method == "interpolate":
            tr.interpolate(sampling_rate=resample_freq, method="lanczos", a=8)
        elif method == "resample":
            tr.resample(sampling_rate=resample_freq, strict_length=True)

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

    :type tr: obspy.core.trace.Trace
    :param tr: The trace to be filtered.
    :type freqmax: float
    :param freqmax: The desired lowpass frequency.
    """
    rp, rs, order = 1, 96, 1e99
    ws = freqmax / (tr.stats.sampling_rate * 0.5)
    wp = ws
    wn = None

    while True:
        if order <= 12:
            break
        wp *= 0.99
        order, wn = signal.cheb2ord(wp=wp, ws=ws, gpass=rp, gstop=rs,
                                    analog=False)
    b, a = signal.cheby2(N=order, rs=rs, Wn=wn, btype="low", analog=0,  # NOQA
                         output="ba")                                   # NOQA
    # Apply twice to get rid of the phase distortion.
    tr.data = signal.filtfilt(b, a, tr.data)


def estimate_prefilter_corners(tr):
    """
    Estimates corners for pre-deconvolution (instrument response removal)
    Essentially band limits the data based on the minimum and maximum allowable
    periods based on the sampling rate and overall waveform length.

    See also: https://ds.iris.edu/files/sac-manual/commands/transfer.html

    Replaces the old `get_pre_filt` function

    :type tr: obspy.core.trace.Trace
    :param tr: trace from which stats are taken
    :rtype: tuple of float
    :return: frequency corners (f0, f1, f2, f3)
    """
    # filtering constants
    fcut1_par = 4.0
    fcut2_par = 0.5

    f1 = fcut1_par / (tr.stats.endtime - tr.stats.starttime)
    nyquist_freq = tr.stats.sampling_rate / 2
    f2 = nyquist_freq * fcut2_par
    f0 = 0.5 * f1
    f3 = 2.0 * f2

    return f0, f1, f2, f3
