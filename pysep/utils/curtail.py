"""
Utilities to curtail station lists based on source receiver parameters or to
curtail streams due to missing data etc.
"""
import numpy as np
from obspy import Stream
from obspy.geodetics import gps2dist_azimuth

from pysep import logger
from pysep.utils.fmt import get_codes


def curtail_by_station_distance_azimuth(event, inv, mindistance_km=0.,
                                        maxdistance_km=1E6, minazimuth=0.,
                                        maxazimuth=360.):
    """
    Remove stations that are greater than a certain distance from event
    Replaces the old `sta_limit_distance` function

    :type event: obspy.core.event.Event
    :param event: Event object to get location from
    :type inv: obspy.core.inventory.Inventory
    :param inv: inventory object to get locations from
    :type mindistance_km: float
    :param mindistance_km: minimum acceptable source-receiver distance in km
    :type maxdistance_km: float
    :param maxdistance_km: maximum acceptable source-receiver distance in km
    :type minazimuth: float
    :param minazimuth: minimum acceptable azimuth in deg
    :type maxazimuth: float
    :param maxazimuth: maximum acceptable azimuth in deg
    :rtype: obspy.core.inventory.Inventory
    :return: a curtailed inventory object which has had stations removed for
        unacceptable distance and azimuth values
    """
    event_latitude = event.preferred_origin().latitude
    event_longitude = event.preferred_origin().longitude

    remove_for_distance, remove_for_azimuth = [], []
    distances, azimuths = {}, {}
    for net in inv:
        for sta in net:
            # e.g., GR.FUR
            netsta_code = f"{net.code}.{sta.code}"
            dist_m, az, baz = gps2dist_azimuth(lat1=event_latitude,
                                               lon1=event_longitude,
                                               lat2=sta.latitude,
                                               lon2=sta.longitude
                                               )
            dist_km = dist_m / 1E3
            if not (mindistance_km <= dist_km <= maxdistance_km):
                remove_for_distance.append(netsta_code)
            elif not (minazimuth <= az <= maxazimuth):
                remove_for_azimuth.append(netsta_code)
            else:
                distances[netsta_code] = dist_km
                azimuths[netsta_code] = az

    logger.info(f"{len(remove_for_distance)} traces outside distance "
                f"bounds [{mindistance_km}, {maxdistance_km}]km"
                )
    if remove_for_distance:
        for remove in remove_for_distance:
            net, sta = remove.split(".")
            inv = inv.remove(network=net, station=sta)
        logger.debug(f"stations removed:\n{remove_for_distance}")

    logger.info(f"{len(remove_for_azimuth)} traces outside azimuth bounds "
                f"[{minazimuth}, {maxazimuth}]deg"
                )
    if remove_for_azimuth:
        for remove in remove_for_azimuth:
            net, sta = remove.split(".")
            inv = inv.remove(network=net, station=sta)
        logger.debug(f"stations removed: {remove_for_azimuth}")

    return inv


def quality_check_waveforms_before_processing(st, remove_clipped=True):
    """
    Quality assurance to deal with bad data before running the
    preprocessing steps. Replaces: `do_waveform_QA`

    :type st: obspy.core.stream.Stream
    :param st: Stream object to pass through QA procedures
    :type remove_clipped: bool
    :param remove_clipped: boolean flag to turn on/off amplitude clipping check
    """
    st_out = st.copy()

    st_out = rename_channels(st_out)
    st_out = remove_stations_for_missing_channels(st_out)  # LL network ONLY!
    st_out = remove_traces_for_bad_data_types(st_out)
    if remove_clipped:
        st_out = remove_for_clipped_amplitudes(st_out)

    return st_out


def quality_check_waveforms_after_processing(st,
                                             remove_insufficient_length=True):
    """
    Quality assurance to deal with bad data after preprocessing, because
    preprocesing step will merge, filter and rotate data.
    Replaces: `do_waveform_QA`

    :type st: obspy.core.stream.Stream
    :param st: Stream object to pass through QA procedures
    :type remove_insufficient_length: bool
    :param remove_insufficient_length: boolean flag to turn on/off insufficient
        length checker
    """
    st_out = st.copy()

    if remove_insufficient_length:
        st_out = remove_stations_for_insufficient_length(st_out)

    return st_out


def remove_traces_for_bad_data_types(st):
    """
    Removed traces from a Stream that have unexpected data types. This might
    occur if e.g., you wildcard the channel and end up grabbing LOG data, which
    uses letters.

    :type st: obspy.core.stream.Stream
    :param st: Stream to check clipping for
    :rtype st: obspy.core.stream.Stream
    :return st: curtailed stream with clipped traces removed
    """
    # NumPy data types smart enough that int32 or int64 will match to 'int'
    acceptable_data_types = [np.integer, np.floating]
    st_out = st.copy()
    for tr in st_out[:]:
        for dtype in acceptable_data_types:
            if np.issubdtype(tr.data.dtype, dtype):
                break
        else:
            logger.warning(f"{tr.get_id()} bad data type: {tr.data.dtype}")
            st_out.remove(tr)
    return st_out


def remove_traces_w_masked_data(st):
    """
    Merge operations may produce masked arrays which are data streams with
    gaps in them. Remove these from the stream
    """
    st_out = st.copy()
    for tr in st_out:
        if np.ma.is_masked(tr.data):
            logger.warning(f"{tr.get_id()} has data gaps, removing")
            st_out.remove(tr)
    return st_out


def remove_for_clipped_amplitudes(st):
    """
    Removed stations with clipped amplitudes
    replaces `clipping_handler.remove_clipped`
    TODO where is that clip factor coming from?

    :type st: obspy.core.stream.Stream
    :param st: Stream to check clipping for
    :rtype st: obspy.core.stream.Stream
    :return st: curtailed stream with clipped traces removed
    """
    st_out = st.copy()
    clip_factor = 0.8 * ((2 ** (24 - 1)) ** 2) ** 0.5  # For a 24-bit signal
    for tr in st_out[:]:
        # Figure out the if any amplitudes are clipped
        if (tr.data[np.abs(tr.data**2)**0.5 > clip_factor]).any():
            logger.info(f"removing {tr.get_id()} for clipped amplitudes")
            st_out.remove(tr)
    return st_out


def rename_channels(st):
    """
    Rename channels which intermix location names with channel names,
    For example: BHX00 -> BHX.00
    We are assuming here that channel codes are either: '00' or '10'
    Historically this is to differentiate STS-1 (00) and STS-2 (10)

    Relevant reading:
    https://ds.iris.edu/ds/newsletter/vol1/no1/1/
                        specification-of-seismograms-the-location-identifier/

    TODO old code strips channels down to 3 letters if they're 4. But
         can't we have 4 letter channel names? NZ does this.
    TODO Do we only expect location codes to be appended to channels?

    :type st: obspy.core.stream.Stream
    :param st: Stream to check incorrect channel naming for
    :rtype st: obspy.core.stream.Stream
    :return st: Stream with renamed channels and locations
    """
    logger.info("cleaning up channel naming")
    st_out = st.copy()
    for tr in st_out:
        for location_code in ["00", "10", "20"]:
            if tr.stats.channel.endswith(location_code):
                channel_code = tr.stats.channel
                tr.stats.location = location_code
                tr.stats.channel = channel_code.strip(location_code)
                logger.debug(f"stripping location code from channel name: "
                             f"{channel_code} -> {tr.stats.channel}")
    return st_out


def remove_stations_for_missing_channels(st, required_number_channels=3,
                                         networks="LL"):
    """
    Remove LLNL stations (network=='LL') with missing channels.

    LLNL data is already problematic, so if there are signs of too many
    issues / problems for a given station then remove that station.

    :type st: obspy.core.stream.Stream
    :param st: Stream to check missing channels for
    :type required_number_channels: int
    :param required_number_channels: expected channels for each station
    :type networks: str
    :param networks: comma-separated list of network codes to check. This
        defaults to 'LL' because this function was meant to parse through
        LLNL data
    """
    st_out = st.copy()
    check_networks = networks.split(",")  # allows function to be more general
    codes = get_codes(st_out, choice="location", up_to=True)

    for code in codes:
        net, sta, loc = code.split(".")
        if net.upper() not in check_networks:
            continue
        st_select = st_out.select(network=net, station=sta,
                                  location=loc, channel="*")
        cha_codes = get_codes(st_select, choice="channel", up_to=False)
        if len(cha_codes) < required_number_channels:
            logger.debug(f"{code}.* returned {len(cha_codes)} channels when "
                         f"the required amount is {required_number_channels}, "
                         f"removing")
            for tr in st_select:
                st_out.remove(tr)

    return st_out


def remove_stations_for_insufficient_length(st):
    """
    Remove stations if the length does not match the mode of all other lengths
    in the stream, which is assumed to be the expected length

    :type st: obspy.core.stream.Stream
    :param st: Stream to check for data gaps and insufficient start and
        end times
    """
    st_out = st.copy()

    # Assuming that the mode of the stream lengths represents the actual value
    stream_lengths = [tr.stats.endtime - tr.stats.starttime for tr in st_out]
    vals, counts = np.unique(stream_lengths, return_counts=True)
    expected_length = vals[np.argmax(counts)]
    logger.debug(f"assuming that the expected stream length is: "
                 f"{expected_length}s")
    for tr, length in zip(st_out[:], stream_lengths):
        if length < expected_length:
            logger.debug(f"{tr.get_id()} has unexpected time length of "
                         f"{length}s, removing")
            st_out.remove(tr)

    return st_out


def subset_streams(st_a, st_b):
    """
    Given two streams of data, check if they have the same length. IF they do,
    return the streams. If they don't, subset the streams so they have the same
    lengths and the same station ids.

    :type st_a: obspy.core.stream.Stream
    :param st_a: stream A to check
    :type st_b: obspy.core.stream.Stream
    :param st_b: stream B to check
    :rtype: tuple of Streams
    :return: curtailed (or not) streams in the same order as input
    """
    def get_ids(st):
        """
        Grab "network.station.location.comp" from obspy stream. because 
        synthetics may have different channels (e.g,. HXZ vs HHZ) in their name 
        so we can't match the whole station ID using Trace.get_id()
        """
        list_out = []
        for tr in st:
            s = tr.stats
            list_out.append(
                    f"{s.network}.{s.station}.{s.location}.{s.component}"
                    )
        return set(list_out)


    if len(st_a) != len(st_b):
        logger.warning(f"stream lengths don't match {len(st_a)} != {len(st_b)} "
                       f"will subset to the shorter length")
    else:
        return st_a, st_b

    st_a_out = Stream()
    st_b_out = Stream()

    # Collect all unique IDs so that we can use them for identification
    sta_ids_a = get_ids(st_a)
    sta_ids_b = get_ids(st_b)
    common_ids = sta_ids_a.intersection(sta_ids_b)

    logger.debug(f"stream subset removes "
                 f"{len(sta_ids_a) - len(common_ids)} traces from `st_a`")
    logger.debug(f"stream subset removes "
                 f"{len(sta_ids_b) - len(common_ids)} traces from `st_b`")

    for station_id in common_ids:
        net, sta, loc, comp = station_id.split(".")
        st_a_out += st_a.select(network=net, station=sta, location=loc, 
                                component=comp)
        st_b_out += st_b.select(network=net, station=sta, location=loc,
                                component=comp)

    assert(len(st_a_out) == len(st_b_out)), (
        f"station subsetting failed to return the same number of streams. " 
        f"check that your data does not contain multiple traces for a single "
        f"component as this can lead to this error message"
    )

    return st_a_out, st_b_out
