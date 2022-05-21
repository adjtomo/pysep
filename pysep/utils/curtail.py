"""
Utilities to curtail station lists based on source receiver parameters or to
curtail streams due to missing data etc.
"""
import numpy as np
from obspy.geodetics import gps2dist_azimuth

from pysep import logger
from pysep.utils.fmt import get_codes


def curtail_by_station_distance_azimuth(event, inv, min_dist=0, max_dist=1E6,
                                        min_az=0, max_az=360):
    """
    Remove stations that are greater than a certain distance from event
    Replaces the old `sta_limit_distance` function
    :return:
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
            if not (min_dist <= dist_km <= max_dist):
                remove_for_distance.append(netsta_code)
            elif not (min_az <= az <= max_az):
                remove_for_azimuth.append(netsta_code)
            else:
                distances[netsta_code] = dist_km
                azimuths[netsta_code] = az

    logger.info(f"{len(remove_for_distance)} traces outside distance "
                f"bounds [{min_dist}, {max_dist}]km"
                )
    if remove_for_distance:
        for remove in remove_for_distance:
            net, sta = remove.split(".")
            inv = inv.remove(network=net, station=sta)
        logger.debug(f"stations removed:\n{remove_for_distance}")

    logger.info(f"{len(remove_for_azimuth)} traces outside azimuth bounds "
                f"[{min_az}, {max_az}]deg"
                )
    if remove_for_azimuth:
        for remove in remove_for_azimuth:
            net, sta = remove.split(".")
            inv = inv.remove(network=net, station=sta)
        logger.debug(f"stations removed: {remove_for_azimuth}")

    return inv, distances, azimuths


def quality_check_waveforms_before_processing(st):
    """
    Quality assurance to deal with bad data before running the
    preprocessing steps. Replaces: `do_waveform_QA`
    """
    st_out = st.copy()

    st_out = rename_channels(st_out)
    st_out = remove_stations_for_missing_channels(st_out)  # LL network ONLY!
    st_out = remove_for_clipped_amplitudes(st_out)

    return st_out


def quality_check_waveforms_after_processing(st):
    """
    Quality assurance to deal with bad data after preprocessing, because
    preprocesing step will merge, filter and rotate data.
    Replaces: `do_waveform_QA`
    """
    st_out = st.copy()

    st_out = remove_stations_for_insufficient_length(st_out)

    return st_out


def remove_for_clipped_amplitudes(st):
    """
    Removed stations with clipped amplitudes

    replaces `clipping_handler.remove_clipped`
    TODO where is that clip factor coming from?
    """
    st_out = st.copy()
    clip_factor = 0.8 * ((2 ** (24 - 1)) ** 2) ** 0.5  # For a 24-bit signal
    for tr in st_out[:]:
        # Figure out the if any amplitudes are clipped
        if len(tr.data[np.abs(tr.data**2)**0.5 > clip_factor]):
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
    """
    st_out = st.copy()

    # Assuming that the mode of the stream lengths represents the actual value
    stream_lengths = [tr.stats.endtime - tr.stats.starttime for tr in st_out]
    vals, counts = np.unique(stream_lengths, return_counts=True)
    expected_length = vals[np.argmax(counts)]
    logger.debug(f"assuming that the expected stream length is: "
                 f"{expected_length}s")

    for tr, length in zip(st_out[:], stream_lengths):
        if length != expected_length:
            logger.debug(f"{tr.get_id()} has insufficient time length of "
                         f"{length}s, removing")
            st_out.remove(tr)

    return st_out

