"""
Utilities to curtail station lists based on source receiver parameters or to
curtail streams due to missing data etc.
"""
import numpy as np
from obspy.geodetics import gps2dist_azimuth

from pysep import logger


def curtail_by_station_distance_azimuth(event, inv, min_dist=0, max_dist=1E6,
                                        min_az=0, max_az=360):
    """
    Remove stations that are greater than a certain distance from event
    Replaces the old `sta_limit_distance` function

    TODO set azimuth values % 360 to not allow negative azimuth values
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
            if max_dist < dist_km < min_dist:
                net.remove(sta)
                remove_for_distance.append(netsta_code)
            elif max_az < az < min_az:
                net.remove(sta)
                remove_for_azimuth.append(netsta_code)
            else:
                distances[netsta_code] = dist_km
                azimuths[netsta_code] = az
        logger.info(f"{len(remove_for_distance)} stations removed for "
                    f"distance outside bounds [{min_dist}, {max_dist}]"
                    )
        logger.debug(f"station list: {[_ for _ in remove_for_distance]}")

        logger.info(f"{len(remove_for_azimuth)} stations removed for "
                    f"distance outside bounds [{min_az}, {max_az}]"
                    )
        logger.debug(f"station list: {[_ for _ in remove_for_azimuth]}")

    return inv, distances, azimuths


def quality_check_waveforms(st):
    """
    Quality assurance to deal with bad data. Replaces: `do_waveform_QA`
    """
    st_out = st.copy()

    st_out = _rename_channels(st_out)
    st_out = _remove_stations_for_missing_channels(st_out)

    logger.debug("filling any data gaps by interpolating missing values")
    st_out.merge(fill_value="interpolate")

    st_out = _remove_stations_for_insufficient_length(st_out)

    return st_out


def _rename_channels(st):
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


def _remove_stations_for_missing_channels(st, required_number_channels=3):
    """
    Remove LLNL stations (network=='LL') with missing channels.

    LLNL data is already problematic, so if there are signs of too many
    issues / problems for a given station then remove that station.
    """
    st_out = st.copy()
    check_networks = ["LL"]  # allows this function to be more general

    # Traverse the stream to get full station codes
    codes = [tr.get_id() for tr in st_out]
    net_codes = set([code.split(".")[0] for code in codes])  # unique networks
    sta_codes = set([code.split(".")[1] for code in codes])  # unique stations
    loc_codes = set([code.split(".")[2] for code in codes])  # unique locations

    for net_code in net_codes:
        if net_code.upper() not in check_networks:
            continue
        for sta_code in sta_codes:
            for loc_code in loc_codes:
                st_select = st_out.select(network=net_code, sttation=sta_code,
                                          location=loc_code, channel="*")
                if len(st_select) <= required_number_channels:
                    logger.debug(f"{net_code}.{sta_code}.{loc_code}.* returned "
                                 f"{len(st_select)} channels when the required "
                                 f"amount is {required_number_channels}, "
                                 f"removing")
                    for tr in st_select:
                        st_out.remove(tr)

    return st_out


def _remove_stations_for_insufficient_length(st):
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

