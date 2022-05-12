"""
Utilities to curtail station lists based on source receiver parameters
"""
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