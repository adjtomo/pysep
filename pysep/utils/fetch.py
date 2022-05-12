"""
Grab information from external webservices or databases
"""
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees, gps2dist_azimuth



def fetch_bulk_station_list_taup(event, inv, distances, phases=None, taupmodel="ak135",
                              reference_time=None, seconds_before_ref=0,
                              seconds_after_ref=0):
    """

    :param stations:
    :param event:
    :param phases:
    :param phase_window:
    :param taupmodel:
    :param reftime:
    :param tbefore_sec:
    :param tafter_sec:
    :return:
    """
    # Find P and S arrival times using TauP
    if phases is not None:
        taup_model = TauPyModel(model=taupmodel)
        for net in inv:
            for sta in net:
                dist_m, az, baz = gps2dist_azimuth(
                    lat1=event.preferred_origin.latitude,
                    lon1=event.preferred_origin.longitude,
                    lat2=sta.latitude, lon2=sta.longitude
                )
                dist_km = dist_m / 1E3
                dist_deg = kilometer2degrees(dist_km, radius=6371)
                depth_km = event.preferred_origin.depth / 1E3
                arrivals = taup_model.get_travel_times(
                    source_depth_in_km=depth_km,
                    distance_in_degree=dist_deg,
                    phase_list=phases
                )
    #
    #
    #             try:
    #                 if Phase2arrivals[0].time < Phase1arrivals[0].time:
    #                     # You are assuming that the first index is the first arrival.  Check this later.
    #                     t1s.append(event.origins[0].time + Phase2arrivals[0].time - tbefore_sec)
    #                     t2s.append(event.origins[0].time + Phase1arrivals[0].time + tafter_sec)
    #                 else:
    #                     t1s.append(event.origins[0].time + Phase1arrivals[0].time - tbefore_sec)
    #                     t2s.append(event.origins[0].time + Phase2arrivals[0].time + tafter_sec)
    #             except:
    #                 t1s.append(reftime - tbefore_sec)
    #                 t2s.append(reftime + tafter_sec)
    #
    # else:
    #     t1s = [reftime - tbefore_sec]
    #     t2s = [reftime + tafter_sec]

    return t1s, t2s