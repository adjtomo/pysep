"""
Grab information from external webservices or databases
"""
from obspy.taup import TauPyModel
from obspy.geodetics import kilometer2degrees, gps2dist_azimuth

from pysep import logger


def get_taup_arrivals(event, inv, phase_list=None, model="ak135", network=None,
                      station=None):
    """
    Retrieve phase arrival times using TauP/TauPy based on event and stations.
    By default only retrieve for first arrivals
    """
    phase_dict = {}
    if phase_list is None:
        phase_list = ["P", "S"]

    taup_model = TauPyModel(model=model)
    logger.info(f"fetching arrivals, phases {phase_list} and model '{model}'")
    for net in inv.select(network=network, station=station):
        for sta in net:
            code = f"{net.code}.{sta.code}"
            dist_m, az, baz = gps2dist_azimuth(
                lat1=event.preferred_origin().latitude,
                lon1=event.preferred_origin().longitude,
                lat2=sta.latitude, lon2=sta.longitude
            )
            dist_km = dist_m / 1E3
            dist_deg = kilometer2degrees(dist_km, radius=6371)
            depth_km = event.preferred_origin().depth / 1E3  # units: m -> km
            arrivals = taup_model.get_travel_times(
                source_depth_in_km=depth_km,
                distance_in_degree=dist_deg,
                phase_list=phase_list
            )
            phase_dict[code] = arrivals

    return phase_dict
