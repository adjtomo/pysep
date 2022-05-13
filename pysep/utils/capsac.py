"""
Utils to honor file formats from SAC and CAP
"""

def append_sac_headers(st, event, inv):
    """
    Append SAC headers to ObsPy streams
    :return:
    """
    for tr in st:
        net_code, sta_code, _, _= tr.get_id().split(".")
        sta = inv.select(network=net_code, station=sta_code)
        sac_header = {
            "evla": event.preferred_origin().latitude,
            "evlo": event.preferred_origin().longitude,
            "evdp": event.preferred_origin().depth / 1E3,  # depth in km
            "mag": event.preferred_magnitude(),
            "stla": sta.latitude,
            "stlo": sta.longitude,
            "stel": sta.elevation / 1E3,  # elevation in km
        }