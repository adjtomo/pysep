"""
Grab information from external webservices or databases
"""
import os
from glob import glob
from pathlib import Path
from obspy.taup import TauPyModel
from obspy.taup.taup_create import build_taup_model
from obspy.geodetics import kilometer2degrees, gps2dist_azimuth

from pysep import logger



def get_taup_arrivals(st=None, event=None, inv=None, phase_list=("ttall",),
                      model="ak135", network=None, station=None):
    """
    Retrieve phase arrival times using TauP/TauPy based on event and station
    locations, either using SAC headers or using event and inventory metadata.

    .. note::

        Requires EITHER `st` with appropriate SAC headers OR
        `event` and `inv` defining event and station locations.

    Available model names can be found here:
        https://docs.obspy.org/packages/obspy.taup.html
    OR internally at:
        pysep/pysep/data/taup_models

    :type st: obspy.core.stream.Stream
    :param st: Stream with appended SAC headers which will be used to gather
        TaupP arrivals from a given `model`
    :type event: obspy.core.event.Event
    :param event: Event object to get location from
    :type inv: obspy.core.inventory.Inventory
    :param inv: inventory object to get locations from
    :type phase_list: list of str
    :param phase_list: phases to search TauP for.  Defaults to 'ttall', which is
        ObsPy's default for getting all phase arrivals
    :type model: str
    :param model: Model to query TauP with. Looks at ObsPy and PySEP models
    :type network: str
    :param network: query only a specific network
    :type station: str
    :param station: query only a specific station
    :rtype: dict
    :return: a dictionary where keys are trace ID consisting of network and
        station (NN.SSS), and values are TauP Arrivals() instances which
        contain information about phase arrivals
    """
    if st is not None:
        logger.debug(f"using SAC headers to fetch phase arrivals '{phase_list}'"
                     f"from model '{model}'")
    else:
        assert(event is not None and inv is not None), (
            f"`get_taup_arrivals` requires `st` w/ SAC headers OR `event` and " 
            f"`inv` defining event and stations."
        )
        logger.debug(f"using event and station metadata to fetch phase "
                     f"arrivals '{phase_list}' from model '{model}'")

    phase_dict = {}
    taup_model = _get_or_build_taup_model(model)
    taup_func = taup_model.get_travel_times

    # OPTION 1: Retrieve event and station metadata from SAC headers
    if st is not None:
        for tr in st:
            # Allow selecting for specific networks and stations
            net, sta, *_ = tr.get_id().split(".")
            trace_id = f"{net}.{sta}"
            if network and net != network:
                continue
            if station and sta != station:
                continue

            sac_header = tr.stats.sac
            if "evdp" not in sac_header or "gcarc" not in sac_header:
                logger.debug(f"skip {tr.get_id()} phase arr., no 'evdp' or "
                             f"'gcarc")
                continue
            depth_km = sac_header["evdp"]  # units: km
            dist_deg = sac_header["gcarc"]
            while True:
                try:
                    arrivals = taup_func(source_depth_in_km=depth_km,
                                         distance_in_degree=dist_deg,
                                         phase_list=phase_list)
                    if not arrivals:
                        logger.debug(f"no TauP arrivals for model {model} and "
                                     f"phases: {phase_list}")
                    phase_dict[trace_id] = arrivals
                    break
                # This will through a ValueError for invalid phase names
                except ValueError as e:
                    # Assuming that the error message lists the phase name last
                    logger.warning(f"'{e}' invalid phase name, removing from "
                                   f"phase list")
                    del_phase = str(e).split(" ")[-1]
                    phase_list.remove(del_phase)
                    pass
    # OPTION 2: Retrieve event and station metadata from metadata objects
    else:
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
                while True:
                    try:
                        arrivals = taup_func(source_depth_in_km=depth_km,
                                             distance_in_degree=dist_deg,
                                             phase_list=phase_list)
                        phase_dict[code] = arrivals
                        break
                    # This will through a ValueError for invalid phase names
                    except ValueError as e:
                        # Assuming that the error message lists phase name last
                        logger.warning(f"'{e}' invalid phase name, removing "
                                       f"from phase list")
                        del_phase = str(e).split(" ")[-1]
                        phase_list.remove(del_phase)
                        pass

    return phase_dict


def _get_or_build_taup_model(model):
    """
    Attempts to retrieve a TauP model from the internal ObsPy directory.
    If nothing is found, will attempt to search the internal PySEP taup_model
    data directory for matching velocity models which can be built

    :type model: str
    :param model: name of model to search for
    :rtype: obspy.taup.TauPyModel
    :return: model matching model name
    :raises FileNotfoundError: If no matching file is found
    """
    try:
        taup_model = TauPyModel(model=model)
    except FileNotFoundError as e:
        logger.debug(f"TauP model {model} not found in ObsPy, searching PySEP")
        # Try to build TauP model from internally stored data. Pathing relative
        # to this file
        pysep_dir = Path(__file__).absolute().parent.parent

        # ObsPy expects model files with extension: 'tvel' or 'nd'
        taup_model_glob = os.path.join(pysep_dir, "data", "taup_models", "*.{}")
        taup_model_files = glob(taup_model_glob.format("tvel"))
        taup_model_files.append(taup_model_glob.format("nd"))

        # Get the base model names, e.g., 'ak_scak'
        model_names = [os.path.splitext(os.path.basename(_))[0] for _ in
                       taup_model_files]
        if model in model_names:
            logger.info(f"matching model found {model}, building TauP model")
            model_file = taup_model_files[model_names.index(model)]
            build_taup_model(model_file)
            taup_model = TauPyModel(model=model)
        else:
            raise FileNotFoundError(f"TauP model {model} not found within "
                                    f"ObsPy or PySEP: {e}")

    return taup_model
