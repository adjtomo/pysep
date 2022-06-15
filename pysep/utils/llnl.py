"""
Specific functions related to Lawrence Livermore National Lab waveform and
catalog data
"""
from obspy import UTCDateTime

from pysep import logger


def return_llnl_catalog_origin_times(event=None):
    """
    Returns a static LLNL catalog which provides origin times and event tags

    :type event: str
    :param event: returns only a specific event based on the event names
    :rtype: dict
    :return: a dictionary where keys are event names (uppercase) and values
        are catalog origin times (UTCDateTime)
    """
    llnl_catalog = {
        # Explosions
        "KERNVILLE": "1988-02-15T18:10:00.09",
        "AMARILLO": "1989-06-27T15:30:00.02",
        "DISKO_ELM": "1989-09-14T15:00:00.10",
        "HORNITOS": "1989-10-31T15:30:00.09",
        "BARNWELL": "1989-12-08T15:00:00.09",
        "METROPOLIS": "1990-03-10T16:00:00.08",
        "BULLION": "1990-06-13T16:00:00.09",
        "AUSTIN": "1990-06-21T18:15:00.00",
        "HOUSTON": "1990-11-14T19:17:00.07",
        "COSO": "1991-03-08T21:02:45.08",
        "BEXAR": "1991-04-04T19:00:00.00",
        "HOYA": "1991-09-14T19:00:00.08",
        "LUBBOCK": "1991-10-18T19:12:00.00",
        "BRISTOL": "1991-11-26T18:35:00.07",
        "JUNCTION": "1992-03-26T16:30:00.00",
        "HUNTERS_TROPHY": "1992-09-18T17:00:00.08",
        "DIVIDER": "1992-09-23T15:04:00.00",
        "ATRISCO": "1982-08-05T14:00:00",  # Explosion not in Ford2009
        # Earthquakes
        "LITTLE_SKULL_MAIN": "1992-06-29T10:14:23.18",
        "LITTLE_SKULL_AFTERSHOCK": "1992-07-05T06:54:13.52",
        "TIMBER_MOUNTAIN": "1995-07-31T12:34:47.35",
        "GROOM_PASS": "1997-04-26T01:49:35.20",
        "CALICO_FAN": "1997-09-12T13:36:54.94",
        "WARM_SPRINGS": "1998-12-12T01:41:32.00",
        "FRENCHMAN_FLAT_1": "1999-01-23T03:00:32.00",
        "FRENCHMAN_FLAT_2": "1999-01-27T10:44:23.30",
        "LITTLE_SKULL": "2002-06-14T12:40:45.36",
        # Earthquakes in Ford, but not in LLNL database
        "AMARGOSA": "1996-09-05T08:16:55.40",
        "INDIAN_SPRINGS": "1997-06-14T19:48:19.45",
        "RALSTON": "2007-01-24T11:30:16.099",
        # Mine collapses
        "ATRISCO_HOLE": "1982-08-05T14:21:38.000",
        "TRONA_MINE_1": "1995-02-03T15:26:10.690",
        "TRONA_MINE_2": "2000-01-30T14:46:51.310"
        }

    llnl_catalog = {key: UTCDateTime(val) for key, val in llnl_catalog.items()}
    if event is None:
        try:
            llnl_catalog = {event.upper(): llnl_catalog[event.upper()]}
        except KeyError:
            raise KeyError(f"Invalid event tag: {event.upper()}; acceptable: "
                           f"{llnl_catalog.keys()}")

    return llnl_catalog


def scale_llnl_waveform_amplitudes(st):
    """
    Rescale amplitudes for LLNL data. The scales are different for different
    channels.

    .. note::
        Old factors
        BB = 4.0e-4
        HF = 1.8e-2
        HF = 1.8e-4

    .. note::
        scales based on tests with events HOYA, BEXAR
        * LL.HOYA..LH? and LL.HOYA..BB? require a sign flip. Otherwise the 
            surface waveform inversions produce a -ISO result.
        * LL.HOYA..LH? and LL.HOYA..BB? body do not require flip.

    :type st: obspy.core.stream.Stream
    :param st: Stream containing LLNL waveforms
    :rtype: obspy.core.stream.Stream
    :return: Stream with scaled LLNL waveforms
    """
    st_out = st.copy()

    # These are STATION dependent amplitude factors
    llnl_scale_factors = {
        "LH": -1 * 1E-2,  # note sign flip here
        "BB": -1 * 1E-9,  # note sign flip here
        "XX": -1 * 1E-9,  # note sign flip here
        "HF": 1E-9,
        "VB": 1E-9,
        "EH": 1E-9,  # originally 1E-7 but some amps too large, e.g,. FF2
        "HH": 1E-9,  # originally 1E-10 but for FF2 some amps too small
        "BH": 1E-9,
        "SH": 1E-9,
        "HG": 1E-5,
    }
    for tr in st_out:
        if tr.stats.network == "LL":
            try:
                channel_prefix = tr.stats.channel[:-1]
                scale_factor_llnl = llnl_scale_factors[channel_prefix]
            except KeyError:
                logger.warning(f"unexpected LLNL station: {tr.stats.station}")
                continue
            tr.data = tr.data * scale_factor_llnl
            if hasattr(tr.stats.sac, "scale"):
                tr.stats.sac["scale"] *= scale_factor_llnl
            else:
                tr.stats.sac["scale"] = scale_factor_llnl

    return st_out

