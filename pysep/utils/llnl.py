"""
Specific functions related to Lawrence Livermore National Lab waveform and
catalog data
"""
from obspy import UTCDateTime

from pysep import logger


def get_llnl_catalog():
    """
    Read in LLNL catalog YAML file and return dict of origintimes and event tags
    """
    pass


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
        * HOYA LL.LH and LL.BB require a sign flip. Otherwise the surface
            waveform inversions produce a -ISO result.
        * HOYA LL.LH and LL.BB body do not require flip.
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
                scale_factor_llnl = llnl_scale_factors[tr.stats.station]
            except KeyError:
                logger.warning(f"unexpected LLNL station: {tr.stats.station}")
                continue
            tr.data = tr.data * scale_factor_llnl
            if hasattr(tr.stats.sac, "scale"):
                tr.stats.sac["scale"] *= scale_factor_llnl
            else:
                tr.stats.sac["scale"] = scale_factor_llnl

    return st_out
