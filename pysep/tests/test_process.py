"""
Test the preprocessing functions which include standardization and rotation

.. note::
    Pysep.preprocess() is almost all ObsPy functions so we don't test that as
    we expect ObsPy to be testing their own functionality
"""
import pytest
import numpy as np
from obspy import read, read_events, read_inventory

from pysep import logger, Pysep
from pysep.utils.cap_sac import append_sac_headers
from pysep.utils.process import (merge_and_trim_start_end_times,
                                 resample_data, format_streams_for_rotation,
                                 rotate_to_uvw, append_back_azimuth_to_stats)

logger.setLevel("DEBUG")  # DELETEME


@pytest.fixture
def test_event():
    """
    Pre-saved event object
    """
    return read_events("./test_data/test_event.xml")[0]


@pytest.fixture
def test_inv():
    """
    Pre-saved event object
    """
    return read_inventory("./test_data/test_inv.xml")


@pytest.fixture
def test_st(test_event, test_inv):
    """
    Pre-saved event object
    """
    st = read("./test_data/test_stream.mseed")
    st = append_sac_headers(st=st, inv=test_inv, event=test_event)
    return st


def test_merge_and_trim_start_end_times(test_st):
    """
    Make sure we can trim per-station waveforms to their shortest component
    and also merge existing waveforms, dropping those with data gaps
    """
    st = merge_and_trim_start_end_times(test_st)
    assert(len(test_st) - len(st) == 4)


def test_resample_data(test_st):
    """
    Resample data to a lower sampling rate
    """
    new_sampling_rate = test_st[0].stats.sampling_rate // 2
    st_int = resample_data(test_st, resample_freq=new_sampling_rate,
                           method="interpolate")
    st_res = resample_data(test_st, resample_freq=new_sampling_rate,
                           method="resample")

    actual_val = test_st[-1].max()
    int_val = st_int[-1].max()
    res_val = st_res[-1].max()

    # Checking that we actually resampled data which will affect amplitude
    assert(actual_val != int_val)
    assert(actual_val != int_val)

    # Make sure we didn't mess up amplitude too badly
    assert(actual_val / int_val > 0.8)
    assert(actual_val / res_val > 0.8)


def test_format_streams_for_rotation(test_st):
    """
    Check that we can add null traces if we're missing data streams
    """
    st = test_st.copy()
    # One horizontal and vertical
    for tr in st.select(station="ATKA", component="E")[:]:
        st.remove(tr)
    # Only one horizontal
    for tr in st.select(station="BESE", component="Z")[:]:
        st.remove(tr)
    for tr in st.select(station="BESE", component="N")[:]:
        st.remove(tr)
    # Only horizontals
    for tr in st.select(station="ALPI", component="Z")[:]:
        st.remove(tr)

    st_out = format_streams_for_rotation(st)
    assert(len(st_out) == 6)
    assert(np.all((st_out.select(station="ATKA", component="N")[0].data == 0)))


def test_rotate_to_uvw(test_st):
    """
    Test orientation UVW rotation
    TODO test function and get better assertions
    """
    st = test_st.select(station="ATKA")
    st.merge()
    st_out = rotate_to_uvw(st)
    for tr in st_out:
        assert(tr.stats.component in ["U", "V", "W"])


def test_append_back_azimuth_to_stats(test_st, test_inv, test_event):
    """
    Check that back azimuth values are appended correctly
    """
    # ensure that test data doesn't have any back azimuth stats
    for tr in test_st:
        del tr.stats.back_azimuth
    assert("back_azimuth" not in test_st[0].stats)

    st = append_back_azimuth_to_stats(st=test_st, inv=test_inv,
                                      event=test_event)
    check_dict = {"ATKA": 45.514, "BESE": 297.444, "ALPI": 335.115}
    for tr in st:
        assert(check_dict[tr.stats.station] ==
               pytest.approx(tr.stats.back_azimuth, 3))


def test_rotate_streams_fail(test_st, test_inv, test_event):
    """
    Ensure that stream rotation does not go ahead if no back azimuth values are
    specified
    """
    sep = Pysep(log_level="CRITICAL", rotate=["ZNE", "RTZ"])
    for tr in test_st:
        del tr.stats.back_azimuth
    assert("back_azimuth" not in test_st[0].stats)

    sep.st = test_st
    sep.inv = test_inv
    sep.event = test_event
    sep.st = sep.preprocess()  # make sure that streams are formatted correctly
    with pytest.raises(TypeError):
        sep.rotate_streams()


def test_rotate_streams(test_st, test_inv, test_event):
    """
    Ensure that stream rotation works as expected after we format SAC headers
    """
    sep = Pysep(log_level="CRITICAL", rotate=["ZNE", "RTZ"])
    for tr in test_st:
        del tr.stats.back_azimuth
    assert("back_azimuth" not in test_st[0].stats)

    sep.st = test_st
    sep.inv = test_inv
    sep.event = test_event
    sep.st = append_sac_headers(sep.st, sep.event, sep.inv)
    sep.st = sep.preprocess()  # make sure that streams are formatted correctly
    st_rotated = sep.rotate_streams()
    assert(len(st_rotated) == 14)
    # Make sure components have been rotated
    for tr in st_rotated:
        assert(tr.stats.component in ["R", "T", "Z", "N", "E"])
