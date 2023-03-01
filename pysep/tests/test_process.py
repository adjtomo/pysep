"""
Test the preprocessing functions which include standardization and rotation

.. note::
    Pysep.preprocess() is almost all ObsPy functions so we don't test that as
    we expect ObsPy to be testing their own functionality
"""
import pytest
import random
import numpy as np
from obspy import read, read_events, read_inventory, Stream

from pysep import logger, Pysep
from pysep.utils.cap_sac import append_sac_headers
from pysep.utils.curtail import remove_traces_w_masked_data
from pysep.utils.process import (merge_gapped_data, trim_start_end_times,
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


def test_trim_start_end_times(test_st):
    """
    Make sure we can trim per-station waveforms to their shortest component
    and also merge existing waveforms, dropping those with data gaps
    """
    st = trim_start_end_times(test_st)
    st = remove_traces_w_masked_data(st)
    assert(len(test_st) - len(st) == 4)


def test_merge_data_gaps(test_st):
    """
    Make sure merging data gaps works for a few different options. Does not
    check values, just runs through the function
    """
    # ATKA already has gappy data on components E and N
    st_gap = test_st.select(station="ATKA")
    assert(len(st_gap) == 5)

    for fill_value in [None, "mean", "interpolate", "latest", 0, 5.5]:
        st = merge_gapped_data(st_gap, fill_value=fill_value)
        st = remove_traces_w_masked_data(st)
        if fill_value in [None, 5.5]:
            assert(len(st) == 1)  # removed E and N from stations
        else:
            assert(len(st) == 3)  # successful merge


def test_merge_data_gap_fraction(test_st):
    """Try out the fraction percentage function"""
    # ATKA already has gappy data on components E and N
    st_gap = test_st.select(station="ATKA")
    assert(len(st_gap) == 5)

    st = merge_gapped_data(st_gap, fill_value=0, gap_fraction=0.01)
    st = remove_traces_w_masked_data(st)
    assert(len(st) == 1)


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
    sep = Pysep(log_level="DEBUG", rotate=["ZNE", "RTZ"])
    for tr in test_st:
        del tr.stats.back_azimuth
    assert("back_azimuth" not in test_st[0].stats)

    sep.st = test_st
    sep.inv = test_inv
    sep.event = test_event
    sep.st = sep.preprocess()  # make sure that streams are formatted correctly

    # No back aizmuth attribute found so streams will NOT be rotated
    st = sep.rotate_streams()
    components = set([tr.stats.component for tr in st])
    assert(not {"R", "T"}.issubset(components))

def test_rotate_streams(test_st, test_inv, test_event):
    """
    Ensure that stream rotation works as expected after we format SAC headers
    """
    sep = Pysep(log_level="DEBUG", rotate=["ZNE", "RTZ"])
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


def test_rotate_streams_enz(test_st, test_inv):
    """
    Test that streams rotate from ENZ even if components are already ENZ
    """
    sep = Pysep(log_level="DEBUG", rotate=["ENZ"])
    sep.st = test_st.select(station="ALPI").copy()
    st_check = test_st.select(station="ALPI").copy()

    sep.inv = test_inv.select(station="ALPI")
    sep.event = test_event

    # This rotates ->ZNE from ZNE. We expect the data to be the same (with some
    # error due to floating point rounding)
    st_rotated = sep.rotate_streams()

    for component in ["E", "N", "Z"]:
        st_rot = st_rotated.select(component=component)
        st_chk = st_check.select(component=component)
        for tr_rot, tr_check in zip(st_rot, st_chk):
            assert(pytest.approx(np.max(tr_rot.data - tr_check.data), 3) == 0)

    # Now we randomly assign azimuths and dips to stations and rotate. We 
    # expect data streams to have been rotated and therefore to be different
    n = 0
    for net in sep.inv:
        s = 0
        for sta in net:
            c = 0
            for cha in sta:
                sep.inv[n][s][c].azimuth = random.randint(1, 89)
                sep.inv[n][s][c].dip = random.randint(1, 89)
                c += 1
            s += 1
        n += 1
    st_rotated = sep.rotate_streams()

    for component in ["E", "N", "Z"]:
        st_rot = st_rotated.select(component=component)
        st_chk = st_check.select(component=component)

        for tr_rot, tr_check in zip(st_rot, st_chk):
            assert(not pytest.approx(
                np.max(tr_rot.data - tr_check.data), 3) != 0
                   )


def test_rotate_streams_12z():
    """
    Test that 12Z components can be rotated to ZNE
    """
    sep = Pysep(rotate=["ENZ"])
    sep.st = read("./test_data/test_12Z_data/stream.ms")
    sep.inv = read_inventory("./test_data/test_12Z_data/inv.xml")

    st_rotate = sep.rotate_streams()
    components = "".join(sorted([_.stats.component for _ in st_rotate]))
    assert(components == "ENZ")

