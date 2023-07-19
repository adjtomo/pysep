"""
Test the general utility functions of Pysep which includes quality
checking waveforms and formatting data with metadata

TODO add function testing clipped station, need to find example data
"""
import os
import pytest
import numpy as np
from glob import glob
from obspy import read, read_events, read_inventory, Stream
from obspy.io.sac.sactrace import SACTrace
from pysep.utils.cap_sac import (append_sac_headers,
                                 format_sac_header_w_taup_traveltimes)
from pysep.utils.curtail import (remove_for_clipped_amplitudes, rename_channels,
                                 remove_stations_for_missing_channels,
                                 remove_stations_for_insufficient_length,
                                 remove_traces_for_bad_data_types,
                                 subset_streams)
from pysep.utils.fmt import format_event_tag, format_event_tag_legacy
from pysep.utils.process import estimate_prefilter_corners
from pysep.utils.plot import plot_source_receiver_map
from pysep.utils.io import read_sem
from pysep.utils.mt import get_gcmt_moment_tensor, get_usgs_moment_tensor



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
def test_st():
    """
    Pre-saved event object
    """
    return read("./test_data/test_stream.mseed")


def test_append_sac_headers(test_st, test_inv, test_event):
    """
    Make sure we can write SAC headers correctly
    """
    st = append_sac_headers(st=test_st, inv=test_inv, event=test_event)
    assert(not hasattr(test_st[0].stats, "sac"))
    assert(hasattr(st[0].stats, "sac"))
    assert(st[0].stats.sac["evla"] == test_event.preferred_origin().latitude)


def test_event_tag_and_event_tag_legacy(test_event):
    """
    Check that event tagging works as expected
    """
    # while were here, make sure event tagging works
    tag = format_event_tag(test_event)
    assert(tag == "2009-04-07T201255_SOUTHERN_ALASKA")

    tag = format_event_tag_legacy(test_event)
    assert(tag == "20090407201255351")


def test_format_sac_headers_w_taup_traveltimes(test_st, test_inv, test_event):
    """
    Make sure we can write SAC headers correctly
    """
    st = append_sac_headers(st=test_st, inv=test_inv, event=test_event)
    st = format_sac_header_w_taup_traveltimes(st=st, model="ak135")
    assert(pytest.approx(st[0].stats.sac["a"], .01) == 224.535)


def test_sac_header_correct_origin_time(tmpdir, test_st, test_inv, test_event):
    """
    Make sure SAC headers are being written with 0 as the event origin time and
    not the trace start time
    """
    st = append_sac_headers(st=test_st, inv=test_inv, event=test_event)
    st[0].write(os.path.join(tmpdir, "test.sac"), format="SAC")  # only write 1
    sac = SACTrace.read(os.path.join(tmpdir, "test.sac"))
    assert(sac.reftime == test_event.preferred_origin().time)


def test_rename_channels(test_st):
    """
    Edit some waveforms to be obviously bad and make sure we can catch it

    TODO add function testing clipped station
    """
    st_edit = test_st.copy()
    # ATKA Location code incorrectly in station name
    st_edit[0].stats.channel = f"{st_edit[0].stats.channel}00"

    st = rename_channels(st_edit)
    assert(st[0].stats.location == "00")
    assert(st_edit[0].stats.location == "")


def test_remove_stations_for_insufficient_length(test_st):
    """
    Some of these data already have insufficient length
    """
    st_out = remove_stations_for_insufficient_length(test_st)
    # Removed 4 traces for insufficient length
    assert(len(test_st) - len(st_out) == 4)


def test_remove_stations_for_missing_channels(test_st):
    """
    Edit some stations to be obviously bad and make sure we can catch it
    """
    st_out = test_st.copy()
    # Dropping single components to test if the whole station drops out
    for tr in st_out.select(station="ATKA", component="Z"):
        st_out.remove(tr)
    for tr in st_out.select(station="ALPI", component="E"):
        st_out.remove(tr)
    st = remove_stations_for_missing_channels(st_out, networks="AK,YV")
    assert(len(st) == 3)
    assert(len(test_st) == 11)


def test_traces_for_bad_data_types(test_st):
    """
    Make sure data types that aren't float or integer arent' allowed through
    Sometimes LOG streams are stored as |S1 dtype which is string of len 1
    """
    st_out = test_st.copy()
    st_out[0].data = np.array(len(test_st[0].data) *["a"], dtype="|S1")
    st_out = remove_traces_for_bad_data_types(st_out)
    # Removed 1 trace for bad data type
    assert(len(test_st) - len(st_out) == 1)


def test_plot_map(test_event, test_inv):
    """
    Make a simple source-receiver map

    TODO This should probably have an image comparison to a baseline but for
        now we just make sure plotting works
    """
    plot_source_receiver_map(inv=test_inv, event=test_event, save=False,
                             show=False)


def test_read_sem():
    """
    Test reading SPECFEM-generated synthetics in as SAC formatted files
    """
    test_synthetics = glob("./test_data/test_synthetics/*.semd")
    test_stations = "./test_data/test_STATIONS"
    test_cmtsolution = "./test_data/test_CMTSOLUTION_2014p715167"

    st = Stream()
    for test_synthetic in test_synthetics:
        st += read_sem(fid=test_synthetic, source=test_cmtsolution,
                       stations=test_stations)

    assert(st[0].stats.sac.evla == -40.5405)


def test_estimate_prefilter_corners(test_st):
    """
    Test prefilter corner estimation based on trace starts time and samp rate
    """
    for tr in test_st:
        f0, f1, f2, f3 = estimate_prefilter_corners(tr)
        assert(f0 == pytest.approx(.0214, 3))
        assert(f1 == pytest.approx(.0427, 3))
        assert(f2 == pytest.approx(12.5, 3))
        assert(f3 == pytest.approx(25.0, 3))


@pytest.mark.skip(reason="need to find correct clipped example data")
def test_remove_for_clipped_amplitudes(test_st):
    """
    TODO
    """
    remove_for_clipped_amplitudes(test_st)


def test_subset_streams(test_st):
    """
    Make sure subsetting streams works by removing additional trace

    TODO this test will FAIL if data has multiple traces for a single component
        present in one stream but not another. This is due to how subset streams
        matches data (only by station ID). Users will need to merge data or take
        care that they do not have data gaps/ multiple traces.
    """
    n = len(test_st) // 2
    test_st_smaller = test_st[:n].copy()
    # test_st.pop(0)  # <- causes subset_streams to fail

    test_st_out, test_st_smaller_out = subset_streams(test_st, test_st_smaller)

    # subset streams should have reduced both streams to the shorter list
    assert(len(test_st_out) == n)
    assert(len(test_st_smaller_out) == n)

    # all station ids should be the same
    for tr_out, tr_smaller_out in zip(test_st_out, test_st_smaller_out):
        assert(tr_out.get_id() == tr_smaller_out.get_id())


def test_get_usgs_moment_tensor():
    """
    Just ensure that getting via USGS works as intended using an example event
    """
    event = test_get_gcmt_moment_tensor()
    del event.focal_mechanisms
    cat = get_usgs_moment_tensor(event=event)
    assert(len(cat) == 1)
    event = cat[0]
    assert hasattr(event, "focal_mechanisms")

    m_rr = event.preferred_focal_mechanism().moment_tensor.tensor.m_rr
    assert(pytest.approx(m_rr, 1E-20) == 4.81E20)


@pytest.mark.skip(reason="test already called by 'test_get_usgs_moment_tensor'")
def test_get_gcmt_moment_tensor():
    """
    Just ensure that getting via GCMT works as intended using an example event
    """
    # Kaikoura Earthquake
    origintime = "2016-11-13T11:02:00"
    magnitude = 7.8

    cat = get_gcmt_moment_tensor(event=None, origintime=origintime,
                                 magnitude=magnitude)
    assert(len(cat) == 1)
    event = cat[0]
    assert hasattr(event, "focal_mechanisms")
    m_rr = event.preferred_focal_mechanism().moment_tensor.tensor.m_rr
    assert(pytest.approx(m_rr, 1E-20) == 3.56E20)

    return event
