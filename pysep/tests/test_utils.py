"""
Test the general utility functions of Pysep which includes quality
checking waveforms and formatting data with metadata

TODO add function testing clipped station, need to find example data
"""
import os
import pytest
from obspy import read, read_events, read_inventory

from pysep.utils.cap_sac import (append_sac_headers,
                                 format_sac_header_w_taup_traveltimes)
from pysep.utils.curtail import (remove_for_clipped_amplitudes, rename_channels,
                                 remove_stations_for_missing_channels,
                                 remove_stations_for_insufficient_length)
from pysep.utils.fmt import format_event_tag
from pysep.utils.plot import plot_source_receiver_map


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

    # while were here, make sure event tagging works
    tag = format_event_tag(test_event)
    assert(tag == "2009-04-07T201255_SOUTHERN_ALASKA")


def test_format_sac_headers_w_taup_traveltimes(test_st, test_inv, test_event):
    """
    Make sure we can write SAC headers correctly
    """
    st = append_sac_headers(st=test_st, inv=test_inv, event=test_event)
    st = format_sac_header_w_taup_traveltimes(st=st, model="ak135")
    assert(pytest.approx(st[0].stats.sac["user4"], .01) == 57.66)


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


def test_plot_map(test_event, test_inv):
    """
    Make a simple source-receiver map

    TODO This should probably have an image comparison to a baseline but for
        now we just make sure plotting works
    """
    plot_source_receiver_map(inv=test_inv, event=test_event, save=False,
                             show=False)


@pytest.mark.skip(reason="need to find correct clipped example data")
def test_remove_for_clipped_amplitudes(test_st):
    """
    TODO
    """
    remove_for_clipped_amplitudes(test_st)

