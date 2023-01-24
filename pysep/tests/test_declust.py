"""
Test the declustering algorithm
"""
import os
import pytest
from obspy import read_events, read_inventory
from pysep import Declust

from pysep import logger
logger.setLevel("DEBUG")


@pytest.fixture
def test_declust():
    """
    Setup a decluster class
    """
    cat = read_events("./test_data/test_nalaska_events.xml")
    inv = read_inventory("./test_data/test_nalaska_inv.xml")

    return Declust(cat=cat, inv=inv)


def test_decluster_events_cartesian(test_declust):
    """
    Test event declustering algorithm
    """
    for select_by in ["magnitude", "magnitude_r", "depth", "depth_r"]:
        cat = test_declust.decluster_events(choice="cartesian", nx=2, ny=2,
                                            min_mags=[4.5], nkeep=4,
                                            select_by=select_by)
        assert(len(cat) == 14)


def test_decluster_events_polar(test_declust):
    """
    Test event declustering algorithm
    """
    for select_by in ["magnitude", "magnitude_r", "depth", "depth_r"]:
        cat = test_declust.decluster_events(choice="polar", nx=2, ny=2,
                                            min_mags=[4.5], nkeep=4,
                                            select_by=select_by)
        assert(len(cat) == 20)


def test_decluster_plot_cartesian(tmpdir, test_declust):
    """
    Test event declustering algorithm
    """
    cat = test_declust.decluster_events(
        choice="cartesian", nx=25, ny=25, zedges=[0, 10, 35],
        min_mags=[4., 5.], nkeep=[4, 4], select_by="magnitude", plot=True,
        plot_dir=tmpdir
    )
    assert(os.path.exists(os.path.join(tmpdir, "decluster_event_catalog.png")))
    assert(os.path.exists(os.path.join(tmpdir, "original_event_catalog.png")))


def test_decluster_plot_polar(tmpdir, test_declust):
    """
    Test event declustering algorithm
    """
    cat = test_declust.decluster_events(
        choice="polar", zedges=[0, 35], min_mags=None, nkeep=5,
        select_by="magnitude_r", plot=True, plot_dir=tmpdir
    )
    assert(os.path.exists(os.path.join(tmpdir, "decluster_event_catalog.png")))
    assert(os.path.exists(os.path.join(tmpdir, "original_event_catalog.png")))


def test_plot(tmpdir, test_declust):
    """

    """
    test_declust.plot(inv=True, show=False,
                      save=os.path.join(tmpdir, "test_plot.png"))
    assert(os.path.exists(os.path.join(tmpdir, "test_plot.png")))