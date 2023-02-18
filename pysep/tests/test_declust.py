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

def test_data_availability(test_declust):
    """
    Make sure data availability calculation returns the same each time
    """
    # Just check one of the events in the catalog
    resource_id = "smi:service.iris.edu/fdsnws/event/1/query?eventid=4599111"
    assert(len(test_declust.data_avail[resource_id]) == 15)


def test_threshold_catalog(test_declust):
    """
    Ensure that thresholding removes the proper number of events
    """
    test_declust.threshold_catalog(zedges=[0, 8, 36, 100], min_data=5,
                                   min_mags=[3, 3, 5])
    assert(len(test_declust.cat) == 215)


def test_decluster_events_cartesian(test_declust):
    """
    Test event declustering algorithm
    """
    for select_by in ["magnitude", "magnitude_r", "depth", "depth_r",
                      "data", "data_r"]:
        cat = test_declust.decluster_events(choice="cartesian", nx=2, ny=2,
                                            min_mags=[4.5], nkeep=4,
                                            select_by=select_by)
        assert(len(cat) == 16)


def test_decluster_events_polar(test_declust):
    """
    Test event declustering algorithm
    """
    for select_by in ["magnitude", "magnitude_r", "depth", "depth_r",
                      "data", "data_r"]:
        cat = test_declust.decluster_events(choice="polar", nx=2, ny=2,
                                            min_mags=[4.5], nkeep=4,
                                            select_by=select_by)
        assert(len(cat) == 24)


def test_decluster_plot_cartesian(tmpdir, test_declust):
    """
    Test event declustering algorithm
    """
    cat = test_declust.decluster_events(
        choice="cartesian", nx=25, ny=25, zedges=[0, 10, 35],
        min_mags=[4., 5.], nkeep=[4, 4], select_by="magnitude", plot=True,
        plot_dir=tmpdir
    )
    assert(os.path.exists(os.path.join(tmpdir, "pre_decluster_crtsn.png")))
    assert(os.path.exists(os.path.join(tmpdir, "declustered_crtsn.png")))
    assert(os.path.exists(os.path.join(tmpdir, "removed_crtsn.png")))


def test_decluster_plot_polar(tmpdir, test_declust):
    """
    Test event declustering algorithm
    """
    cat = test_declust.decluster_events(
        choice="polar", zedges=[0, 35], min_mags=None, nkeep=5,
        select_by="magnitude_r", plot=True, plot_dir=tmpdir
    )
    assert(os.path.exists(os.path.join(tmpdir, "pre_decluster_plr.png")))
    assert(os.path.exists(os.path.join(tmpdir, "decluster_plr.png")))
    assert(os.path.exists(os.path.join(tmpdir, "removed_plr.png")))


def test_plot(tmpdir, test_declust):
    """
    Test the plot function on its own
    """
    test_declust.plot(inv=True, show=False, color_by="data", cmap="inferno",
                      save=os.path.join(tmpdir, "test_plot.png"))
    assert(os.path.exists(os.path.join(tmpdir, "test_plot.png")))


def test_srcrcv_weight(tmpdir, test_declust):
    """test srcrv weight calculation"""
    test_declust.calculate_srcrcv_weights(
        write=os.path.join(tmpdir, "weights.txt"), plot=True, show=False,
        save=os.path.join(tmpdir, "srcrcvwght.png")
    )

