"""
Test the declustering algorithm
"""
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


def test_decluster_events(test_declust):
    """
    Test event declustering algorithm
    """
    for select_by in ["magnitude", "magnitude_r", "depth", "depth_r"]:
        cat = test_declust.decluster_events(nx=2, ny=2, min_mags=[4.5],
                                            nkeep=4, select_by=select_by)
        assert(len(cat) == 14)


def test_plot(test_declust):
    """

    """
    test_declust.plot(inv=True)