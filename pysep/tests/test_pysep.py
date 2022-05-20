"""
Test the general functionality of Pysep
"""
import pytest
from pysep import Pysep


@pytest.fixture
def config_file():
    """
    Test config file
    :return:
    """
    return "./test_data/test_config_anchorage.yaml"


@pytest.fixture
def test_pysep(config_file):
    """
    Test Pysep instance
    :param config_file:
    :return:
    """
    pysep = Pysep(config_file=config_file)
    pysep.load_config()
    pysep.check()
    return pysep


def test_get_client(test_pysep):
    """

    :param test_pysep:
    :return:
    """
    client = test_pysep.get_client()
    assert(client.base_url == "http://service.iris.edu")


def test_get_event(test_pysep):
    """

    :param test_pysep:
    :return:
    """
    event = test_pysep.get_event()
    assert(event.preferred_magnitude().mag == 4.6)


def test_get_inventory(test_pysep):
    """

    :param test_pysep:
    :return:
    """
    test_pysep._client = test_pysep.get_client()
    inv = test_pysep.get_stations()
    assert(f"{inv[0].code}{inv[0][0].code}{inv[0][0][0].code}" == "AKATKABHE")


def test_get_waveform(test_pysep):
    """

    :param test_pysep:
    :return:
    """
    test_pysep._client = test_pysep.get_client()
    test_pysep.inv = test_pysep.get_stations()
    st = test_pysep.get_waveforms()
    assert(len(st) == 11)
    assert(st[-1].get_id() == "YV.ALPI..BHZ")
