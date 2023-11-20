"""
Test the general functionality of Pysep
"""
from re import L
import pytest
import os
import yaml
from glob import glob
from pysep import Pysep, logger


@pytest.fixture
def config_file():
    """
    Test config file for a 2009 Anchorage earthquake
    """
    return "./test_data/test_config_anchorage.yaml"


@pytest.fixture
def test_pysep(config_file):
    """
    Test Pysep instance based on a Config file, gather event and StationXML
    """
    pysep = Pysep(config_file=config_file, log_level=None)
    pysep.load()
    pysep.check()
    pysep.c = pysep.get_client()
    pysep.event = pysep.get_event()
    pysep.inv = pysep.get_stations()

    return pysep


def test_get_client(test_pysep):
    """
    Simple check that we're looking at IRIS for our client
    """
    assert(test_pysep.c.base_url == "http://service.iris.edu")


def test_get_event(test_pysep):
    """
    Simple check that the event created by User parameters is right
    """
    assert(test_pysep.event.preferred_magnitude().mag == 4.6)


def test_get_inventory(test_pysep):
    """
    Check that the collected stations are correct
    :param test_pysep:
    :return:
    """
    inv = test_pysep.inv
    assert(len(inv.get_contents()["channels"]) == 9)
    assert(f"{inv[0].code}{inv[0][0].code}{inv[0][0][0].code}" == "AKATKABHE")


def test_get_waveform(test_pysep):
    """
    Get waveforms for a few test stations from IRIS
    :param test_pysep:
    :return:
    """
    st = test_pysep.get_waveforms()
    assert(len(st) == 11)
    assert(st[-1].get_id() == "YV.ALPI..BHZ")


def test_curtail_stations(test_pysep):
    """
    Exclude a station based on maximum distance getting exceeded
    """
    inv = test_pysep.curtail_stations()

    nsta_pre_curtail = len(test_pysep.inv.get_contents()["channels"])
    nsta_post_curtail = len(inv.get_contents()["channels"])

    assert(nsta_pre_curtail - nsta_post_curtail == 6)


def test_config_files_have_correct_parameters(tmpdir, test_pysep):
    """
    Check that all the config files have the correct parameters based on what
    PySEP writes to ensure that the config files are up to date with codebase
    """
    # Write a config file from the Pysep instance and read in the yaml file
    fid = os.path.join(tmpdir, "test_config.yaml")
    test_pysep.write_config(fid=fid)
    config = set(yaml.safe_load(open(fid)).keys())

    # Get the directory containing the example config files
    config_dir = "../configs/"

    # Iterate over all the yaml files in the directory
    for dir_ in glob(os.path.join(config_dir, "*")):
        if not os.path.isdir(dir_):
            continue
        for fid in glob(os.path.join(dir_, "*")):
            if not fid.endswith(".yaml"):
                continue
            with open(fid, "r") as f:
                config_check = set(yaml.safe_load(f).keys())
                assert(config.issubset(config_check)), \
                    f"Config file {fid} does not have the correct parameters"