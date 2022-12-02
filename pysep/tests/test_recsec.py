"""
Test the general functionality of the RECord SECtion plotter.

.. note::
    This test suite simply creates some plots to make sure the general
    machinery is working as expected, but otherwise has no real solid checks on
    the outputs. If this code grows more cumbersome, we will want to do figure
    matching checks or something more sophisticated.
"""
import os
import pytest
import numpy as np
from copy import copy
from pysep import RecordSection


@pytest.fixture
def recsec(tmpdir):
    """Initiate a RecordSection instance"""
    return RecordSection(pysep_path="./test_data/test_SAC", show=False,
                         save=os.path.join(tmpdir, "recsec.png"))


@pytest.fixture
def recsec_w_synthetics(tmpdir):
    """Initiate a RecordSection instance"""
    return RecordSection(pysep_path="./test_data/test_SAC",
                         syn_path="./test_data/test_synthetics",
                         cmtsolution="./test_data/test_CMTSOLUTION_2014p715167",
                         stations="./test_data/test_STATIONS",
                         show=False, save=os.path.join(tmpdir, "recsec.png"))


def test_plot_recsec(recsec):
    """Simply test out the functinoality of plotw_rs"""
    recsec.process_st()
    recsec.get_parameters()
    recsec.plot()


def test_plot_recsec_w_synthetics(recsec_w_synthetics):
    """Simply test out the functinoality of plotw_rs"""
    recsec_w_synthetics.process_st()
    recsec_w_synthetics.get_parameters()
    recsec_w_synthetics.plot()


def test_plot_recsec_sort_by(recsec):
    """Simply test out the functinoality of plotw_rs"""
    acceptable_sort_by = ["default", "azimuth", "backazimuth",
                          "distance", "alphabetical", "abs_azimuth",
                          "abs_distance"]
    for sort_by in acceptable_sort_by:
        recsec_test = copy(recsec)
        recsec_test.sort_by = sort_by
        recsec_test.process_st()
        recsec_test.get_parameters()
        recsec_test.save = f"{recsec.save}_{sort_by}.png"
        recsec_test.plot()


def test_plot_recsec_scale_by(recsec):
    """scale by option testing"""
    acceptable_scale_by = ["normalize", "global_norm"]
                           # "geometric_spreading"]
    for scale_by in acceptable_scale_by:
        recsec_test = copy(recsec)
        recsec_test.scale_by = scale_by
        recsec_test.process_st()
        recsec_test.get_parameters()
        recsec_test.save = f"{recsec.save}_{scale_by}.png"
        recsec_test.plot()


def test_plot_recsec_time_shift(recsec):
    """apply a whole bunch of time shift elements"""
    recsec.time_shift_s = 100
    recsec.zero_pad_s = [200, 500]
    recsec.move_out = 4
    recsec.process_st()
    recsec.get_parameters()
    recsec.plot()


def test_plot_recsec_preprocess(recsec):
    """preprocess and filter the record section"""
    recsec.min_period_s = 10
    recsec.max_period_s = 30
    recsec.integrate = 2
    recsec.components = "ZN"
    recsec.process_st()
    recsec.get_parameters()
    recsec.plot()


def test_plot_recsec_y_label_loc(recsec):
    """test all available y_label locations options"""
    acceptable_y_label_loc = ["default", "y_axis", "y_axis_right", "x_min",
                              "x_max", None]
    recsec.process_st()
    recsec.get_parameters()
    for y_label_loc in acceptable_y_label_loc:
        recsec.y_label_loc = y_label_loc
        recsec.plot()


def test_plot_recsec_distance_units(recsec):
    """test all available distance units"""
    acceptable_distance_units = ["km", "km_utm", "deg"]
    for distance_units in acceptable_distance_units:
        recsec.distance_units = distance_units
        recsec.process_st()
        recsec.get_parameters()
        recsec.plot()


def test_plot_recsec_plot_options(recsec):
    """test plotting options"""
    recsec.y_axis_spacing = 3.5
    recsec.amplitude_scale_factor = 8.2
    recsec.azimuth_start_deg = 102.1
    recsec.process_st()
    recsec.get_parameters()
    recsec.plot()


def test_recsec_calc_time_offset(recsec_w_synthetics):
    """testing that synthetics and data which do not share origin time 
    plot together correctly by checking that the time offsets are calced"""
    # Pad 100s zeros to data and shift starttime to match
    for tr in recsec_w_synthetics.st:
        tr.data = np.append(np.zeros(int(100 * tr.stats.sampling_rate)), 
                            tr.data)
        tr.stats.starttime -= 100

    recsec_w_synthetics.process_st()
    recsec_w_synthetics.get_parameters()
    for tr in recsec_w_synthetics.st:
        assert(tr.stats.time_offset == -100)
