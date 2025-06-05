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


# For debugging, to show figures turn SHOW=True
SHOW=False

@pytest.fixture
def recsec(tmpdir):
    """Initiate a RecordSection instance"""
    return RecordSection(pysep_path="./test_data/test_SAC", 
                         scale_by="normalize", show=SHOW,
                         save=os.path.join(tmpdir, "recsec.png"))


@pytest.fixture
def recsec_w_synthetics(tmpdir):
    """Initiate a RecordSection instance"""
    return RecordSection(pysep_path="./test_data/test_SAC",
                         syn_path="./test_data/test_synthetics",
                         source="./test_data/test_CMTSOLUTION_2014p715167",
                         stations="./test_data/test_STATIONS",
                         scale_by="normalize",  # for visual checks
                         show=SHOW, save=os.path.join(tmpdir, "recsec.png"))


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


def test_plot_recsec_time_shift_array(recsec):
    """apply an array of time shifts to shift each trace differently"""
    recsec.time_shift_s = [50, 125, 200]
    recsec.zero_pad_s = [200, 500]
    recsec.move_out = 4
    recsec.process_st()
    recsec.get_parameters()
    recsec.plot()

def test_plot_recsec_time_shift_syn_same_as_obs(tmpdir):
    """apply the same time shift to data and synthetics"""
    # Cannot use fixtures because the time shift values need to be set at init
    # for this specific use case. This is just a caveat of the test case, Users
    # will not run into this issue
    recsec_w_synthetics = RecordSection(
        pysep_path="./test_data/test_SAC", 
        syn_path="./test_data/test_synthetics",
        source="./test_data/test_CMTSOLUTION_2014p715167",
        stations="./test_data/test_STATIONS",
        scale_by="normalize",  # for visual checks
        time_shift_s=88., zero_pad_s=[200, 500],
        show=SHOW, save=os.path.join(tmpdir, "recsec.png")
        )

    recsec_w_synthetics.process_st()
    recsec_w_synthetics.get_parameters()
    recsec_w_synthetics.plot()


def test_plot_recsec_time_shift_syn_array(recsec_w_synthetics):
    """apply an array of time shifts to shift each trace differently"""
    recsec_w_synthetics.time_shift_s = [50, 125, -200]
    recsec_w_synthetics.time_shift_s_syn = [-200, -100, 25]
    recsec_w_synthetics.zero_pad_s = [200, 500]

    recsec_w_synthetics.process_st()
    recsec_w_synthetics.get_parameters()
    recsec_w_synthetics.plot()


def test_plot_recsec_time_shift_syn(recsec_w_synthetics):
    """apply different time shift to data and synthetics"""
    recsec_w_synthetics.time_shift_s = 88.
    recsec_w_synthetics.time_shift_s_syn = -51.
    recsec_w_synthetics.zero_pad_s = [200, 500]

    recsec_w_synthetics.process_st()
    recsec_w_synthetics.get_parameters()
    recsec_w_synthetics.plot()


def test_plot_recsec_time_shift_syn_zero(recsec_w_synthetics):
    """make sure zero synthetic time shift works"""
    recsec_w_synthetics.time_shift_s = 88.
    recsec_w_synthetics.time_shift_s_syn = 0.
    recsec_w_synthetics.zero_pad_s = [200, 500]

    recsec_w_synthetics.process_st()
    recsec_w_synthetics.get_parameters()
    recsec_w_synthetics.plot()


def test_plot_recsec_time_shift_zero_w_syn(recsec_w_synthetics):
    """make sure zero time shift and nonzero synthetic time shift works"""
    recsec_w_synthetics.time_shift_s = 0.
    recsec_w_synthetics.time_shift_s_syn = 69.
    recsec_w_synthetics.zero_pad_s = [200, 500]

    recsec_w_synthetics.process_st()
    recsec_w_synthetics.get_parameters()
    recsec_w_synthetics.plot()


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


def test_recsec_calc_time_offset_no_trim(recsec_w_synthetics):
    """testing that synthetics and data which do not share origin time 
    plot together correctly by checking that the time offsets are calced"""
    # Turn off trim so that time shifts are retained
    recsec_w_synthetics.trim = False

    # Pad 100s zeros to data and shift starttime to match
    for tr in recsec_w_synthetics.st:
        tr.data = np.append(np.zeros(int(100 * tr.stats.sampling_rate)), 
                            tr.data)
        tr.stats.starttime -= 100

    recsec_w_synthetics.process_st()
    recsec_w_synthetics.get_parameters()

    for tr in recsec_w_synthetics.st:
        assert(tr.stats.time_offset == -100)


def test_recsec_calc_time_offset_w_trim(recsec_w_synthetics):
    """testing that synthetics and data which do not share origin time 
    plot together correctly by checking that trim is applied"""
    # Pad 100s zeros to data and shift starttime to match
    for tr in recsec_w_synthetics.st:
        tr.data = np.append(np.zeros(int(100 * tr.stats.sampling_rate)), 
                            tr.data)
        tr.stats.starttime -= 100

    recsec_w_synthetics.process_st()
    recsec_w_synthetics.get_parameters()

    for tr, tr_syn in zip(recsec_w_synthetics.st, recsec_w_synthetics.st_syn):
        tdiff = tr.stats.starttime - tr_syn.stats.starttime
        assert(tdiff < tr.stats.delta)


def test_recsec_zero_amplitude(recsec):
    """
    waveforms that have zero amplitude and are normalized should be able 
    to bypass normalizations which lead to weird plotting (see #131).

    .. note::

        This does not really test that the method is working correctly because
        dividing a NumPy array by zero leads to NaNs in the array which just
        won't plot. This is more of a visual test to make sure that the 
        zero amplitude is plotting correctly, look for green lines
    """
    recsec.kwargs.scale_by = "normalize"
    recsec.kwargs.obs_color = "green"
    recsec.linewidth = 30
    for tr in recsec.st:
        tr.data *= 0
    recsec.process_st()
    recsec.get_parameters()
    recsec.plot()


def test_recsec_mismatched_locations_fails(recsec_w_synthetics):
    """
    Test that when RecSec encounters stations with different location codes, 
    it fails the parameter check
    """
    # Ensure that the location codes for `st` and `st_syn` differ 
    for tr in recsec_w_synthetics.st:
        tr.stats.location = "11"
    for tr in recsec_w_synthetics.st_syn:
        tr.stats.location = "00"

    with pytest.raises(SystemExit):
        recsec_w_synthetics.check_parameters()

        
def test_recsec_mismatched_locations(recsec_w_synthetics):
    """
    Test when RecSec encounters stations that match all their station code 
    except for location code, which is sometimes used as a free-floating 
    parameter to indicate e.g., synthetics or something and doesn't always
    determine that waveforms do not match.
    """
    # Ensure that the location codes for `st` and `st_syn` differ 
    for tr in recsec_w_synthetics.st:
        tr.stats.location = "11"
    for tr in recsec_w_synthetics.st_syn:
        tr.stats.location = "00"

    recsec_w_synthetics.remove_locations = True
    recsec_w_synthetics.check_parameters()
    recsec_w_synthetics.process_st()
    recsec_w_synthetics.get_parameters()
    recsec_w_synthetics.plot()
    
    
@pytest.mark.skip(reason="visual test for zero pad acceptance, unskip to view")
def test_plot_recsec_zero_pad(recsec):
    """
    Ensure that zero pad only adds zeros to the start and end of the waveform
    but does NOT time shift it. 
    
    .. note::
        This is a visual check for now so unskip this test and set constant 
        SHOW==True to have a look at the waveforms that get plotted. In the 
        future we could do a check on the time of the max amplitude but I  
        did not think it was worth the investment to write that test.
    """
    recsec.scale_by = "normalize"
    recsec.min_period_s = 1/3
    recsec.max_period_s = 10.
    recsec.xlim_s = [-25, 50]
    recsec.tmarks = [20]

    recsec_copy = copy(recsec)

    # Zero pad, preprocess, and find the index of the maximum amplitude
    recsec.process_st()
    recsec.get_parameters()
    recsec.plot()

    # Same processing but add zero pad
    recsec_zeropad = copy(recsec_copy)
    recsec_zeropad.zero_pad_s = [10, 5]
    recsec_zeropad.process_st()
    recsec_zeropad.get_parameters()
    recsec_zeropad.plot()
