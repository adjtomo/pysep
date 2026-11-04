"""
Test the general functionality of the PrettyPlot waveform plotter (pplot).

.. note::
    Most tests only make figures to check that the machinery runs. Merging,
    response removal and the command line entry point are checked against
    expected values.
"""
import os
import sys
import pytest
import numpy as np
import matplotlib

matplotlib.use("Agg")

from obspy import read, read_inventory  # NOQA
from pysep import PrettyPlot  # NOQA
from pysep.pplot import fast_remove_response, merge_stream, main  # NOQA


TEST_DATA = os.path.join(os.path.dirname(__file__), "test_data")
MSEED = os.path.join(TEST_DATA, "test_stream.mseed")
INV = os.path.join(TEST_DATA, "test_inv.xml")
SYN = os.path.join(TEST_DATA, "test_synthetics", "NZ.BFZ.BXZ.semd")


@pytest.fixture
def fid_bhz(tmpdir):
    """
    Write a single gap-free channel to its own file

    :type tmpdir: py.path.local
    :param tmpdir: pytest temporary directory
    :rtype: str
    :return: path to the MiniSEED file holding AK.BESE..BHZ
    """
    fid = os.path.join(tmpdir, "AK.BESE..BHZ.mseed")
    read(MSEED).select(id="AK.BESE..BHZ").write(fid, format="MSEED")
    return fid


@pytest.fixture
def fids_bhz(tmpdir):
    """
    Write the three gap-free BHZ channels to one file each

    :type tmpdir: py.path.local
    :param tmpdir: pytest temporary directory
    :rtype: list of str
    :return: paths to the three MiniSEED files
    """
    fids = []
    for tr in read(MSEED).select(channel="BHZ"):
        fid = os.path.join(tmpdir, f"{tr.get_id()}.mseed")
        tr.write(fid, format="MSEED")
        fids.append(fid)
    return fids


def _run(fids, tmpdir, **kwargs):
    """
    Run PrettyPlot without showing the figure and check the saved figure

    :type fids: list of str
    :param fids: files to plot
    :type tmpdir: py.path.local
    :param tmpdir: pytest temporary directory to save the figure in
    :rtype: pysep.pplot.PrettyPlot
    :return: the PrettyPlot instance after running
    """
    save = os.path.join(tmpdir, "pplot.png")
    pp = PrettyPlot(fids, show=False, save=save, **kwargs)
    pp.run()
    assert os.path.exists(save)
    return pp


def test_default_plot(fid_bhz, tmpdir):
    """Plot raw data on a relative time axis, trimmed with `xlim`"""
    pp = _run([fid_bhz], tmpdir, xlim=["50", "150"])
    assert pp.ax.get_xlim() == (50., 150.)


def test_relative_time_units(fid_bhz, tmpdir):
    """`xlim` is given in the units of `time`, here minutes"""
    pp = _run([fid_bhz], tmpdir, time="m", xlim=["1", "2"])
    assert pp.ax.get_xlim() == (1., 2.)
    assert pp._xvals.min() < 1 < 2 < pp._xvals.max()


def test_response_removal_and_spectrogram(fid_bhz, tmpdir):
    """Remove the response, filter and add a spectrogram in local time"""
    pp = _run([fid_bhz], tmpdir, inv=INV, rr_pre_filt=[0.01, 0.02, 10, 20],
              fmin=0.1, fmax=5, time="a-09", spectrogram=True,
              xlim=["2009-04-07T11:13:00", "2009-04-07T11:15:00"])
    assert pp.ylabel == "Velocity [m/s]"
    assert pp.ax_spectra is not None
    # Response removal must undo the time zone shift afterwards
    assert pp.st[0].stats.starttime < pp.st[0].stats.endtime


def test_response_removal_inventory_object(fid_bhz, tmpdir):
    """`inv` also accepts an ObsPy Inventory"""
    pp = _run([fid_bhz], tmpdir, inv=read_inventory(INV), rr_output="DISP")
    assert pp.ylabel == "Displacement [m]"


def test_fast_remove_response_matches_obspy():
    """The interpolated response matches ObsPy's to a small tolerance"""
    st = read(MSEED).select(id="AK.BESE..BHZ")
    inv = read_inventory(INV)
    kwargs = dict(output="VEL", pre_filt=[0.01, 0.02, 10, 20])
    tr_obspy = st[0].copy().remove_response(inventory=inv, **kwargs)
    # Use a small `n_eval` so this short trace goes through the fast path
    tr_fast = fast_remove_response(st[0].copy(), inv, n_eval=2 ** 10,
                                   **kwargs)
    misfit = np.abs(tr_fast.data - tr_obspy.data).max()
    assert misfit / np.abs(tr_obspy.data).max() < 1E-4


def test_merge_gappy_file(tmpdir):
    """`merge` joins the segments of a gappy channel into one trace"""
    fid = os.path.join(tmpdir, "AK.ATKA..BHE.mseed")
    st = read(MSEED).select(id="AK.ATKA..BHE")
    assert len(st) == 2
    st.write(fid, format="MSEED")

    pp = PrettyPlot([fid], show=False)
    assert len(pp.st) == 1  # without merge only the first segment is kept
    assert pp.st[0].stats.npts == st[0].stats.npts

    pp = PrettyPlot([fid], merge=True, show=False)
    assert len(pp.st) == 1
    assert pp.st[0].stats.endtime == st[1].stats.endtime


def test_merge_stream_fill_values():
    """`merge_stream` accepts numbers and rejects unknown fill methods"""
    st = read(MSEED).select(id="AK.ATKA..BHE")
    merged = merge_stream(st, fill="0")
    assert len(merged) == 1
    assert not np.ma.isMaskedArray(merged[0].data)
    with pytest.raises(ValueError):
        merge_stream(st, fill="bad_value")


@pytest.mark.parametrize("wf_type", ["default", "stack", "recsec",
                                     "recsec_stack"])
def test_wf_types(fids_bhz, tmpdir, wf_type):
    """Overlain, stacked and record section plots of several files"""
    _run(fids_bhz, tmpdir, wf_type=wf_type, wf_norm=True, colors="C?",
         fmin=0.1, fmax=5)


def test_grid_plot(fids_bhz, tmpdir):
    """Grid plot with per-subplot y-limits"""
    ylim_grid = [-1, 1, -2, 2, -3, 3]
    pp = _run(fids_bhz, tmpdir, wf_type="grid", wf_norm=True,
              ylim_grid=ylim_grid)
    assert len(pp.axs) == 3
    for ax, (ymin, ymax) in zip(pp.axs, np.reshape(ylim_grid, (3, 2))):
        assert ax.get_ylim() == (ymin, ymax)


def test_grid_with_spectrogram_fails(fid_bhz):
    """Spectrograms cannot be combined with grid plots"""
    with pytest.raises(AssertionError):
        PrettyPlot([fid_bhz], wf_type="grid", spectrogram=True)


def test_colors_colormap(fids_bhz, tmpdir):
    """A colormap name is sampled into one color per trace"""
    pp = _run(fids_bhz, tmpdir, colors="viridis")
    assert len(pp._colors) == 3


def test_tmarks_taup_group_vels(fid_bhz, tmpdir):
    """Time marks, TauP arrivals and group velocities on relative time"""
    _run([fid_bhz], tmpdir, fmin=0.1, fmax=5, tmarks=["60", "90"],
         tmarks_c=["r", "b"], tp_phases=["P", "S"], tp_dist_km=300,
         tp_depth=10, tp_start="20", group_vels=[3.5, 3.0])


def test_add_trace(fid_bhz, tmpdir):
    """Overlay a second time series on a twin y-axis"""
    pp = _run([fid_bhz], tmpdir, time="a", add_trace=MSEED,
              tr_label="other", tr_ylabel="Counts")
    assert pp.twax is not None


def test_read_specfem_synthetic(tmpdir):
    """SPECFEM two-column ASCII is read with PySEP's `read_sem`"""
    pp = _run([SYN], tmpdir)
    assert pp.st[0].stats.channel == "BXZ"


def test_command_line(fid_bhz, tmpdir, monkeypatch):
    """The `pplot` entry point turns the negative flags into keywords"""
    save = os.path.join(tmpdir, "pplot_cli.png")
    monkeypatch.setattr(sys, "argv", [
        "pplot", fid_bhz, "--inv", INV, "--fmin", "0.1", "--noshow",
        "--no_legend", "--sp_dbscale_off", "--spectrogram", "--save", save
    ])
    main()
    assert os.path.exists(save)
