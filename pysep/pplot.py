"""
PrettyPlot (pplot): pretty plotting of ObsPy Streams with options for
spectrograms, TauP arrivals, plotting other time-series data, record sections,
grid plotting, instrument response removal and processing.

Run from the command line with ``pplot -h``, or from Python with
``from pysep.pplot import PrettyPlot``.
"""
import argparse
import sys
import os
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.fft import next_fast_len
from scipy.ndimage import uniform_filter1d
from scipy.signal import hilbert

from dateutil.rrule import MINUTELY, SECONDLY
from matplotlib import mlab
from matplotlib.colors import Normalize
from matplotlib.dates import date2num, AutoDateLocator
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable

from obspy import read, read_inventory, UTCDateTime, Stream, Inventory
from obspy.taup import TauPyModel
from obspy.geodetics import kilometers2degrees, degrees2kilometers
from obspy.imaging.cm import obspy_sequential
from obspy.imaging.util import ObsPyAutoDateFormatter

# `logger` is defined in pysep/__init__.py before its submodule imports, so
# this is safe. `read_sem` comes from the submodule, not the `pysep`
# namespace, so that importing this module from pysep/__init__.py does not
# create a circular import
from pysep import logger
from pysep.utils.io import read_sem

SECONDS_PER_DAY = 3600.0 * 24.0

# Central control for all figure fontsizes. Change a value here to update it
# everywhere it is used in the script.
FONTSIZES = {
    "title": 14,
    "subplot_label": 27.5,
    "axis_label": 12.,
    "xtick_label": 12.,
    "ytick_label": 12.,
    "legend": 10,
    "yaxis_offset": 12,
    "taup_phase_label": 12,
    "group_velocity_label": 12,
    "spectrogram_colorbar_label": 10,
    "spectrogram_trace_id": 8,
    "xtick_date_label": "small",
    "waveform_index_label": "large",
    "waveform_index_label_highlight": "large",
}

def parse_args():
    """All modifications are accomplished with command line arguments"""
    parser = argparse.ArgumentParser()

    # Waveform Processing
    parser.add_argument("fids", nargs="+", 
                        help="required, file ID(s), wildcards accepted")
    parser.add_argument("--merge", action="store_true",
                        help="merge traces with the same ID into one, "
                             "e.g., a gappy file or consecutive day files "
                             "of one channel. Gaps are filled according to "
                             "`merge_fill`. Every channel found is plotted "
                             "and `idx` is ignored")
    parser.add_argument("--merge_fill", type=str, default="interpolate",
                        help="how `merge` fills gaps: 'interpolate' "
                             "(default, linear between the samples on "
                             "either side), 'latest' (repeat the last "
                             "sample) or a number, e.g. 0")
    parser.add_argument("--idx", nargs="?", type=int, default=0,
                        help="if the file in `fids` read has more than one " 
                             "trace in it, allows user to choose index of " 
                             "which trace to plot. May have unintended " 
                             "consequences if you have multiple `fids`")
    parser.add_argument("-tp", "--taper", nargs="?", type=float, default=0,
                        help="optional taper percentange")
    parser.add_argument("-f1", "--fmin", nargs="?", type=float, default=None,
                        help="optional filtering freqmin in Hz")
    parser.add_argument("-f2", "--fmax", nargs="?", type=float, default=None,
                        help="optional filtering freqmax in Hz")
    parser.add_argument("-z", "--zerophase", action="store_true", default=False,
                        help="apply zerophase filter or not")
    parser.add_argument("-C", "--corners", type=int, default=4,
                        help="number of corners to apply on ")
    parser.add_argument("-r", "--resample", type=float, default=None,
                        help="resample the data prior to processing")
    parser.add_argument("-t0", "--t0", nargs="?", type=float, default=0,
                        help="SPECFEM USER_T0 if synthetics")
    parser.add_argument("-ts", "--tstart", nargs="?", type=float, default=0,
                        help="for relative time axis, set t0 value, defaults 0")
    parser.add_argument("--detrend", nargs="?", type=str, default=None,
                        help="demean the data prior to other processing",
                        choices=["simple", "linear", "demean"])
    parser.add_argument("--integrate", nargs="?", type=int, default=0,
                        help="Integrate the time series, will demean and taper,"
                             " value for integrate will be the number of times")
    parser.add_argument("--differentiate", nargs="?", type=int, default=0,
                        help="Diff. the time series, will demean and taper,"
                             " value for integrate will be the number of times")
    parser.add_argument("--trim_pct", nargs="?", type=float, default=1.0,
                        help="Buffer amount in time to trim the data if using " 
                             "`xlim` to prevent processing the tails of the " 
                             "cut data, given in percentage. That is if the " 
                             "requested cut length is 30s, an additional 30s " 
                             "on either side of the trace will be cut and " 
                             "processed")
    
    # Instrument response removal (all parameters start with rr_), these
    # map onto ObsPy's Stream.remove_response()
    parser.add_argument("--inv", nargs="+", type=str, default=None,
                        help="StationXML file(s) used to remove the "
                             "instrument response from the raw data. "
                             "Response removal happens after trimming and "
                             "before all other processing. See all 'rr_*' "
                             "parameters to control it")
    parser.add_argument("--rr_output", type=str, default="VEL",
                        choices=["DISP", "VEL", "ACC", "DEF"],
                        help="output units after response removal: "
                             "displacement, velocity, acceleration or "
                             "'DEF' for the response's default units")
    parser.add_argument("--rr_water_level", type=_float_or_none, default=60.,
                        help="water level in dB for the deconvolution, "
                             "'none' to turn it off. Defaults to 60")
    parser.add_argument("--rr_pre_filt", nargs=4, type=float, default=None,
                        help="four corner frequencies (f1 f2 f3 f4) in Hz "
                             "of a cosine taper applied in the frequency "
                             "domain before deconvolution")
    parser.add_argument("--rr_zero_mean_off", action="store_true",
                        help="do not demean the data before deconvolution")
    parser.add_argument("--rr_taper_off", action="store_true",
                        help="do not taper the data before deconvolution")
    parser.add_argument("--rr_taper_fraction", type=float, default=0.05,
                        help="taper fraction of the cosine taper applied "
                             "before deconvolution. Defaults to 0.05")
    parser.add_argument("--rr_exact", action="store_true",
                        help="use ObsPy's remove_response() as is, which "
                             "evaluates the response at every FFT "
                             "frequency. Slower for long windows; by "
                             "default the response is evaluated on a "
                             "coarse grid and interpolated")

    # Waveform options
    parser.add_argument("--wf_type", nargs="?", type=str, default="default",
                        help="Option for how to plot the waveforms:\n"
                             "-'default': plot all `fids` on top of each "
                             "other with absolute amplitudes.\n"
                             "-'stack': stack all `fids` on top of each "
                             "other w some transparency, and then plot the "
                             "mean value of all traces as a dark line on top "
                             "of all\n"
                             "-'recsec': plot all waveforms with relative "
                             "amplitudes on a record section style plot so "
                             "that each waveform is visually distinct from the "
                             "others\n"
                             "-'recsec_stack', recsec but add stack\n"
                             "-'grid': plot subplots for each waveform in " 
                             "`fids`. cant be used with spectrogram\n",
                        choices=["default", "stack", "recsec", "recsec_stack",
                                 "grid"])
    parser.add_argument("--wf_norm", action="store_true", 
                        help="normalize traces to their individual max")
    parser.add_argument("--wf_scale", type=float, default=1,
                        help="scale the amplitudes by an arbitrary value, used "
                             "to emphasize small arrivals")
    parser.add_argument("--wf_abs", action="store_true", 
                        help="plot the absolute value of the waveform, does" 
                             "not affect processing, just final plotting")
    parser.add_argument("--wf_env", action="store_true",
                        help="overlay the envelope of each waveform (the "
                             "magnitude of its analytic signal from the "
                             "Hilbert transform) as a solid line")
    parser.add_argument("--wf_env_smooth", type=float, default=0.,
                        help="length in seconds of a centred running mean "
                             "applied to the envelope, e.g. 60 to compare "
                             "with 1-min mean data. 0 (default) for none")
    parser.add_argument("--wf_env_color", type=str, default=None,
                        help="color of the envelope line(s). Defaults to "
                             "'C0' for a black waveform, otherwise the "
                             "waveform color")
    parser.add_argument("--wf_recsec_spacing", nargs="?", type=float, default=1,
                        help="if `wf_type` is some form of recsec, choose the " 
                             "spacing modifier between each of the waveforms. " 
                             "By default it matches the order of magnitude so " 
                             "this modifier set to 1 is a good starting guess "
                             "but you can move it up or down by steps of 0.1 " 
                             "to get waveforms closer or further")
    parser.add_argument("--wf_spacing_exact", action="store_true", 
                        help="For `wf_recsec_spacing` the value is used exact "
                             "rather than as an order of magnitude scale")
    parser.add_argument("--wf_order", type=int, nargs="+", default=None,
                        help="allow custom order when plotting waveforms, " \
                             "based on the alphabetical index. First you need "
                             "to plot the normal sorted one and then you can" \
                             "rearrange based on that order")
    parser.add_argument("--wf_grid", nargs=2, type=int, default=None,
                        help="(nrow, ncol) If `wf_type` is 'grid', determine " \
                            "the number of rows and columns of that grid. " \
                            "For example, '--wf_grid 2 3' will give you 2 rows "
                            "3 columns. Default is one column")
    parser.add_argument("--nrows", type=int, default=None,
                        help="Number of rows for subplots grid using gridspec")
    parser.add_argument("--ncols", type=int, default=None,
                        help="Number of columns for subplots grid using "
                             "gridspec")
    
    # Plot Aesthetics for Waveforms
    parser.add_argument("-c", "--colors", nargs="+", type=str, default="k",
                        help="how to select color for waveforms:\n" 
                             "'k', 'red', '#d62728': one color for all "
                             "waveforms\n" 
                             "['k', 'r']: make a list of colors for each trace "
                             "must match the length of traces\n" 
                             "'C?': use the C0, C1, ... color cycle of MPL\n" 
                             "'viridis': use a colormap and plot discrete " 
                             "colors from the map based on length of stream"
                             )
    parser.add_argument("--alphas", nargs="+", type=float, default=None,
                        help="alpha of the time series line, number of inputs "
                             "must match the length of `fid`")
    parser.add_argument("-l", "--labels", nargs="+", type=str, default=None,
                        help="optional labels legend, must match len of `fid`")
    parser.add_argument("-lw", "--linewidth", nargs="?", type=float, 
                        default=0.5, help="linewidth of the time series line")
    parser.add_argument("--ylabel", nargs="?", type=str, default="Amplitude",
                        help="label for units, defaults to displacement")
    parser.add_argument("-y", "--ylim", nargs="+", type=float, default=None,
                        help="amplitude axis limits in s")
    parser.add_argument("--ylim_grid", nargs="+", type=float, default=None,
                        help="Only used if `wf_type`=='grid'. Set a different "
                             "y-axis limit for each subplot, given as a flat "
                             "list of ymin/ymax pairs, one pair per subplot, "
                             "e.g., '--ylim_grid -1 1 -2 2' sets the first "
                             "subplot to [-1, 1] and the second to [-2, 2]. "
                             "Must provide a pair for each subplot in the "
                             "grid. Overrides `--ylim` for grid plots")

    # Time axis (X-axis)
    parser.add_argument("-t", "--time", nargs="?", type=str, default="s",
                    help="units for x-axis/time axis. choice: 's'econds "
                            "(default), 'm'inutes, 'h'ours, 'a'bsolute (wip)."
                            "If using 'a' you may add '+i' or '-i' to "
                            "time shift the array, e.g., to go from UTC to "
                            "local time. E.g., 'a-7' will subtract 7 hours.")
    parser.add_argument("--minticks", type=int, default=3, 
                        help="min ticks if --time='a'")
    parser.add_argument("--maxticks", type=int, default=6, 
                        help="max ticks if --time='a'")
    parser.add_argument("-x", "--xlim", nargs="+", default=None, type=str,
                        help="time axis limits in s or if `time`=='a' then "
                             "values should be in datetime, see tmarks")
    parser.add_argument("-tm", "--tmarks", nargs="+", 
                        help="plot vertical lines at given relative times, "
                             "should match the units of `time`. If `time`=='a' "
                             "then each tmark should be a datetime "
                             "YYYY-MM-DDTHH:MM:SS")
    parser.add_argument("--tmarks_c", nargs="+", default="k",
                        help="colors for each of the time marks, should "
                             "either be single letter for all marks or match " 
                             "length of tmarks for individual colors")
    parser.add_argument("--group_vels", nargs="+", default=None, type=float,
                        help="plot time markers for given group velocities in "
                             "units km/s. Arrival times are calculated using "
                             "`tp_dist_deg` or `tp_dist_km` as well as "
                             "`tp_start`")

    # TauP Phase Arrivals
    parser.add_argument("--tp_phases", nargs="+", type=str, default=None,
                        help="taup phase arrivals, requires all 'tp_*' params")
    parser.add_argument("--tp_model", nargs="?", type=str, default="iasp91",
                        help="taup model, defaults to 'iasp91'")
    parser.add_argument("--tp_dist_km", nargs="?", type=float, default=None,
                        help="TauP source receiver distance in km")
    parser.add_argument("--tp_dist_deg", nargs="?", type=float, default=None,
                        help="TauP source receiver distance in degrees")
    parser.add_argument("--tp_depth", nargs="?", type=float, default=None,
                        help="TauP source depth km")
    parser.add_argument("--tp_start", nargs="?", type=str, default=None,
                        help="origintime of the event, must match the format "
                             "of --time. If not given, assumes first sample "
                             "of --xlim defines origintime")
    parser.add_argument("--tp_cmap", nargs="?", type=str, default="rainbow_r",
                        help="Matplotlib colormap that defines how to color "
                             "each of the TauP phases")
    parser.add_argument("--tp_plot_rays", action="store_true", default=False,
                        help="Create a separate plot of raypaths for the given "
                             "source parameters and phase list")
    
    # Spectrogram (all parameters end with _s)
    parser.add_argument("--spectrogram", action="store_true", default=False,
                        help="plot spectrogram of the raw trace. See all " 
                             "'sp_*' parameters to control the look")
    parser.add_argument("--sp_cmap", nargs="?", type=str, 
                        default="viridis",
                        help="colormap of the spectrogram")
    parser.add_argument("--sp_numcol", nargs="?", type=int, default=256,
                        help="number of colors in colormap of the spectrogram")
    parser.add_argument("--sp_dbscale_off", action="store_true",
                        help="Turn off dB scaling for colors/amplitudes in " 
                             "spectrogram")    
    parser.add_argument("--sp_logscale", action="store_true",
                        help="turn on log scale for spectrogram y-axis")
    parser.add_argument("--sp_clip", default=[0., 1.], type=float, nargs="+",
                        help=" adjust cmap to clip at lower and/or upper end")
    parser.add_argument("--sp_vmax", default=None, type=float,
                        help="maximum value threshold for spectrogram")
    parser.add_argument("--sp_interp", default="antialiased", type=str,
                        help="image interpolation of the spectrogram, "
                             "'antialiased' (default) or e.g. 'nearest'. "
                             "Not used with `sp_logscale`")
    parser.add_argument("--sp_vmin", default=None, type=float,
                        help="minimum value threshold for spectrogram, use "
                             "with `sp_vmax` to fix the color scale across "
                             "figures")
    parser.add_argument("--sp_wlen", type=float, default=None,
                        help="spectrogram window length in seconds. By "
                             "default the shortest window (>=128 samples) "
                             "giving at most ~2000 time columns is used, so "
                             "longer time ranges get longer windows")
    parser.add_argument("--sp_overlap", type=float, default=0.9,
                        help="fractional overlap of spectrogram windows, "
                             "0 to 1. Defaults to 0.9")
    parser.add_argument("--sp_idx", type=int, default=0,
                        help="Iff multiple waveforms plotted, choose which of "
                             "them to use for the spectrogram. Defaults to 0. " 
                             "order is alphabetical")
    
    # Plot additional time series
    parser.add_argument("--add_trace", nargs="?", type=str, default=None,
                        help="Path to plot an additional time series that does "
                             "not have the same sampling rate . "
                             "All --tr_* pars are associated. Data should be " 
                             "in 2-column ASCII or an ObsPy readable format")
    parser.add_argument("--tr_time", nargs="?", type=str, default="a",
                        help="Time of the time series [a]bsolute (UTC) or " \
                             "[r]elative to xlim[0] or start of main trace",
                        choices=["a", "r"])
    parser.add_argument("--tr_ylabel", type=str, default="",
                        help="Y-label text")
    parser.add_argument("--tr_label", type=str, default="",
                        help="Label text for shared legend")
    parser.add_argument("--tr_color", type=str, default="C0",
                       help="color of the time series")
    parser.add_argument("--tr_ylim", type=float, nargs=2, default=None,
                        help="y-limits (min max) of the additional time "
                             "series axis")
    parser.add_argument("--tr_front", action="store_true", default=False,
                        help="draw the additional time series in front of "
                             "the waveform instead of behind it")
    
    # Misc plotting options
    parser.add_argument("--fig_size", type=float, nargs="+", default=None,
                        help="Figure size")
    parser.add_argument("--dpi", type=float, nargs="?", default=None,
                        help="Dots per inch, defaults to Matplotlib's "
                             "default (100)")
    parser.add_argument("--fig_len", type=float, nargs="?", default=8,
                        help="Figure length, use with 'fig_asp', overrides "
                             "'--fig_size' if both given")
    parser.add_argument("--fig_asp", type=float, nargs="?", default=None,
                        help="Figure aspect ratio, use with 'fig_len', "
                             "overrides '--fig_size' if given")
    parser.add_argument("--no_legend", action="store_true",
                        help="Turn off waveform legend")
    parser.add_argument("--ncol_legend", type=int, default=1,
                        help="Number of columns in legend")
    parser.add_argument("--title", nargs="?", type=str, default=None,
                        help="title of the figure, defaults to ID and fmin/max")
    parser.add_argument("--title_append", nargs="?", type=str, 
                        default="", help="append to default title")
    parser.add_argument("--title_prepend", nargs="?", type=str, 
                        default="", help="prepend to default title")
    parser.add_argument("--subplot_label", nargs="?", type=str, default=None,
                        help="For publications, annotate a subplot label to "
                             "top right corner of the figure")
    parser.add_argument("--frameoff", action="store_true",
                        help="Turn everything off except the waveform")
    parser.add_argument("-s", "--save", type=str, default=None,
                        help="filename to save figure")
    parser.add_argument("-o", "--output", type=str, default=None,
                        help="optional path to output processed seismograms")
    parser.add_argument("--transparent", action="store_true", default=False,
                        help="Saved figure transparent/alpha background. If " 
                             "not set then background will be white")
    parser.add_argument("--noshow", action="store_true", default=False,
                        help="dont show the figure, default behavior will show")

    return parser.parse_args()


def _float_or_none(value):
    """
    Argparse type that accepts a float or the string 'none'

    :type value: str
    :param value: command line value
    :rtype: float or None
    :return: the value as a float, or None if 'none' was given
    """
    if str(value).lower() == "none":
        return None
    return float(value)


def find_nearest(array, value):
    """
    Find the nearest index in a consecutive (time) array given a chosen value
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def _nearest_pow_2(x):
    """
    Find power of two nearest to x

    >>> _nearest_pow_2(3)
    2.0
    >>> _nearest_pow_2(15)
    16.0

    :type x: float
    :param x: Number
    :rtype: int
    :return: Nearest power of 2 to x
    """
    a = math.pow(2, math.ceil(np.log2(x)))
    b = math.pow(2, math.floor(np.log2(x)))
    if abs(a - x) < abs(b - x):
        return a
    else:
        return b
    

def envelope(data, samp_rate, smooth=0.):
    """
    Envelope of a signal: the magnitude of its analytic signal from the
    Hilbert transform, optionally smoothed with a centred running mean

    :type data: np.array
    :param data: signal
    :type samp_rate: float
    :param samp_rate: sampling rate in Hz
    :type smooth: float
    :param smooth: running mean length in seconds, 0 for none
    :rtype: np.array
    :return: envelope, same length as `data`
    """
    data = np.asarray(data, dtype=float)
    # Pad the FFT to a fast length; the analytic signal is cut back after
    env = np.abs(hilbert(data, N=next_fast_len(len(data))))[:len(data)]
    n = int(round(smooth * samp_rate))
    if n > 1:
        env = uniform_filter1d(env, size=n, mode="nearest")
    return env


def spectrogram(ax, xvals, data, samp_rate, xrange=None, per_lap=0.9,
                wlen=None, log=False, dbscale=False, mult=8., cmap=None,
                ncolors=256, zorder=None, clip=[0.0, 1.0], vmax=None,
                vmin=None, max_cols=2000, max_rows=1024,
                interpolation="antialiased"):
    """
    Computes and plots the spectrogram of a time series on a shared time axis.

    .. note::
        Modified from ObsPy's spectrogram function
        https://docs.obspy.org/_modules/obspy/imaging/spectrogram.html\
            #spectrogram

        Only the displayed time range is computed, and the number of
        columns and rows is capped so that hours of data stay fast to
        compute and draw. A linear frequency axis is drawn with imshow;
        a log frequency axis uses pcolormesh on log-spaced rows.

    :type ax: matplotlib.axes.Axes
    :param ax: axis to plot into, sharing its x-axis with the waveform
    :type xvals: np.ndarray
    :param xvals: x-axis value of every sample in `data` (seconds, minutes,
        hours or Matplotlib date numbers)
    :type data: np.ndarray
    :param data: time series
    :type samp_rate: float
    :param samp_rate: sampling rate of `data` in Hz
    :type xrange: tuple of float or None
    :param xrange: (xmin, xmax) in `xvals` units to compute the spectrogram
        over, e.g. the plot's x-limits. None uses all of `data`
    :type per_lap: float
    :param per_lap: fraction of overlap of the sliding window, 0 to 1
    :type wlen: float or None
    :param wlen: window length in seconds. If None, the shortest power-of-two
        window of at least 128 samples is used that gives no more than
        `max_cols` columns, so long time ranges get longer windows (and
        finer frequency resolution) instead of millions of columns
    :type log: bool
    :param log: logarithmic frequency axis if True, linear otherwise
    :type dbscale: bool
    :param dbscale: if True 10 * log10 of the power is plotted, if False the
        square root
    :type mult: float
    :param mult: zero-pad each window to `mult` times its length for a
        smoother look, limited so there are at most `max_rows` frequencies
    :type cmap: str or None
    :param cmap: colormap name, defaults to ObsPy's sequential colormap
    :type ncolors: int
    :param ncolors: number of discrete colors in the colormap
    :type zorder: float
    :param zorder: zorder of the image
    :type clip: [float, float]
    :param clip: clip the color range at these fractions of the amplitude
        range (linear or dB depending on `dbscale`)
    :type vmax: float or None
    :param vmax: absolute upper color limit, overrides the upper `clip`
    :type vmin: float or None
    :param vmin: absolute lower color limit, overrides the lower `clip`
    :type max_cols: int
    :param max_cols: maximum number of time columns
    :type max_rows: int
    :param max_rows: maximum number of frequency rows
    :type interpolation: str
    :param interpolation: Matplotlib imshow interpolation of the linear
        frequency axis. 'antialiased' (default) averages rows and columns
        when there are more than screen pixels; 'nearest' draws only some
        of them, which looks speckled
    :rtype: matplotlib.cm.ScalarMappable
    :return: the plotted image, for a colorbar
    """
    if not cmap:
        cmap = obspy_sequential
    else:
        cmap = plt.get_cmap(cmap, ncolors)
    samp_rate = float(samp_rate)
    xvals = np.asarray(xvals)
    data = np.asarray(data, dtype=float)
    # x units per second, from the full span to avoid date round-off
    units_per_sec = (xvals[-1] - xvals[0]) / ((len(xvals) - 1) / samp_rate)

    # Samples inside the displayed range
    if xrange is not None:
        i0 = max(0, int(np.searchsorted(xvals, xrange[0])))
        i1 = min(len(data), int(np.searchsorted(xvals, xrange[1])) + 1)
    else:
        i0, i1 = 0, len(data)
    npts_show = i1 - i0

    # Window length: user choice, or grow it until the column count fits
    if wlen:
        nfft = int(_nearest_pow_2(wlen * samp_rate))
    else:
        nfft = 128
        while npts_show / (nfft * (1 - per_lap)) > max_cols and \
                nfft < npts_show // 4:
            nfft *= 2
    nlap = int(nfft * float(per_lap))
    # A user-set window may still give too many columns; skip ahead instead
    step = nfft - nlap
    if (npts_show - nfft) / step > max_cols:
        nlap = max(0, nfft - int(np.ceil((npts_show - nfft) / max_cols)))
        step = nfft - nlap

    # Pad half a window on each side (where data exists) so the first and
    # last columns are centred at the edges of the displayed range
    i0 = max(0, i0 - nfft // 2)
    i1 = min(len(data), i1 + nfft // 2)
    data = data[i0:i1]
    npts = len(data)
    if npts < nfft:
        msg = (f"Input signal too short ({npts} samples, window length "
               f"{nfft / samp_rate} seconds, nfft {nfft} samples, sampling "
               f"rate {samp_rate} Hz)")
        raise ValueError(msg)

    pad_to = nfft
    if mult is not None:
        pad_to = int(_nearest_pow_2(mult)) * nfft
        while pad_to // 2 > max_rows and pad_to > nfft:
            pad_to //= 2

    specgram, freq, time = mlab.specgram(data - data.mean(), Fs=samp_rate,
                                         NFFT=nfft, pad_to=pad_to,
                                         noverlap=nlap)
    if len(time) < 2:
        msg = (f"Input signal too short ({npts} samples, window length "
               f"{nfft / samp_rate} seconds, nfft {nfft} samples, {nlap} "
               f"samples window overlap, sampling rate {samp_rate} Hz)")
        raise ValueError(msg)
    logger.info(f"spectrogram: {nfft / samp_rate:g} s windows, "
                f"{len(time)} columns x {len(freq) - 1} rows")

    # db scale and remove zero/offset for amplitude
    if dbscale:
        with np.errstate(divide="ignore"):
            specgram = 10 * np.log10(specgram[1:, :])
        finite = np.isfinite(specgram)
        specgram[~finite] = specgram[finite].min()
    else:
        specgram = np.sqrt(specgram[1:, :])
    freq = freq[1:]

    vmin_frac, vmax_frac = clip
    if vmin_frac < 0 or vmax_frac > 1 or vmin_frac >= vmax_frac:
        msg = "Invalid parameters for clip option."
        raise ValueError(msg)
    _range = float(specgram.max() - specgram.min())
    if vmin is None:
        vmin = specgram.min() + clip[0] * _range
    if vmax is None:
        vmax = specgram.min() + vmax_frac * _range
    norm = Normalize(vmin, vmax, clip=True)

    # Bin edges: mlab gives window centres (seconds after the first sample
    # of `data`) and bin-centre frequencies
    halfbin_time = (time[1] - time[0]) / 2.0
    halfbin_freq = (freq[1] - freq[0]) / 2.0
    t_edges = np.concatenate((time - halfbin_time,
                              [time[-1] + halfbin_time]))
    x_edges = xvals[i0] + t_edges * units_per_sec
    f_edges = np.concatenate((freq - halfbin_freq,
                              [freq[-1] + halfbin_freq]))

    kwargs = {"cmap": cmap, "zorder": zorder, "norm": norm}
    if log:
        # Sample the rows on a log-spaced grid so pcolormesh stays small
        nrow = min(len(freq), 400)
        # Frequencies below 1 / window length are not resolved (only zero
        # padding), so the log axis starts there
        f_low = max(f_edges[0], samp_rate / nfft)
        f_log = np.geomspace(f_low, f_edges[-1], nrow + 1)
        centres = np.sqrt(f_log[:-1] * f_log[1:])
        rows = np.clip(np.searchsorted(f_edges, centres) - 1, 0,
                       len(freq) - 1)
        im = ax.pcolormesh(x_edges, f_log, specgram[rows], rasterized=True,
                           **kwargs)
        ax.set_yscale("log")
    else:
        im = ax.imshow(specgram, origin="lower", aspect="auto",
                       interpolation=interpolation,
                       extent=(x_edges[0], x_edges[-1], f_edges[0],
                               f_edges[-1]), **kwargs)
    return im


def fast_remove_response(tr, inventory, output="VEL", water_level=60.,
                         pre_filt=None, zero_mean=True, taper=True,
                         taper_fraction=0.05, n_eval=2 ** 15):
    """
    Remove the instrument response of a Trace in place, like ObsPy's
    Trace.remove_response(), but much faster for long traces.

    ObsPy evaluates the full response with evalresp at every FFT frequency
    (about one million points for an hour of 100 Hz data), which dominates
    the run time. The response is smooth, so here it is evaluated at
    about `n_eval` frequencies (evenly spaced, plus log-spaced points at low
    frequency) and interpolated (real and imaginary parts) onto the FFT
    grid. Every other step (demean, taper, pre-filter,
    water level, inverse FFT) follows ObsPy 1.5 exactly. Traces short
    enough to need fewer than `n_eval` frequencies, and polynomial
    responses, are passed straight to ObsPy.

    :type tr: obspy.core.trace.Trace
    :param tr: trace to correct, modified in place
    :type inventory: obspy.core.inventory.Inventory
    :param inventory: inventory containing the trace's response
    :type output: str
    :param output: 'DISP', 'VEL', 'ACC' or 'DEF'
    :type water_level: float or None
    :param water_level: water level in dB, None to turn it off
    :type pre_filt: list of float or None
    :param pre_filt: four corner frequencies of the frequency-domain taper
    :type zero_mean: bool
    :param zero_mean: demean the data before deconvolution
    :type taper: bool
    :param taper: apply a cosine taper before deconvolution
    :type taper_fraction: float
    :param taper_fraction: fraction of the trace tapered at each end
    :type n_eval: int
    :param n_eval: number of frequencies at which the response is evaluated
    :rtype: obspy.core.trace.Trace
    :return: the corrected trace
    """
    from obspy.core.inventory import PolynomialResponseStage
    from obspy.signal.invsim import (cosine_taper, cosine_sac_taper,
                                     invert_spectrum)
    from obspy.signal.util import _npts2nfft

    kwargs = dict(output=output, water_level=water_level, pre_filt=pre_filt,
                  zero_mean=zero_mean, taper=taper,
                  taper_fraction=taper_fraction)
    response = tr._get_response(inventory)
    npts = tr.stats.npts
    nfft = _npts2nfft(npts)
    if (not response.response_stages
            or isinstance(response.response_stages[0],
                          PolynomialResponseStage)
            or nfft // 2 + 1 <= n_eval):
        return tr.remove_response(inventory=inventory, **kwargs)

    data = tr.data.astype(np.float64)
    if zero_mean:
        data -= data.mean()
    if taper:
        data *= cosine_taper(npts, taper_fraction, sactaper=True,
                             halfcosine=False)
    data = np.fft.rfft(data, n=nfft)
    freqs = np.fft.rfftfreq(nfft, d=tr.stats.delta)

    # Coarse response, interpolated onto the FFT frequencies. Linear
    # spacing covers the passband; extra log-spaced points resolve the steep
    # response near 0 Hz, which matters for DISP without `pre_filt`
    f_eval = np.union1d(np.linspace(0, freqs[-1], n_eval),
                        np.geomspace(freqs[1], freqs[-1], n_eval // 4))
    resp_eval = response.get_evalresp_response_for_frequencies(
        f_eval, output=output)
    freq_response = (np.interp(freqs, f_eval, resp_eval.real)
                     + 1j * np.interp(freqs, f_eval, resp_eval.imag))

    if pre_filt:
        data *= cosine_sac_taper(freqs, flimit=pre_filt)
    if water_level is None:
        freq_response[0] = 0.0
        freq_response[1:] = 1.0 / freq_response[1:]
    else:
        invert_spectrum(freq_response, water_level)
    data *= freq_response
    data[-1] = abs(data[-1]) + 0.0j
    tr.data = np.fft.irfft(data)[0:npts]
    return tr


def merge_stream(st, fill="interpolate"):
    """
    Merge traces with the same ID into one continuous trace each, e.g., the
    segments of a gappy file, or consecutive day files of one channel.
    Traces keep the order in which their IDs first appear in `st`.

    :type st: obspy.core.stream.Stream
    :param st: stream to merge
    :type fill: str or float
    :param fill: how to fill gaps: 'interpolate' (linear between the samples
        on either side), 'latest' (repeat the last sample before the gap),
        or a number to fill with
    :rtype: obspy.core.stream.Stream
    :return: merged stream
    """
    if fill not in ("interpolate", "latest"):
        try:
            fill = float(fill)
        except (TypeError, ValueError):
            raise ValueError(f"`merge_fill` must be 'interpolate', 'latest' "
                             f"or a number, not {fill!r}") from None

    order = []
    for tr in st:
        if tr.id not in order:
            order.append(tr.id)

    gaps = st.get_gaps()
    st = st.copy()
    # A numeric fill is a float, which ObsPy cannot write into integer data
    # (e.g., raw counts), so make the data float first
    if not isinstance(fill, str):
        for tr in st:
            tr.data = tr.data.astype(np.float64)
    st.merge(method=1, fill_value=fill)
    if gaps:
        total = sum(max(gap[6], 0) for gap in gaps)
        logger.info(f"merged {len(gaps)} gap(s)/overlap(s), {total:.2f} s "
                    f"of gaps filled with {fill!r}")

    return Stream([tr for id_ in order for tr in st.select(id=id_)])


def convert_timezone(code, st):
    """
    When plotting in absolute time (args.time == "a*"), we allow time 
    shifting  by time zone to get to the correct time local time. Returns 
    streams with converted time

    .. note::

        Confusingly, trying to input time zones in UTCDateTime like
        UTCDateTime("YYYY-MM-DDTHH:MM:SS-09) 
        assumes that the given date is in local time, and that the modifier -09
        is the conversion to UTC. But since seismic data is IN UTC, we cannot
        use this because it shifts us in the wrong way.

    :type code: str
    :param code: e.g., +08 to shift forward by 8 hours
    """
    assert(len(code) == 3), f"must be +?? or -??, not {code}"
    assert(code[0] in ["+", "-"])

    sign = code[0]
    shift = int(code[1:]) * 60 * 60  # hours -> seconds

    for tr in st:
        if sign == "+":
            tr.stats.starttime += shift
        elif sign == "-":
            tr.stats.starttime -= shift

    return st

    
def set_plot_aesthetic(
        ax, xtick_labels=True, ytick_labels=True,
        ytick_fontsize=FONTSIZES["ytick_label"],
        xtick_fontsize=FONTSIZES["xtick_label"], tick_linewidth=1.5,
        tick_length=5., tick_direction="in", ytick_format="sci",
        xlabel_fontsize=FONTSIZES["axis_label"],
        ylabel_fontsize=FONTSIZES["axis_label"], axis_linewidth=1.5,
        spine_zorder=500, title_fontsize=FONTSIZES["title"],
        spine_top=True, spine_bot=True, spine_left=True, spine_right=True, 
        xtick_minor=None, xtick_major=None, ytick_minor=None, ytick_major=None,
        xgrid_major=True, xgrid_minor=True, ygrid_major=True, ygrid_minor=True,
        **kwargs):
    """
    Set a uniform look for figures. Adapted from
    :func:`pysep.utils.plot.set_plot_aesthetic` but with PrettyPlot's own
    defaults (see `FONTSIZES`), so the two are kept separate.
    """
    ax.title.set_fontsize(title_fontsize)
    ax.tick_params(axis="both", which="both", width=tick_linewidth,
                   direction=tick_direction, length=tick_length)
    ax.tick_params(axis="x", labelsize=xtick_fontsize)
    ax.tick_params(axis="y", labelsize=ytick_fontsize)
    ax.xaxis.label.set_size(xlabel_fontsize)
    ax.yaxis.label.set_size(ylabel_fontsize)

    # Thicken up the bounding axis lines
    for axis, flag in zip(["top", "bottom", "left", "right"],
                          [spine_top, spine_bot, spine_left, spine_right]):
        # Deal with the case where command line users are inputting strings
        if isinstance(flag, str):
            flag = bool(flag.capitalize() == "True")
        ax.spines[axis].set_visible(flag)
        ax.spines[axis].set_linewidth(axis_linewidth)

    # Set spines above azimuth bins
    for spine in ax.spines.values():
        spine.set_zorder(spine_zorder)

    # Scientific format for Y-axis
    if ytick_format == "sci":
        try:
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
            ax.yaxis.get_offset_text().set_fontsize(FONTSIZES["yaxis_offset"])
        except AttributeError:
            # If we are in log format axis this will not work
            pass
    else:
        try:
            ax.ticklabel_format(axis="y", style=ytick_format)
        except AttributeError:
            pass

    # Set xtick label major and minor which is assumed to be a time series
    if xtick_major:
        ax.xaxis.set_major_locator(MultipleLocator(float(xtick_major)))
    if xtick_minor:
        ax.xaxis.set_minor_locator(MultipleLocator(float(xtick_minor)))
    if ytick_major:
        ax.yaxis.set_major_locator(MultipleLocator(float(ytick_major)))
    if ytick_minor:
        ax.yaxis.set_minor_locator(MultipleLocator(float(ytick_minor)))

    plt.sca(ax)

    if xgrid_major:
        plt.grid(visible=True, which="major", axis="x", alpha=0.2, linewidth=1)
    if xgrid_minor:
        plt.grid(visible=True, which="minor", axis="x", alpha=0.2, linewidth=.5)
    if ygrid_major:
        plt.grid(visible=True, which="major", axis="y", alpha=0.2, linewidth=1)
    if ygrid_minor:
        plt.grid(visible=True, which="minor", axis="y", alpha=0.2, linewidth=.5)


def _set_xaxis_obspy_dates(ax, ticklabels_small=True, minticks=3, maxticks=6):
    """
    Set Formatter/Locator of x-Axis to use ObsPyAutoDateFormatter and do some
    other tweaking.

    In contrast to normal matplotlib ``AutoDateFormatter`` e.g. shows full
    timestamp on first tick when zoomed in so far that matplotlib would only
    show hours or minutes on all ticks (making it impossible to tell the date
    from the axis labels) and also shows full timestamp in matplotlib figures
    info line (mouse-over info of current cursor position).

    :type ax: :class:`matplotlib.axes.Axes`
    :rtype: None
    """
    ax.xaxis_date()
    locator = AutoDateLocator(minticks=minticks, maxticks=maxticks)
    locator.intervald[MINUTELY] = [1, 2, 5, 10, 15, 30]
    locator.intervald[SECONDLY] = [1, 2, 5, 10, 15, 30]
    ax.xaxis.set_major_formatter(ObsPyAutoDateFormatter(locator))
    ax.xaxis.set_major_locator(locator)
    if ticklabels_small:
        plt.setp(ax.get_xticklabels(), fontsize=FONTSIZES["xtick_date_label"])

def cmaphex(nvals, cmap="seismic"):
    """
    Return a list of hex codes for `nvals` of a `cmap` which can be used to
    provide a sequential colormap made up of discrete color values
    """
    cmap = plt.get_cmap(cmap, nvals)
    hex_out = []
    for i in range(cmap.N):
        rgba = cmap(i)
        # rgb2hex accepts rgb or rgba
        hex_out.append(mpl.colors.rgb2hex(rgba))
    return hex_out


class PrettyPlot():
    """
    Command line and scriptable waveform and spectrogram plotter. Every
    command line option of ``pplot`` is a keyword argument; call
    :meth:`PrettyPlot.run` to make the figure.
    """
    def __init__(self, fids, idx=0,
                 # Processing
                 taper=0, fmin=None, fmax=None, zerophase=False, corners=4,
                 resample=False, t0=0, tstart=0, detrend=False, integrate=0,
                 merge=False, merge_fill="interpolate",
                 differentiate=0, trim_pct=1.,
                 # Response removal
                 inv=None, rr_output="VEL", rr_water_level=60.,
                 rr_pre_filt=None, rr_zero_mean=True, rr_taper=True,
                 rr_taper_fraction=0.05, rr_exact=False,
                 # Waveform plotting
                 wf_abs=False, wf_type="default", wf_norm=False, wf_scale=1,
                 wf_env=False, wf_env_smooth=0., wf_env_color=None,
                 wf_recsec_spacing=1, wf_spacing_exact=False, wf_order=None,
                 wf_grid=None, nrows=None, ncols=None,
                 # Plotting Aesthetics
                 colors="k", alphas=None, labels=None, linewidth=0.5, 
                 ylabel="Amplitude", ylim=None, ylim_grid=None, xlim=None,
                 # Time Axis
                 time="s", minticks=3, maxticks=6, tmarks=None, 
                 tmarks_c="k", group_vels=None,
                 # TauP 
                 tp_phases=None, tp_model="iasp91", tp_dist_km=None,
                 tp_dist_deg=None, tp_depth=None, tp_start=None, 
                 tp_cmap="rainbow_r", tp_plot_rays=False,
                 # Spectrogram
                 spectrogram=False, sp_cmap="viridis", sp_numcol=256,
                 sp_dbscale=True, sp_logscale=False, sp_clip=None,
                 sp_vmax=None, sp_vmin=None, sp_interp="antialiased",
                 sp_idx=0, sp_wlen=None,
                 sp_overlap=0.9,
                 # Additional Time Series
                 add_trace=None, tr_time="a", tr_label="", tr_ylabel="", 
                 tr_color="C0", tr_ylim=None, tr_front=False,
                 # Misc.
                 fig_size=None, dpi=100, fig_len=8, fig_asp=None,
                 legend=True, ncol_legend=1, title=None, title_append="",
                 title_prepend="", subplot_label=None,
                 save=None, output=None, show=True, frameoff=False, 
                 transparent=False,
                 **kwargs
                 ):
        """Input parameters, see argparser for descriptions"""
        # Waveforms
        self.fids = fids
        self.idx = idx
        self.merge = merge
        self.merge_fill = merge_fill

        # Preprocessing
        self.taper = taper
        self.fmin = fmin
        self.fmax = fmax
        self.zerophase = zerophase
        self.corners = corners
        self.resample = resample
        self.t0 = t0
        self.tstart = tstart
        self.detrend = detrend
        self.integrate = integrate
        self.differentiate = differentiate
        self.trim_pct = trim_pct

        # Instrument response removal, `inv` may be StationXML path(s) or an
        # ObsPy Inventory
        self.inv = inv
        self.rr_output = rr_output
        self.rr_water_level = rr_water_level
        self.rr_pre_filt = rr_pre_filt
        self.rr_zero_mean = rr_zero_mean
        self.rr_taper = rr_taper
        self.rr_taper_fraction = rr_taper_fraction
        self.rr_exact = rr_exact

        # Plot structure control
        self.wf_abs = wf_abs
        self.wf_env = wf_env
        self.wf_env_smooth = wf_env_smooth
        self.wf_env_color = wf_env_color
        self.wf_type = wf_type
        self.wf_norm = wf_norm
        self.wf_scale = wf_scale
        self.wf_recsec_spacing = wf_recsec_spacing
        self.wf_spacing_exact = wf_spacing_exact
        self.wf_order = wf_order
        self.nrows = nrows
        self.ncols = ncols
        # `wf_grid` is shorthand for (nrows, ncols)
        if wf_grid:
            assert len(wf_grid) == 2, \
                f"`wf_grid` must be two integers (nrows ncols), not {wf_grid}"
            self.nrows, self.ncols = [int(_) for _ in wf_grid]

        # Waveform style parameters
        self.colors = colors
        self.alphas = alphas
        self.labels = labels
        self.linewidth = linewidth
        self.ylim = ylim
        self.ylim_grid = ylim_grid
        self.ylabel = ylabel

        # Time axis
        self.time = time
        self.maxticks = maxticks
        self.minticks = minticks
        self.xlim = xlim
        self.tmarks = tmarks
        self.tmarks_c = tmarks_c
        self.group_vels = group_vels

        # (optional) TauP arrival times
        self.tp_phases = tp_phases
        self.tp_model = tp_model
        self.tp_dist_km = tp_dist_km
        self.tp_dist_deg = tp_dist_deg
        self.tp_depth = tp_depth
        self.tp_start = tp_start
        self.tp_cmap = tp_cmap
        self.tp_plot_rays = tp_plot_rays

        # (optional) Spectrograms
        self.spectrogram = spectrogram
        self.sp_cmap = sp_cmap
        self.sp_numcol = sp_numcol
        self.sp_dbscale = sp_dbscale
        self.sp_logscale = sp_logscale
        self.sp_clip = sp_clip or [0., 1.]
        self.sp_vmax = sp_vmax
        self.sp_vmin = sp_vmin
        self.sp_interp = sp_interp
        self.sp_idx = sp_idx
        self.sp_wlen = sp_wlen
        self.sp_overlap = sp_overlap

        # Additional time series on shared-X axis
        self.add_trace = add_trace
        self.tr_time = tr_time
        self.tr_label = tr_label
        self.tr_ylabel = tr_ylabel
        self.tr_color = tr_color
        self.tr_ylim = tr_ylim
        self.tr_front = tr_front

        # Figure output control
        self.fig_size = fig_size
        self.dpi = dpi
        self.fig_len = fig_len
        self.fig_asp = fig_asp
        self.legend = legend
        self.ncol_legend = ncol_legend
        self.title = title
        self.title_append = title_append
        self.title_prepend = title_prepend
        self.subplot_label = subplot_label
        self.save = save
        self.output = output
        self.show = show
        self.frameoff = frameoff
        self.transparent = transparent

        # Internal attributes for plotting functions
        self._xvals = None
        self._data = None
        self._labels = None 
        self._colors = None
        self._alphas = None
        self._ylims = None
        self._ref_time = None
        self._tz_shift = 0.

        # Parameter conflicts
        if self.wf_type == "grid":
            assert(not self.spectrogram), \
                f"spectrogram option cannot be used with grid style plots"

        # Populate Stream object
        self.st = Stream()
        for fid in self.fids:
            try:
                st = read(fid)
            except TypeError:
                # Not ObsPy-readable, assume SPECFEM two-column ASCII
                st = read_sem(fid)
            logger.info(f"reading {fid}")
            if len(st) > 1 and not self.merge:
                logger.warning(f"{len(st)} traces found, taking index "
                               f"{self.idx}. Use `merge` to join gappy data")
                st = st[self.idx]
            self.st += st
        if self.merge:
            self.st = merge_stream(self.st, fill=self.merge_fill)
        logger.info(f"read {self.st.__str__(extended=True)}")


    def setup_plot(self):
        """
        Set up the plot based on input parameters, choose how the figure and
        axes visually appears
        """
        self.f = None
        self.ax = None
        self.axs = []
        self.ax_spectra = None
        self.twax = None

        # Spectrogram + waveform plot
        if self.spectrogram:
            if self.fig_len and self.fig_asp:
                figsize = (self.fig_len, self.fig_len * self.fig_asp)
            else:
                figsize = self.fig_size or (8, 6)  # default

            logger.info("setting up waveform and spectrogram plot")
            self.f, axs = plt.subplots(2, dpi=self.dpi, figsize=figsize, 
                                       sharex=True)
            self.f.subplots_adjust(hspace=0)
            self.ax_spectra, self.ax = axs  # waveform on the bottom
            self.axs = [self.ax] * len(self.st)  # for consistency with grid
        # Grid style with subplots, no `ax`
        elif self.wf_type == "grid":
            if self.wf_type == "grid":
                self.nrows = self.nrows or len(self.st)
                self.ncols = self.ncols or 1
            figsize = self.fig_size or (7 * self.ncols, 4 * self.nrows)
            self.f = plt.figure(figsize=figsize, dpi=self.dpi)
            gs = gridspec.GridSpec(self.nrows, self.ncols, figure=self.f, 
                                   hspace=0.2, wspace=0.2)
            num_subplots = min(self.nrows * self.ncols, len(self.st))
            # Set the first axis then share x axis with all others
            ax0 = self.f.add_subplot(gs[0])
            self.axs.append(ax0)
            for i in range(1, num_subplots):
                ax = self.f.add_subplot(gs[i // self.ncols, i % self.ncols],
                                        sharex=ax0)
                self.axs.append(ax)
            self.ax = self.axs[0]  # for consistency with other plots
        # Waveform only plot
        else:
            if self.fig_len and self.fig_asp:
                figsize = (self.fig_len, self.fig_len * self.fig_asp)
            else:
                figsize = self.fig_size  or (7, 4)  # default

            logger.info("setting up waveform plot")
            self.f, self.ax = plt.subplots(dpi=self.dpi, figsize=figsize)
            self.axs = [self.ax] * len(self.st)  # for consistency with grid


    def trim_waveform(self):
        """
        Trim the waveform data so we don't plot the whole thing. Has some 
        functionality to add buffer to the trimmed tail ends of the data for
        FFT style processing which requires stable ends.
        """         
        if self.resample:
            logger.info(f"resampling to {self.resample}")
            self.st.resample(sampling_rate=self.resample)
       
        # Shift to correct time axis
        if self.time.startswith("a") and len(self.time) > 1:
            _start = self.st[0].stats.starttime
            self.st = convert_timezone(code=self.time[1:], st=self.st)
            # Kept so response removal can match inventory epochs in UTC
            self._tz_shift = self.st[0].stats.starttime - _start

        # Relative time axes (and `xlim`) are measured from the first sample
        # of the first trace as read, so remember it before any trimming
        self._ref_time = self.st[0].stats.starttime

        # Trim data shorter so we don't process the entire waveform but add some
        # buffer so that preprocessing on the tails of the data doesn't show up.
        # the tail data will stay there but we will set the xlim of the plot to 
        # not show it
        if self.xlim:    
            if self.time.startswith("a"):
                start = UTCDateTime(self.xlim[0])
                end = UTCDateTime(self.xlim[1])
                buffer = (end - start) * self.trim_pct  # some % of the record
                self.st.trim(starttime=start - buffer, endtime=end + buffer)
            else:
                # `xlim` is in the units of `time` and includes the t0 and
                # tstart offsets, so undo those to get seconds after the
                # reference time
                x0, x1 = [self._rel_to_seconds(float(_)) for _ in self.xlim]
                buffer = (x1 - x0) * self.trim_pct
                self.st.trim(starttime=self._ref_time + (x0 - buffer),
                             endtime=self._ref_time + (x1 + buffer))
            logger.info(f"trimming with {buffer}s buffer on either end")

            if not self.st:
                logger.critical("Trimming removed all data, please check "
                                "`xlim` values")
                sys.exit()

    def _seconds_per_unit(self):
        """
        Return the number of seconds in one unit of a relative time axis

        :rtype: float
        :return: 1 for 's', 60 for 'm', 3600 for 'h'
        """
        return {"s": 1., "m": 60., "h": 3600.}.get(self.time, 1.)

    def _rel_to_seconds(self, x):
        """
        Convert a value on the relative time axis to seconds after the
        reference time (first sample of the first trace as read)

        :type x: float
        :param x: time in the units of `time`, including t0/tstart offsets
        :rtype: float
        :return: seconds after `_ref_time`
        """
        return (x + self.t0 - self.tstart) * self._seconds_per_unit()

    def _default_origin(self):
        """
        Origin time used for TauP arrivals and group velocities when
        `tp_start` is not given: the start of `xlim` if set, otherwise the
        first plotted sample. Matches the units of the time axis.

        :rtype: UTCDateTime or float
        :return: UTCDateTime if `time` is absolute, else relative time
        """
        if self.time.startswith("a"):
            if self.xlim:
                return UTCDateTime(self.xlim[0])
            return self.st[0].stats.starttime
        if self.xlim:
            return float(self.xlim[0])
        return float(self._xvals[0])

    def remove_response(self):
        """
        Remove the instrument response with ObsPy's Stream.remove_response()
        using the StationXML given by `inv`. Runs after trimming, so only the
        (buffered) plot window is deconvolved, and before all other
        processing. If no `ylabel` was set, it is updated to match the
        output units.
        """
        if isinstance(self.inv, Inventory):
            inv = self.inv
        else:
            fids = [self.inv] if isinstance(self.inv, str) else self.inv
            inv = Inventory()
            for fid in fids:
                inv += read_inventory(fid)

        logger.info(f"removing response, output={self.rr_output}, "
                    f"water_level={self.rr_water_level}, "
                    f"pre_filt={self.rr_pre_filt}")

        # Inventory epochs are in UTC, so undo any time zone shift while
        # looking up the response
        for tr in self.st:
            tr.stats.starttime -= self._tz_shift
        kwargs = dict(output=self.rr_output, water_level=self.rr_water_level,
                      pre_filt=self.rr_pre_filt, zero_mean=self.rr_zero_mean,
                      taper=self.rr_taper,
                      taper_fraction=self.rr_taper_fraction)
        try:
            for tr in self.st:
                if self.rr_exact:
                    tr.remove_response(inventory=inv, **kwargs)
                else:
                    fast_remove_response(tr, inv, **kwargs)
        finally:
            for tr in self.st:
                tr.stats.starttime += self._tz_shift

        # Label the amplitude axis with physical units unless the user
        # already chose a label
        units = {"DISP": "Displacement [m]", "VEL": "Velocity [m/s]",
                 "ACC": r"Acceleration [m/s$^2$]"}
        if self.rr_output in units and \
                str(self.ylabel).lower() == "amplitude":
            self.ylabel = units[self.rr_output]

    def process_waveforms(self):
        """
        Process the waveforms based on input parameters with options for 
        filtering, detrending, tapering etc.
        """
        # Downsample all waveforms to the lowest sampling rate
        lowest_samp_rate = min([tr.stats.sampling_rate for tr in self.st])
        for tr in self.st:
            if tr.stats.sampling_rate > lowest_samp_rate:
                logger.info(f"resampling {tr.get_id()} to {lowest_samp_rate}")
                tr.resample(lowest_samp_rate)

        # Trim start and endtimes
        latest_starttime = max([tr.stats.starttime for tr in self.st])
        earliest_endtime = min([tr.stats.endtime for tr in self.st])
        self.st.trim(latest_starttime, earliest_endtime)

        if self.detrend:
            logger.info(f"detrending with '{self.detrend}'")
            self.st.detrend(type=self.detrend)
        taper = self.taper
        if self.integrate or self.differentiate:
            if self.taper == 0:
                taper = 0.05
                logger.info(f"setting taper pct {taper}")
        if taper:
            self.st.taper(taper)
            logger.info(f"tapering trace {taper * 100}%")
        if self.integrate:
            for i in range(self.integrate):
                logger.info("integrating trace")
                self.st.integrate()
        if self.differentiate:
            for i in range(self.differentiate):
                logger.info("differentiating trace")
                self.st.differentiate()
        if self.fmin and self.fmax:
            logger.info(f"filtering between {self.fmin} and {self.fmax} Hz, "
                        f"zerophase={self.zerophase}, corners={self.corners}")
            self.st.filter(
                "bandpass", freqmin=self.fmin, freqmax=self.fmax,
                zerophase=self.zerophase, corners=self.corners
                )
        elif self.fmin:
            logger.info(f"highpass filtering at {self.fmin} Hz, "
                        f"zerophase={self.zerophase}, corners={self.corners}")
            self.st.filter(
                "highpass", freq=self.fmin, zerophase=self.zerophase,
                corners=self.corners
                )
        elif self.fmax:
            logger.info(f"lowpass filtering at {self.fmax} Hz, "
                        f"zerophase={self.zerophase}, corners={self.corners}")
            self.st.filter(
                "lowpass", freq=self.fmax, zerophase=self.zerophase,
                corners=self.corners
                )
        
    def get_plot_parameters(self):
        """
        Establish what will be plotted and the styles for each of them. The data
        itself will control things like stacking or record sections. The final
        internal attributes created are: xvals, data, labels, colors, alphas
        """
        # TIME AXIS
        # assumes that all traces have the same start and end times
        times = self.st[0].times()
        if self.time.startswith("a"):
            xvals = ((times / SECONDS_PER_DAY) + 
                     date2num(self.st[0].stats.starttime.datetime))
            for ax in self.axs:
                _set_xaxis_obspy_dates(ax, minticks=self.minticks, 
                                       maxticks=self.maxticks)
        else:
            # Measure from the reference time so trimming (by `xlim` or to
            # a common window) does not reset the axis to zero
            if self._ref_time is None:
                self._ref_time = self.st[0].stats.starttime
            offset = self.st[0].stats.starttime - self._ref_time
            xvals = (times + offset) / self._seconds_per_unit()
            xvals -= self.t0
            xvals += self.tstart

        # DATA ASSIMILATION
        # Determine the order in which to plot the waveforms, either by stream
        # order (alphabetical) or by user-defined through `wf_order`
        if not self.wf_order:
            data = np.array([tr.data for tr in self.st])
            labels = [tr.get_id() for tr in self.st]
            srs = [tr.stats.sampling_rate for tr in self.st]
        else:
            # Allow arbitrary order set by User
            data, labels, srs = [], [], []
            for i in self.wf_order:
                data.append(self.st[i].data)
                labels.append(self.st[i].get_id())
                srs.append(self.st[i].stats.sampling_rate)

        # DATA MANIPULATION
        # Normalize data by its own max so everything falls between [-1, 1]
        if self.wf_norm:
            for i, d in enumerate(data[:]):
                data[i] = d / np.abs(d).max()
        # Envelopes are taken before `wf_abs`, which would distort them
        self._env = None
        if self.wf_env:
            if "recsec" in self.wf_type:
                logger.warning("`wf_env` is drawn without record section "
                               "offsets")
            self._env = [envelope(d, sr, self.wf_env_smooth) * self.wf_scale
                         for d, sr in zip(data, srs)]
        # Absolute value of each trace if we just want to look at rel. amplitude
        if self.wf_abs:
            data = [np.abs(d) for d in data]
            labels = [f"abs {l}" for l in labels]
        # Scale the amplitudes of each of the data arrays by some value
        data = [d * self.wf_scale for d in data]

        n = len(data)

        # RECSEC/STACK PLOTTING
        # Stack plotting means add the summation of all waveforms as new trace 
        if self.wf_type.endswith("stack"):
            n += 1
            wf_stack_type = "absmean"  # hard-coded for now
            # Determine how to stack
            if wf_stack_type == "mean":
                stacked_data = np.mean(data, axis=0)
                labels.append("stacked mean")
            # Pseudo envelope by meaning all abs values
            elif wf_stack_type == "absmean":
                stacked_data = np.abs(np.mean(np.abs(data), axis=0))
                labels.append("abs stacked mean")

            # Add stacked data to the list of plottable data
            data = np.vstack((data, stacked_data))
            self.axs += [self.ax]
        # Record section means we offset each new data strean by a constant val
        if self.wf_type.startswith("recsec"):
            logger.info("waveform option `recsec`, scaling y-axis")
            if self.wf_spacing_exact:
                # Directly set the spacing with a vlaue
                spacing = self.wf_recsec_spacing
            else:
                # Set an arbitrary spacing based on the order of magnitude 
                oom = np.floor(np.log10(np.abs(self.st[0].data).max()))
                spacing = self.wf_recsec_spacing* 10 ** oom
            logger.info(f"recsec spacing is {spacing}")

            # Visually stack the waveform data by adding offset values to data
            for i, d in enumerate(data[:]):
                data[i] = d + (i * spacing)

        # WAVEFORM STYLES
        # Determine the waveform styling including colors and opacity (alpha)
        alphas = self.alphas or [1] * n
        # User-selection of colors. A single value can be 'C?' (color cycle),
        # any Matplotlib color (e.g. 'k', 'red', '#d62728'), or a colormap
        # name, which is sampled into `n` discrete colors
        _colors = self.colors
        if isinstance(_colors, str):
            _colors = [_colors]
        if len(_colors) > 1:
            colors = list(_colors)
            if len(colors) < n:
                logger.warning(f"{len(colors)} colors for {n} traces, "
                               f"repeating colors")
                colors = [colors[i % len(colors)] for i in range(n)]
        elif _colors[0] in ("C?", "CN"):
            colors = [f"C{i}" for i in range(n)]
        elif mpl.colors.is_color_like(_colors[0]):
            colors = [_colors[0]] * n
        elif _colors[0] in mpl.colormaps:
            cmap = mpl.colormaps[_colors[0]]
            colors = cmap(np.linspace(0.1, 1, n))
        else:
            raise ValueError(f"`colors` value {_colors[0]!r} is not a "
                             f"Matplotlib color or colormap name")
        # User-overwrite of the waveform labels put in legends
        if self.labels:
            labels = self.labels 

        # Set internal attributes for plotting functions
        self._xvals = xvals
        self._data = data
        self._labels = labels 
        self._colors = colors
        self._alphas = alphas

    def plot_waveforms(self):
        """Plot the waveforms based on plot parameters defined before"""
        # Grid style plot or not
        if self.axs:
            axs = self.axs
        else:
            axs = [self.ax] * len(self._data)

        for i, (data, ax, c, l, a) in enumerate(
            zip(self._data, axs,  self._colors, self._labels, self._alphas)):
                # Create time axis on the fly based on the axis chosen but 
                # assuming all data share the same start and end times
                # ax.plot(self._xvals, data,  c=c,  lw=self.linewidth,
                #         zorder=6+i, label=l, alpha=a
                #         )
                ax.plot(self._xvals, data,  c=c,  lw=self.linewidth,
                    zorder=6+i, alpha=a, label=l, 
                    )
                if self.wf_env:
                    if self.wf_env_color:
                        env_c = self.wf_env_color
                    else:
                        env_c = "C0" if c in ["k", "black"] else c
                    _l = l.removeprefix("abs ") if l else ""
                    env_l = f"{_l} envelope" if _l else "envelope"
                    if self.wf_env_smooth:
                        env_l += f" ({self.wf_env_smooth:g} s mean)"
                    ax.plot(self._xvals, self._env[i], c=env_c, lw=1.5,
                            zorder=6 + len(self._data) + i, label=env_l)

    def plot_additional_traces(self):
        """
        Plot additional time series with a common X axis but with a different Y 
        axis. This can be used to plot other geophysical or physical data.
        Required parameters would be `add_trace` which should be in an ObsPy 
        readable format. Spectrogram of this other data can be plotted, too.

        Does not work with ``wf_type='grid'`` option
        """
        twax = self.ax.twinx()
        st_new = read(self.add_trace).merge()
        data_new = st_new[0].data

        # Shift to correct time axis
        if self.time.startswith("a") and len(self.time) > 1:
            st_new = convert_timezone(code=self.time[1:], st=st_new)

        # Custom time axis for the new X-axis so it can still be plotted with
        if self.tr_time.startswith("a"):
            # Set xvalues to datetime objects
            xvals_new = ((st_new[0].times() / SECONDS_PER_DAY) +
                            date2num(st_new[0].stats.starttime.datetime))
            _set_xaxis_obspy_dates(twax, maxticks=self.maxticks)
        else:
            xvals_new = st_new[0].times() 

            # Set time axis
            if self.tr_time == "s":
                xvals_new /= 1  # not necessary but for consistency
            elif self.tr_time == "m":
                xvals_new /= 60
            elif self.tr_time == "h": 
                xvals_new /= 60 ** 2
            else:
                logger.warning("unknown time axis choice, default to 's'econds")
                xvals_new /= 1

            # Offset time axis based on user defined criteria
            xvals_new -= self.t0
            xvals_new += self.tstart

        # Plot on the same axis as the waveform
        if self.tr_label:
            tr_label = self.tr_label
        else:
            tr_label = st_new[0].get_id()
        twax.plot(xvals_new, data_new, lw=1, c=self.tr_color, 
                  label=tr_label, zorder=5, markersize=1.25, alpha=0.8)

        twax.set_ylabel(self.tr_ylabel, rotation=-90, labelpad=20)
        if self.tr_ylim:
            twax.set_ylim(self.tr_ylim)

        self.twax = twax
        # Draw the waveform axis above the twin axis so the added trace sits
        # behind the waveform; the waveform axis background is hidden so the
        # added trace still shows through. `tr_front` keeps Matplotlib's
        # default, where the twin axis is drawn on top
        if not self.tr_front:
            self.ax.set_zorder(twax.get_zorder() + 1)
            self.ax.patch.set_visible(False)
        if self.spectrogram:
            _size = 2.5
            _pad = 0.05
            div = make_axes_locatable(twax)
            cax = div.append_axes("right", size=f"{_size}%", pad=_pad)
            cax.set_axis_off()

    def plot_spectrogram(self):    
        """
        Plot spectrogram on a separate axis above the main waveform plot. Plots
        the spectrogram for the `sp_idx` index of the Stream. Cannot be used
        in conjunction with ``wf_type='grid'`` option
        """
        # Only the displayed range is computed, not the trimming buffers
        xrange = None
        if self.xlim:
            if self.time.startswith("a"):
                xrange = [date2num(UTCDateTime(_).datetime)
                          for _ in self.xlim]
            else:
                xrange = [float(_) for _ in self.xlim]
        tr = self.st[self.sp_idx]
        im = spectrogram(
            ax=self.ax_spectra, xvals=self._xvals, data=tr.data,
            samp_rate=tr.stats.sampling_rate, xrange=xrange,
            wlen=self.sp_wlen, per_lap=self.sp_overlap,
            log=self.sp_logscale, dbscale=self.sp_dbscale,
            cmap=self.sp_cmap, ncolors=self.sp_numcol,
            clip=self.sp_clip, vmax=self.sp_vmax, vmin=self.sp_vmin,
            interpolation=self.sp_interp,
            )
        self.ax_spectra.set_ylabel("Freq. [Hz]")
        self.ax_spectra.axis("tight")
        self.ax_spectra.grid(False) 

        # SPECTROGRAM COLORBAR
        #  Pushes over the spectrogram in the same axis
        _size = 2.5
        _pad = 0.05
        div_spectra = make_axes_locatable(self.ax_spectra)
        cax_spectra = div_spectra.append_axes("right", size=f"{_size}%",
                                                pad=_pad)
        cbar = self.f.colorbar(im, cax=cax_spectra, orientation="vertical")
        # Change the plot aeshtetic of the colorbar, should actually be in
        # set_plot_aesthetic() but I'm lazy
        for spine in cbar.ax.spines.values():
            spine.set_linewidth(1.5)
        cbar.ax.set_frame_on(True)
        # Padding was determined by trial and error
        if self.sp_dbscale:
            # _label = r"Rel. PSD [$10\log_{10}((m/s^2)/Hz)$]"
            _label = f"Power [dB]"
            _labelpad = -37.5
        else:
            _label = r"PSD [$(m/s^2/Hz)$]"
            _labelpad = -25
        cbar.ax.set_ylabel(_label, rotation=270, labelpad=_labelpad,
                           fontsize=FONTSIZES["spectrogram_colorbar_label"])

        # HACKY: pushes over the waveform plot over by the same amount but '
        # then turns the space invisible to preserve the shared x-axis
        div = make_axes_locatable(self.ax)
        cax = div.append_axes("right", size=f"{_size}%", pad=_pad)
        cax.set_axis_off()

    def plot_taup_arrivals(self):
        """
        Get and plot TauP arrivals based on theoretical model, event distance
        and depth. Phases can be input. Relevant parameters are:

        Parameters
        ----------
        - tp_phases (list of str): phases that must match TauP phase list
        - tp_model (str): TauP model to be used
        - tp_dist_km (float): Hypo distance in km (see also `tp_dist_deg`). 
            Choose one of `tp_dist_km` or `tp_dist_deg`
        - tp_dist_deg (float): Hypo distance in deg (see also `tp_dist_km`)
        - tp_depth (float): Depth of event in km
        - tp_start (float): If time=='a', the actual origin time of the event
            to match the data. TauP arrivals are given in relative timing
        """
        # Get phase arrivals from TauP if requested
        arrivals = None
        assert(self.tp_dist_km is not None or self.tp_dist_deg is not None)
        assert(self.tp_depth is not None)
        if self.tp_dist_km:
            dist_deg = kilometers2degrees(self.tp_dist_km)
            self.tp_dist_deg = dist_deg
        else:
            dist_deg = self.tp_dist_deg

        model = TauPyModel(model=self.tp_model)
        tp_arrivals = model.get_travel_times(
            source_depth_in_km=self.tp_depth,
            distance_in_degree=dist_deg,
            phase_list=self.tp_phases
            )

        # If some arrivals have multiple entires, only take first and last to 
        # get a range which we will plot with a window
        arrivals = {arrival.name: [] for arrival in tp_arrivals}
        logger.info(f"TauP Arrivals for {self.tp_model}")
        for i, arrival in enumerate(tp_arrivals):
            arrivals[arrival.name].append(arrival.time)
            logger.info(f"{arrival.name} = {arrival.time:.2f} s")

        if not arrivals:
            logger.warning(f"No arrivals found for given "
                           f"depth={self.tp_depth}km and distance "
                           f"{dist_deg:.2f}deg")
           
        # Determine the time series starttime of the event because the TauP
        # times are relative to an origin time
        if self.tp_start is None:
            tp_start = self._default_origin()
        else:
            if self.time.startswith("a"):
                tp_start = UTCDateTime(self.tp_start)
            else:
                tp_start = float(self.tp_start)  # seconds

        # Offset by USER_T0 if SPECFEM synthetics
        tp_start += self.t0

        # Discrete colormap for each unique arrivals
        cvals = cmaphex(nvals=len(arrivals), cmap=self.tp_cmap)

        plotted_names = []
        for i, (name, times) in enumerate(arrivals.items()):
            # Convert arrival times to time series reference
            if self.time.startswith("a"):
                times = [date2num((tp_start + time).datetime) for time in times]
            else:
                times = [tp_start + time for time in times]

            # Plot a vertical line for each of the arrivals and add to legend
            for ax_i, ax in enumerate(self.axs):
                ylim = self._ylims[ax_i] or self.ylim
                for time in times:
                    ax.axvline(time, alpha=1,  ls="-", color=cvals[i],
                               zorder=100)

                    # TauP often returns multiple arrivals for the same phase
                    # name, only label the first one, the colors will match tho
                    if name not in plotted_names:
                        ax.text(x=time, y=ylim[0]*0.8, s=name,
                                c="k", zorder=125,
                                fontsize=FONTSIZES["taup_phase_label"])
                        plotted_names.append(name)

        # Plot the global ray paths for sanity check
        # https://github.com/obspy/obspy/blob/master/obspy/taup/tau.py#L250
        if self.tp_plot_rays:
            tp_rays = model.get_ray_paths(
                source_depth_in_km=self.tp_depth,
                distance_in_degree=dist_deg,
                phase_list=self.tp_phases,
                )
            # Modify types here
            tp_rays.plot_rays(legend=True, indicate_wave_type=True,
                              plot_type="cartesian", show=False)

    def plot_tmarks(self):
        """
        Plot time marks if requested
        """
        if len(self.tmarks_c) == 1:
            colors = self.tmarks_c * len(self.tmarks)
        else:
            colors = self.tmarks_c

        for tmark, c in zip(self.tmarks, colors):
            # If time is absolute, tmarks are given as UTCDateTime strings
            if self.time.startswith("a"):
                tmark = date2num(UTCDateTime(tmark).datetime)
                # date2num(st[0].stats.starttime.datetime) + tmark
            # Else time is given as a floating point value in s/m/h
            else:
                tmark = float(tmark)

            for ax in self.axs:
                ax.axvline(tmark, c=c, lw=1, zorder=100, 
                           ls="--", alpha=0.8)

    def plot_group_vels(self, c="r", alpha=1, ls="-", lw=2.,
                        fontsize=FONTSIZES["group_velocity_label"],
                        zorder=100):
        """
        Create time mark lines for a given set of group velocities `group_vel`
        given in units of km/s
        Requires values from TauP arrivals, namely: `tp_dist` and `tp_start`
        to get the correct arrival time
        Annotates the value of the group velocity at the max amplitude of the
        waveform.
        """
        # Get phase arrivals from TauP if requested
        assert(self.tp_dist_km is not None or self.tp_dist_deg is not None)
        if self.tp_dist_deg:
            self.tp_dist_km = degrees2kilometers(self.tp_dist_deg)

        # e.g., {5 km/s: 10s}
        arrivals = {str(vel): self.tp_dist_km / vel for vel in self.group_vels}

        # Determine the time series starttime of the event because the TauP
        # times are relative to an origin time
        if self.tp_start is None:
            tp_start = self._default_origin()
        else:
            if self.time.startswith("a"):
                tp_start = UTCDateTime(self.tp_start)
            else:
                tp_start = float(self.tp_start)  # seconds

        # Offset by USER_T0 if SPECFEM synthetics
        tp_start += self.t0

        for i, (name, time) in enumerate(arrivals.items()):
            # Convert arrival times to time series reference
            if self.time.startswith("a"):
                time = date2num((tp_start + time).datetime) 
            else:
                time += tp_start 

            for ax_i, ax in enumerate(self.axs):
                ylim = self._ylims[ax_i] or self.ylim
                # Plot discrete arrivals
                ax.axvline(time, ymin=0, ymax=0.05, alpha=alpha, ls=ls, lw=lw,
                           c=c, zorder=zorder)
                ax.text(x=time, y=self.st[0].max(), s=f"{name}km/s",
                             fontsize=fontsize, color=c, alpha=alpha)
                    
    def set_ylim(self):
        """Establish ylims either by user set or from the natural plot lims"""
        # Establish ylim as [ymin, ymax]
        user_ylim = bool(self.ylim)
        if self.ylim:
            if len(self.ylim) == 1:
                ymax = abs(self.ylim[0])
                ymin = -1 * ymax
                self.ylim = [ymin, ymax]
        # Use natural plot limits
        else:
            self.ylim = plt.ylim()

        # Allow a different ylim for each subplot when using grid-style plots
        if self.wf_type == "grid" and self.ylim_grid:
            assert(len(self.ylim_grid) == 2 * len(self.axs)), (
                f"`ylim_grid` must provide a [ymin, ymax] pair for each of "
                f"the {len(self.axs)} grid subplots, i.e., "
                f"{2 * len(self.axs)} values, not {len(self.ylim_grid)}"
                )
            self._ylims = [self.ylim_grid[2 * i: 2 * i + 2]
                           for i in range(len(self.axs))]
        elif self.wf_type == "grid" and not user_ylim:
            # Each subplot keeps its own natural limits
            self._ylims = [ax.get_ylim() for ax in self.axs]
        else:
            self._ylims = [self.ylim] * len(self.axs)

    def set_pplot_aesthetics(self):
        """
        Set plot aesthetics to give a standard look to all the waveform figures
        """
        # For grid plots, only the bottom-most subplot in each column should
        # carry an x-axis label
        if self.wf_type == "grid":
            bottom_axes = set()
            for col in range(self.ncols):
                col_idxs = [i for i in range(len(self.axs))
                            if i % self.ncols == col]
                if col_idxs:
                    bottom_axes.add(col_idxs[-1])
        else:
            bottom_axes = set(range(len(self.axs)))

        for i, ax in enumerate(self.axs):
            if self.legend:
                handles, labels = ax.get_legend_handles_labels()
                # Include the added trace (`add_trace`) in the same legend
                if self.twax is not None and ax is self.ax:
                    _h, _l = self.twax.get_legend_handles_labels()
                    handles, labels = handles + _h, labels + _l
                leg = ax.legend(handles, labels,
                                loc="upper right", bbox_to_anchor=(1,1),
                                bbox_transform=ax.transAxes,
                                prop={"size": FONTSIZES["legend"]},
                                ncol=self.ncol_legend)
                leg.set(zorder=100)

            # Grid subplots share the x-axis, so only the bottom one in each
            # column needs tick labels
            if i not in bottom_axes:
                ax.tick_params(axis="x", labelbottom=False)
            if i in bottom_axes:
                if self.time.startswith("a"):
                    ax.set_xlabel(f"Time [UTC{self.time[1:]}]")
                else:
                    ax.set_xlabel(f"Time [{self.time}]")
            if self.ylabel:
                ax.set_ylabel(self.ylabel)
            else:
                ax.set_ylabel("Amplitude")

            # Subset x axis
            xstart, xend = self._xvals.min(), self._xvals.max()
            if self.xlim:
                if self.time.startswith("a"):
                    xstart = date2num(UTCDateTime(self.xlim[0]).datetime)
                    xend = date2num(UTCDateTime(self.xlim[1]).datetime)
                else:
                    xstart, xend = [float(_) for _ in self.xlim]
            ax.set_xlim(xstart, xend)  

            # Subset y axis for waveform plot
            if self._ylims[i]:
                ax.set_ylim(self._ylims[i][0], self._ylims[i][1])

            # Spectrogram annotation if there are multiple waveforms plotted
            # put it outside the axis so it shows up on white background
            if self.spectrogram and len(self.st) > 1:
                self.ax.text(0.99, 1.01, self.st[self.sp_idx].get_id(),
                                horizontalalignment="right",
                                verticalalignment="bottom",
                                transform=self.ax.transAxes,
                                fontsize=FONTSIZES["spectrogram_trace_id"])

            # Finish off by setting plot aesthetics
            set_plot_aesthetic(ax, ytick_labels=False)

        if self.spectrogram:
            set_plot_aesthetic(self.ax_spectra, ytick_format="plain")

        if self.title is None:
            title = f"{self.st[0].stats.starttime.year}."
            title += f"{self.st[0].stats.starttime.julday:0>3}"
            if self.st[0].stats.starttime.julday != \
                self.st[0].stats.endtime.julday:
                    title += f"-{self.st[0].stats.endtime.julday}"
            if self.fmin or self.fmax:
                _fmin = self.fmin or 0
                _fmax = self.fmax or self.st[0].stats.sampling_rate / 2
                title += f" [{_fmin}, {_fmax}]Hz"
            # Append some information on the TauP arrivals
            if self.tp_phases:
                title += (f"\n(TauP={self.tp_model}; $\\Delta$="
                          f"{self.tp_dist_deg:.2f}deg; Z={self.tp_depth}km)")
            if self.title_append:
                title += f"\n{self.title_append}"
            if self.title_prepend:
                title = f"{self.title_prepend}{title}"
        else:
            title = self.title

        # With a spectrogram the title goes on the top (spectrogram) axis
        # rather than between the spectrogram and the waveform
        title_ax = self.ax_spectra if self.spectrogram else self.ax
        title_ax.set_title(title, fontsize=FONTSIZES["title"])
        # plt.suptitle(title)

        # Brute force turn off everything
        if self.frameoff:
            self.ax.set_axis_off()
            self.ax.get_legend().remove()
            title_ax.set_title("")

        plt.tight_layout()

    def annotate_subplot_label(self):
        """
        Annotate subplot label. Only works in a select amount of cases
        """
        # Get xlim and ylim from the axis, assuming everything has been set
        plt.sca(self.axs[0])
        plt.figtext(0.001, .99, s=self.subplot_label, c="k",
                    fontsize=FONTSIZES["subplot_label"],
                    zorder=250, ha="left", va="top")

    def finalize(self):
        """
        Final plot adjustments
        """
        if self.save:
            if self.save == "auto":
                _fid_out = f"{self.fids[0]}.png"
            else:
                _fid_out = self.save
            logger.info(f"saving to {_fid_out}")
            plt.savefig(_fid_out, transparent=self.transparent)
            
            # bbox_inches='tight' screws with the actual figure dimensions
                        #bbox_inches="tight")

        if self.show:
            plt.show()
        plt.close("all")

        if self.output:
            if not os.path.exists(self.output):
                os.makedirs(self.output)
            for tr in self.st:
                fid_out = os.path.join(self.output, f"{tr.get_id()}.mseed")
                logger.info(f"writing {fid_out}")
                tr.write(fid_out, format="MSEED")

    def run(self):
        """
        Main plotting function
        """
        self.setup_plot()
        self.trim_waveform()
        if self.inv:
            self.remove_response()
        self.process_waveforms()
        self.get_plot_parameters()
        self.plot_waveforms()
        self.set_ylim()
        if self.add_trace:
            self.plot_additional_traces()
        if self.spectrogram:
            self.plot_spectrogram()
        if self.tmarks:
            self.plot_tmarks()
        if self.group_vels:
            self.plot_group_vels()
        if self.tp_phases:
            self.plot_taup_arrivals()
        if self.subplot_label:
            self.annotate_subplot_label()
        self.set_pplot_aesthetics()
        self.finalize()


def main():
    """
    Command line entry point for ``pplot``: parse arguments, convert the
    negative (``*_off``, ``no*``) flags to keyword arguments and run
    PrettyPlot
    """
    args = parse_args()
    args_dict = vars(args).copy()
    fids = args_dict.pop("fids")
    args_dict["legend"] = not args_dict.pop("no_legend", False)
    args_dict["show"] = not args_dict.pop("noshow", False)
    args_dict["sp_dbscale"] = not args_dict.pop("sp_dbscale_off", False)
    args_dict["rr_zero_mean"] = not args_dict.pop("rr_zero_mean_off", False)
    args_dict["rr_taper"] = not args_dict.pop("rr_taper_off", False)
    pp = PrettyPlot(fids, **args_dict)
    pp.run()


if __name__ == "__main__":
    main()
