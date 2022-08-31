#!/usr/bin/evn python3
"""
RECord SECtion plotting tool for seismic waveforms (observed and synthetic)

This is a refactor of Pysep's Python utility `plotw_rs`, a record section
plotting script. The intent of this script is to plot multiple time series'
based on source-receiver characteristics (i.e., src-rcv distance, backazimuth).

.. note:: Code History
    Written by Carl Tape (11/21/2011) and Yun Wang (11/2011) in Matlab
    Translated to Python by Nealy Sims (1/2021)
    Upgraded by Aakash Gupta (9/2021)
    Refactored by Bryant Chow (3/2022)

.. requires::
    obspy >= 1.2 (expected to bring in numpy and matplotlib)

.. rubric:: Examples
    1) Print the help message to see available options and flags
    $ python recsec.py -h

    2) From the command line: The following example code blocks work with pysep 
        to download data for a southern California event.

        # cd path/to/pysep
        $ python rungetwaveform.py event_input_mtuq2022 2  # 20200404015318920

    a) Plot a record section for the socal event with default values

        $ python recsec.py --pysep_path 20200404015318920 

    b) Plot high-passed data with 7 km/s move out (focused on direct arrivals),
        show direct arrivals through to surface waves for all traces, thin lines
        to accentuate high frequencies. Overwrite previous figure and split
        figure onto multiple pages

        $ python recsec.py --pysep_path 20200404015318920 \
            --move_out 7 --min_period_s 1 --xlim_s 100 175 \
            --linewidth .25 --max_traces_per_rs 60 --overwrite

    c) Plot bandpassed data with 4 km/s move out (focused on surface waves),
        horizontal components only (radial, transverse),
        thicken up default linewidth and increase spacing between adjacent 
        seismograms 

        $ python recsec.py --pysep_path 20200404015318920 \
            --components RT --move_out 4 --min_period_s 2 --max_period_s 50 \
            --xlim_s 50 200 --y_axis_spacing 3 --linewidth 1 \
            --amplitude_scale_factor 4 --overwrite 

    d) Plot bandpassed transverse component, sort by azimuth, scale by the 
        maximum amplitude of ALL traces shown on the figure. 
        Scale amplitudes by factor 2 for better visualization
        and start azimuth plotting at 180* (as opposed to default 0*).

        $ python recsec.py --pysep_path 20200404015318920 \
            --sort_by azimuth --scale_by global_norm --components T \
            --min_period_s 2 --max_period_s 30 --move_out 4 \
            --amplitude_scale_factor 2 --azimuth_start_deg 180 \
            --linewidth 1 --overwrite

    e) Plot bandpassed vertical components sorted by absolute distance.
        Reduce amplitudes by 1/4 (0.25). Set y label fontsize smaller than 
        default and at the max X value of the figure (far to the right)

        $ python recsec.py --pysep_path 20200404015318920 \
            --sort_by abs_distance_r --components Z --min_period_s 2 \
            --max_period_s 50 --amplitude_scale_factor 0.25 \
            --y_label_loc x_max --y_label_fontsize 7 --overwrite \
            --linewidth 1

    3) From the Python interpreter: this is an example code block that can be
        written into a Python script or run from inside a Python interactive 
        environment. Code block assumes we have downloaded event data with Pysep

        >>> import os
        >>> from glob import glob
        >>> from obspy import read
        >>> from recsec import plotw_rs
        >>> st = Stream()
        >>> for fid in glob(os.path.join("20200404015318920", "*.?")):
        >>>     st += read(fid)
        >>> plotw_rs(st=st, sort_by="distance_r")
"""
import os
import sys
import argparse
import textwrap
import numpy as np
import matplotlib.pyplot as plt

from glob import glob
from datetime import datetime
from matplotlib.ticker import MultipleLocator
from obspy import read, Stream
from obspy.geodetics import (kilometers2degrees, gps2dist_azimuth)

from pysep import logger
from pysep.utils.io import read_synthetics

# Unicode degree symbol for plot text
DEG = u"\N{DEGREE SIGN}"


def myround(x, base=5, choice="near"):
    """
    Round value x to nearest base, round 'up','down' or to 'near'est base

    :type x: float
    :param x: value to be rounded
    :type base: int
    :param base: nearest integer to be rounded to
    :type choice: str
    :param choice: method of rounding, 'up', 'down' or 'near'
    :rtype roundout: int
    :return: rounded value
    """
    if choice == "near":
        roundout = int(base * round(float(x)/base))
    elif choice == "down":
        roundout = int(base * np.floor(float(x)/base))
    elif choice == "up":
        roundout = int(base * np.ceil(float(x)/base))
    return roundout


class Dict(dict):
    """Simple dictionary overload for nicer get/set attribute characteristics"""
    def __setattr__(self, key, value):
        self[key] = value

    def __getattr__(self, key):
        return self[key]


class RecordSection:
    """
    Record section plotting tool which takes ObsPy streams and 
    1) preprocesses and filters waveforms,
    2) sorts source-receiver pairs based on User input, 
    3) produces record section waveform figures.
    """
    def __init__(self, pysep_path=None, syn_path=None, stations=None,
                 cmtsolution=None, st=None, st_syn=None, sort_by="default",
                 scale_by=None, time_shift_s=None, zero_pad_s=None,
                 move_out=None, min_period_s=None, max_period_s=None,
                 preprocess="st", max_traces_per_rs=None, integrate=0,
                 xlim_s=None, components="ZRTNE12", y_label_loc="default",
                 y_axis_spacing=1, amplitude_scale_factor=1,
                 azimuth_start_deg=0., distance_units="km", 
                 geometric_spreading_factor=0.5, geometric_spreading_k_val=None, 
                 figsize=(9, 11), show=True, save="./record_section.png", 
                 overwrite=False, **kwargs):
        """
        Set the default record section plotting parameters and enforce types.
        Run some internal parameter derivation functions by manipulating input
        data and parameters.

        :type pysep_path: str
        :param pysep_path: path to Pysep output, which is expected to contain
            trace-wise SAC waveform files which will be read
        :type syn_path: str
        :param syn_path: path to SPECFEM generated ASCII synthetics.
        :type st: obspy.core.stream.Stream
        :param st: Stream objects containing observed time series to be plotted
            on the record section. Can contain any number of traces
        :type st_syn: obspy.core.stream.Stream
        :param st_syn: Stream objects containing synthetic time series to be
            plotted on the record section. Must be same length as `st`
        :type sort_by: str
        :param sort_by: How to sort the Y-axis of the record section, available:
            - 'default': Don't sort, just iterate directly through Stream
            - 'alphabetical': sort alphabetically A->Z. Components sorted 
                separately with parameter `components`
            - 'azimuth': sort by source-receiver azimuth (deg) with constant
                vertical spacing on the y-axis. Requires `azimuth_start_deg`
            - 'backazimuth': sort by source-receiver backazimuth (deg) with
                constant vertical spacing. Requires `azimuth_start_deg`
            - 'distance': sort by source-receiver distance (km) with constant
                vertical spacing. Requires `distance_units`
            - 'abs_distance': absolute vertical spacing of waveforms defined by
                source-receiver distance. Requires `distance_units`
            - 'abs_azimuth': absolute vertical spacing of waveforms defined
                by source-receiver azimuth (deg).
            - 'abs_backazimuth': absolute vertical spacing of waveforms by
                source-receiver backazimuth (deg).
            - '*_r': Add a '_r' to any of the values about to REVERSE the sort,
                e.g., 'alphabetical_r' sort will go Z->A
        :type scale_by: list
        :param scale_by: scale amplitude of waveforms by available:
            - None: Not set, no amplitude scaling, waveforms shown raw
            - 'normalize': scale each trace by the maximum amplitude,
                i.e., > a /= max(abs(a))  # where 'a' is time series amplitudes
            - 'global_norm': scale by the largest amplitude to be displayed on
                the screen. Will not consider waveforms which have been 
                excluded on other basis (e.g., wrong component)
                i.e., > st[i].max /= max([max(abs(tr.data)) for tr in st])
            - 'geometric_spreading': TODO this is not implemented yet.
                scale amplitudes globally by predicting the
                expected geometric spreading amplitude reduction and correcting
                for this factor. Requires `geometric_spreading_factor`, and 
                optional `geometric_spreading_k_val`
        :type time_shift_s: float OR list of float
        :param time_shift_s: apply static time shift to waveforms, two options:
            1. float (e.g., -10.2), will shift ALL waveforms by
                that number (i.e., -10.2 second time shift applied)
            2. list (e.g., [5., -2., ... 11.2]), will apply individual time
                shifts to EACH trace in the stream. The length of this list MUST
                match the number of traces in your input stream.
        :type zero_pad_s: list
        :param zero_pad_s: zero pad data in units of seconsd. applied after
            tapering and before filtering. Input as a tuple of floats,
            * (start, end): a list of floats will zero-pad the START and END
                of the trace. Either can be 0 to ignore zero-padding at either
                end
        :type amplitude_scale_factor: float OR list of float
        :param amplitude_scale_factor: apply scale factor to all amplitudes.
            Used as a dial to adjust amplitudes manually. Defaults to 1.
            Two options:
            1. float (e.g., 1.2), will multiply ALL waveforms by that number
            2. list (e.g., [5., -2., ... 11.2]), will apply individual amplitude
                scale to EACH trace in the stream. The length of this list MUST
                match the number of traces in your input stream.
        :type move_out: float
        :param move_out: Optional. A velocity value that will be used to
            calculate move out, which will time shift seismograms based on
            their source receiver distance. This parameter will be ADDED
            to time_shift_s (both float and list), if it is provided.
            Should be in units of `distance_units`/s
        :type min_period_s: float
        :param min_period_s: minimum filter period in seconds
        :type max_period_s: float
        :param max_period_s: maximum filter period in seconds
        :type preprocess: str
        :param preprocess: choose which data to preprocess, options are:
            - 'st': process waveforms in st (default)
            - 'st_syn': process waveforms in st_syn. st still must be given
            - 'both': process waveforms in both st and st_syn
            - None: do not run preprocessing step (including filter)
        :type max_traces_per_rs: int
        :param max_traces_per_rs: maximum number of traces to show on a single
            record section plot. Defaults to all traces in the Stream
        :type xlim_s: list of float
        :param xlim_s: [start, stop] in units of time, seconds, to set the
            xlimits of the figure
        :type components: str
        :param components: a sequence of strings representing acceptable
            components from the data. Also determines the order these are shown
            EVEN when sorted by other variables. For example, components=='ZR'
            would only display Z and R components, and Z components would be
            should BEFORE R components for the SAME station.
        :type integrate: int
        :param integrate: apply integration `integrate` times on all traces.
            acceptable values [-inf, inf], where positive values are integration
            and negative values are differentiation
            e.g., if integrate == 2,  will integrate each trace twice.
            or    if integrate == -1, will differentiate once
            or    if integrate == 0,  do nothing (default)
        :type y_axis_spacing: float
        :param y_axis_spacing: spacing between adjacent seismograms applied to
            Y-axis on relative (not absolute) scales. Defaults to 1.
        :type y_label_loc: str
        :param y_label_loc: Location to place waveform labels on the y-axis
            - 'default': auto choose the best location based on `sort_by`
            - 'y_axis': Replace tick labels on the y-axis (left side of figure),
                This won't work if using absolute sorting and will be over-
                written by 'default'
            - 'y_axis_right': Replace tick labels on the right side of the 
                y-axis. This option won't work with absolute sorting
            - 'x_min': Place labels on top of the waveforms at the global min
                x-value on the figure
            - 'x_min': Place labels on top of the waveforms at the global max
                x-value on the figure
            - None: Don't plot any text labels
        :type azimuth_start_deg: float
        :param azimuth_start_deg: If sorting by azimuth, this defines the 
            azimuthal angle for the waveform at the top of the figure.
            Set to 0 for default behavior
        :type distance_units: str
        :param distance_units: Y-axis units when sorting by epicentral distance
            'km': kilometers on the sphere
            'deg': degrees on the sphere
            'km_utm': kilometers on flat plane, UTM coordinate system
        :type geometric_spreading_factor: float
        :param geometric_spreading_factor: geometrical spreading factor when
            using the `scale_by` parameter. Defaults to 0.5 for surface waves.
            Use values of 0.5 to 1.0 for regional surface waves
        :type geometric_spreading_k_val: float
        :param geometric_spreading_k_val: K value used to scale the geometric
            spreading factor (TODO figure out what this actually is)
        :type figsize: tuple of float
        :param figsize: size the of the figure, passed into plt.subplots()
        :type show: bool
        :param show: show the figure as a graphical output
        :type save: str
        :param save: path to save output figure, will create the parent
            directory if it doesn't exist. If None, will not save.
        :type overwrite: bool
        :param overwrite: if the path defined by `save` exists, will overwrite
            the existing figure
        :raises AssertionError: if any parameters are set incorrectly
        """
        # Read files from path if provided
        if pysep_path is not None and os.path.exists(pysep_path):
            # Expecting to find SAC files labelled as such
            fids = glob(os.path.join(pysep_path, "*.sac"))
            if not fids:
                # Or if legacy naming schema, assume that the sac files have a
                # given file format: event_tag.NN.SSS..LL.CC.c
                fids = glob(os.path.join(pysep_path, "*.*.*.*.*.?"))
            if fids:
                logger.info(f"Reading {len(fids)} files from: {pysep_path}")
                # Overwrite stream, so reading takes precedence
                st = Stream()
                for fid in fids:
                    st += read(fid)

        # Read in SPECFEM generated synthetics and generate SAC headed streams
        if syn_path is not None and os.path.exists(syn_path):
            assert(cmtsolution is not None and os.path.exists(cmtsolution))
            assert(stations is not None and os.path.exists(stations))
            fids = glob(os.path.join(syn_path, "??.*.*.sem?*"))
            if fids:
                logger.info(f"Reading {len(fids)} synthetics from: {syn_path}")
                st_syn = Stream()
                for fid in fids:
                    st_syn += read_synthetics(fid=fid, cmtsolution=cmtsolution,
                                              stations=stations)
        # Allow plotting ONLY synthetics and no data
        if st is None:
            st = st_syn.copy()
            st_syn = None
        assert st, ("Stream object not found, please check inputs `st` "
                    "and `pysep_path")

        # User defined parameters, do some type-setting
        self.st = st.copy()
        try:
            self.st_syn = st_syn.copy()
        except AttributeError:
            self.st_syn = None

        # Y-Axis sorting parameters
        self.sort_by = sort_by.lower()
        self.y_axis_spacing = float(y_axis_spacing)
        self.azimuth_start_deg = float(azimuth_start_deg)
        self.components = str(components)

        # Amplitude scaling parameters
        self.scale_by = scale_by
        self.amplitude_scale_factor = amplitude_scale_factor
        self.geometric_spreading_factor = float(geometric_spreading_factor)
        self.geometric_spreading_k_val = geometric_spreading_k_val

        # Time shift parameters
        self.move_out = move_out
        self.time_shift_s = time_shift_s
        self.zero_pad_s = zero_pad_s

        # Filtering parameters
        self.min_period_s = min_period_s
        self.max_period_s = max_period_s
        self.preprocess = preprocess
        self.max_traces_per_rs = max_traces_per_rs
        self.integrate = int(integrate)

        # Plotting parameters
        self.xlim_s = xlim_s
        self.distance_units = distance_units.lower()
        self.y_label_loc = y_label_loc
        self.figsize = figsize
        self.show = bool(show)
        self.save = save
        self.overwrite = bool(overwrite)
        self.kwargs = Dict(kwargs)

        # Run checks to ensure that all the parameters are set properly
        # And get some new, initial parameters that are required for other
        # 'get' functions
        self.check_parameters()
        self.stats = self.get_srcrcv_stats()
        self.distances, self.azimuths, self.backazimuths = \
            self.get_srcrcv_dist_az_baz()

        # Internally used parameters that will be filled out by other functions
        self.f = None
        self.ax = None
        self.idx = []
        self.station_ids = []
        self.max_amplitudes = []
        self.amplitude_scaling = []
        self.y_axis = []
        self.xlim = []  # unit: samples
        self.sorted_idx = []

    def check_parameters(self):
        """
        Check that parameters are set properly and in line with how they
        are expected by the program

        .. note::
            Not using assertions here because we want all incorrect parameters
            to be evaluated together and displayed at once, that way the user
            doesn't have to run this function multiple times to figure out how
            to set their parameters correctly

        :raises AssertionError: If any parameters are not set as expected by
            plotw_rs functionality
        """
        logger.info("checking parameter acceptability")

        # Used to keep track of which parameters failed and in what way
        err = Dict()

        # Check to make sure there is data
        if not bool(self.st):
            err.st = f"stream has {len(self.st)} traces, no data"

        # Check that stream has SAC headers if we want to sort by dist or (b)az.
        # Pysep should have created these
        if self.sort_by != "default" and \
                any(_ in self.sort_by for _ in ["azimuth", "distance"]):
            _idx = []
            for i, tr in enumerate(self.st):
                if not hasattr(tr.stats, "sac"):
                    _idx.append(i)
            if _idx:
                err.st = (f"{len(_idx)} traces have no SAC header, recsec "
                          f"expects SAC headers for sorting. Trace indexes "
                          f"are: {_idx}")

        if self.st_syn is not None:
            if len(self.st) != len(self.st_syn):
                err.st_syn = f"length must match `st` (which is {len(self.st)})"

        # Check the `sort_by` sorting parameter options
        acceptable_sort_by = ["default", "azimuth", "backazimuth",
                              "distance", "alphabetical", "abs_azimuth",
                              "abs_distance"]
        # Allow reverse sorts
        acceptable_sort_by += [f"{_}_r" for _ in acceptable_sort_by]
        if self.sort_by not in acceptable_sort_by:
            err.sort_by = f"must be in {acceptable_sort_by}"

        if "azimuth" in self.sort_by:
            if not (0 <= self.azimuth_start_deg <= 360):
                err.azimuth_start_deg = f"0 < azi < 360"

        acceptable_distance_units = ["km", "km_utm", "deg"]
        if ("distance" in self.sort_by) and \
                (self.distance_units not in acceptable_distance_units):
            err.azimuth_start_deg = \
                f"must be in {acceptable_distance_units}"

        if self.scale_by is not None:
            acceptable_scale_by = ["normalize", "global_norm",
                                   "geometric_spreading"]
            if self.scale_by not in acceptable_scale_by:
                err.scale_by = f"must be in {acceptable_scale_by}"

        if self.time_shift_s is not None:
            acceptable_time_shift_s = [int, float, list]
            if type(self.time_shift_s) not in acceptable_time_shift_s:
                err.time_shift_s = f"must be in {acceptable_time_shift_s}"
            if isinstance(self.time_shift_s, list) and \
                    len(self.time_shift_s) != len(self.st):
                err.time_shift_s = f"must be list of length {len(self.st)}"

        if self.zero_pad_s is not None:
            assert(isinstance(self.zero_pad_s, list)), (
                f"`zero_pad_s` must be a list of floats representing the " 
                f"[`start`, `end`] amount to pad around trace, units seconds"
            )
            assert(len(self.zero_pad_s) == 2), (
                f"`zero_pad_s` must be a list of floats length 2"
            )
            _start, _end = self.zero_pad_s
            assert(_start >= 0. and _end >= 0.), \
                f"`zero_pad_s` values must be >= 0"

        # Defaults to float value of 1
        acceptable_scale_factor = [int, float, list]
        if type(self.amplitude_scale_factor) not in acceptable_scale_factor:
            err.amplitude_scale_factor = f"must be in {acceptable_scale_factor}"
        if isinstance(self.amplitude_scale_factor, list) and \
                len(self.amplitude_scale_factor) != len(self.st):
            err.amplitude_scale_factor = f"must be list length {len(self.st)}"

        if self.min_period_s is not None and self.max_period_s is not None:
            if self.min_period_s >= self.max_period_s:
                err.min_period_s = "must be less than `max_period_s`"

        if self.preprocess is not None:
            acceptable_preprocess = ["both", "st", "st_syn"]
            if self.preprocess not in acceptable_preprocess:
                err.preprocess = f"must be in {acceptable_preprocess}"
        if self.preprocess in ["st_syn", "both"]:
            assert(self.st is not None and self.st_syn is not None), (
                f"`preprocess` choice requires both `st` & `st_syn` to exist."
                f"If you only have one or the other, set: `preprocess`=='st'"
            )

        # Overwrite the max traces per record section, enforce type int
        if self.max_traces_per_rs is None:
            self.max_traces_per_rs = int(len(self.st))
        else:
            self.max_traces_per_rs = int(self.max_traces_per_rs)

        if self.xlim_s is not None:
            if len(self.xlim_s) != 2:
                err.xlim_s = f"must be of length 2, [start, stop]"
            elif self.xlim_s[0] > self.xlim_s[1]:
                err.xlim_s = f"start time must be less than stop time"

        acceptable_y_label_loc = ["default", "y_axis", "y_axis_right", "x_min",
                                  "x_max", None]
        if self.y_label_loc not in acceptable_y_label_loc:
            err.y_label_loc = f"must be in {acceptable_y_label_loc}"
        if "abs" in self.sort_by and "y_axis" in self.sort_by:
            err.y_label_loc = (f"absolute sorting means 'y_axis' label loc is" 
                               f"not available")

        if len(self.figsize) != 2:
            err.figsize = "must be tuple defining (horizontal, vertical) extent"

        if os.path.exists(self.save) and not self.overwrite:
            err.save = (f"path {self.save} already exists. Use '--overwrite' " 
                        f"flag to save over existing figures.")

        _dirname = os.path.abspath(os.path.dirname(self.save))
        if not os.path.exists(_dirname):
            logger.info(f"creating output directory {_dirname}")
            os.makedirs(_dirname)

        if err:
            out = "ERROR - Parameter errors, please make following changes:\n"
            out += "\n".join([f"`{key}`: {val}" for key, val in err.items()])
            logger.info(out)
            sys.exit(-1)

    def get_skip_idx(self):
        """
        Get a list of any traces that don't adhere to user-defined boundaries
        such as dist, az, baz, id, or component matches. Don't actually remove
        the traces from the stream but rather just collect indices we can use
        to skip when plotting.

        TODO add distance, azi and backazi skip criteria

        :rtype: np.array
        :return: returns an indexing list which can be used to skip over
            traces that don't adhere to certain criteria
        """
        skip_idx = []
        for idx in self.idx:
            tr = self.st[idx]
            # Component-wise removal
            if tr.stats.component not in self.components:
                skip_idx.append(idx)
            # !!! Add more here
        logger.info(f"criteria check will remove "
                    f"{len(skip_idx)}/{len(self.st)} traces")

        return np.array(skip_idx)

    def get_parameters(self):
        """
        Calculate parameters in a specific order and based on the user-defined
        information.

        .. note::
            The order of function calls here is important! Some of the 'get' 
            functions require the results of other 'get' functions. 

        Calculated Parameters
        ::
            np.array idx:
                a linear indexing of all the traces in the stream
            np.array station_ids:
                an ordered list of station ids, used to get station names
                that match the index defined in `idx`
            np.array max_amplitudes:
                abs max amplitudes of each trace, used for normalization
            np.array amplitude_scaling:
                An array to scale amplitudes based on user choices
            np.array time_shift_s:
                An array to time shift time series based on user choices
            np.array y_axis:
                Y-Axis values based on sorting algorithm, used for plotting
            np.array distances:
                source-receiver distances in `distance_units` units
            np.array azimuths:
                source-receiver azimuths in degrees
            np.array backazimuths:
                source-receiver backazimuths in degrees
            np.array sorted_idx:
                sorted indexing on all the traces of the stream based on the
                chosen sorting algorithm
        """
        # Extract basic information from the Stream
        self.idx = np.arange(0, len(self.st), 1)
        self.station_ids = np.array([tr.get_id() for tr in self.st])

        self.time_shift_s = self.get_time_shifts()  # !!! OVERWRITES user input
        self.xlim = self.get_xlims()

        # Max amplitudes should be RELATIVE to what were showing (e.g., if
        # zoomed in on the P-wave, max amplitude should NOT be the surface wave)
        for tr, xlim in zip(self.st, self.xlim):
            start, stop = xlim
            self.max_amplitudes.append(max(abs(tr.data[start:stop])))
        self.max_amplitudes = np.array(self.max_amplitudes)

        # Figure out which indices we'll be plotting
        sorted_idx = self.get_sorted_idx()
        skip_idx = self.get_skip_idx()
        # Remove skip indexes from sorted index to get the final ordered list
        self.sorted_idx = np.array([_ for _ in sorted_idx if _ not in skip_idx])

        # Figure out how to manipulate each of the traces in the Stream
        self.y_axis = self.get_y_axis_positions()
        self.amplitude_scaling = self.get_amplitude_scaling()

    def get_xlims(self):
        """
        The x-limits of each trace depend on the overall time shift (either 
        static or applied through move out), as well as the sampling rate of
        each trace (which can vary). Retrieve an index-dependent list of
        x-limits which can be used to truncate the time series during plotting.

        .. note::
            Requires that get_time_shifts() has already been run

        :rtype: np.array
        :return: an array of tuples defining the start and stop indices for EACH
            trace to be used during plotting. Already includes time shift 
            information so xlim can be applied DIRECTLY to the time shifted data
        """
        xlim = []
        if self.xlim_s is None:
            # None's will index the entire trace
            xlim = np.array([(None, None) for _ in range(len(self.st))])
        else:
            # Looping to allow for delta varying among traces,
            # AND apply the time shift so that indices can be used directly in
            # the plotting function
            for tr, tshift in zip(self.st, self.time_shift_s):
                start, stop = [int(_/tr.stats.delta) for _ in self.xlim_s]
                sshift = int(tshift / tr.stats.delta)  # unit: samples
                xlim.append((start-sshift, stop-sshift))

        xlim = np.array(xlim)
        return xlim

    def get_srcrcv_stats(self):
        """
        Get source receiver information such as min max values, and
        count-related numbers (e.g., num stations) to be used mainly for print
        statements and text information

        Stats Arguments
        ::
            np.array event_names:
                unique event names taken from the SAC header
            int nevents:
                number of unique events in the stream
            np.array unique_sta_ids:
                unique station codes taken from trace stats
            int nstation_ids:
                number of unique station codes
            np.array network_codes:
                unique network codes taken from station ids
            int nnetwork:
                number of unique network codes
            np.array station_codes:
                unique station codes taken from station ids
            int nstation:
                number of unique station codes
            np.array location_codes:
                unique location codes taken from station ids
            int nlocation:
                number of unique location codes
            np.array channel_codes:
                unique channel codes taken from station ids
            int nchannel:
                number of unique channel codes
            bool reverse_sort:
                determine if the user wants to reverse their sort, they do this
                by appending '_r' to the end of the `sort_by` argument
        """
        logger.info("getting source-receiver stats")

        def _unique(list_):
            """return a unique numpy array derived from a list"""
            return np.unique(np.array(list_, dtype=str))

        stats = Dict()
        stats.event_names = _unique([tr.stats.sac.kevnm for tr in self.st])
        stats.nevents = len(stats.event_names)
        stats.unique_sta_ids = _unique([tr.get_id() for tr in self.st])
        stats.longest_id = max([len(_) for _ in stats.unique_sta_ids])
        stats.nstation_ids = len(stats.unique_sta_ids)

        # Get unique network, station, location and channel codes. Also numbers
        for name in ["network", "station", "location", "channel", "component"]:
            stats[f"{name}_codes"] = _unique(
                [getattr(tr.stats, name) for tr in self.st]
            )
            stats[f"n{name}"] = len(stats[f"{name}_codes"])

        # We use `not` in `reverse_sort` to make the top of the y-axis the
        # starting point, which seems more intuitive for record sections, but
        # is opposite the behavior when you increment from 0
        stats.reverse_sort = not bool("_r" in self.sort_by)

        # Initiate empty lists for _plot_trace() to fill with min and max data
        # values which can be used for global plotting parameters like xlims
        stats.xmin, stats.xmax, stats.ymin, stats.ymax = [], [], [], []

        return stats

    def get_time_shifts(self):
        """
        Very simple function which allows float inputs for time shifts and
        ensures that time shifts are always per-trace arrays
        Applies the move out by calculating a time shift using src-rcv distance

        :rtype: np.array
        :return: a stream-lengthed array of time shifts that can be applied
            per trace
        """
        # No user input means time shifts will be 0, so nothing happens
        time_shift_arr = np.zeros(len(self.st))
        if self.time_shift_s is not None:
            # User inputs a static time shift
            if isinstance(self.time_shift_s, (int, float)):
                time_shift_arr += self.time_shift_s
            # User input an array which should have already been checked for len
            else:
                time_shift_arr = self.time_shift_s
        time_shift_arr = np.array(time_shift_arr)

        # Further change the time shift if we have move out input
        if self.move_out:
            logger.info(f"apply {self.move_out} {self.distance_units}/s "
                        f"move out")
            move_out_arr = self.distances / self.move_out
            time_shift_arr -= move_out_arr

        return time_shift_arr

    def get_srcrcv_dist_az_baz(self):
        """
        Convenience function to wrap _get_srcrcv_dist_az_baz_trace into a loop 
        over the whole stream and return lists of distances, azimuths, and 
        backazimuths

        :rtype distances: np.array
        :return distances: source-receiver distances in user-defined units in
            the original order of Stream
        :rtype azimuths: np.array
        :return azimuths: source-receiver azimuths (deg) in the original
            order of Stream
        :rtype backazimuths: np.array
        :return backazimuths: source-receiver azimuths (deg) in the original
            order of Stream
        """
        logger.info("calculating source-receiver distance and (back)azimuths")

        distances, azimuths, backazimuths = [], [], []
        for tr in self.st:
            gcd, az, baz = self._get_srcrcv_dist_az_baz_trace(tr=tr)
            distances.append(gcd)
            azimuths.append(az)
            backazimuths.append(baz)

        distances = np.array(distances)
        azimuths = np.array(azimuths)
        backazimuths = np.array(backazimuths)

        return distances, azimuths, backazimuths

    def _get_srcrcv_dist_az_baz_trace(self, tr=None, idx=0):
        """
        Check the source-receiver characteristics such as src-rcv distance,
        azimuth, backazimuth for a given trace.

        .. note::
            This function ASSUMES that SAC headers have been written to the
            traces. Otherwise we will need more complicated ways to get
            event lat and lon

        :type tr: obspy.core.trace.Trace
        :param tr: trace to get srcrcv information for. If None, will use `idx`
        :type idx: int
        :param idx: if `tr`==None, user can provide index of self.st (Stream)
            defaults to index 0
        :rtype gcdist: float
        :return gcdist: great circle distance in units specified by
            `distance_units`, can be 'km' or 'deg'
        :rtype az: float
        :return az: azimuth (degrees) between source and receiver
        :rtype baz: float
        :return baz: azimuth (degrees) between source and receiver
        """
        # Default to the first trace in the Stream
        if tr is None:
            tr = self.st[int(idx)]
        try:
            dist = tr.stats.sac.dist  # units: km
            az = tr.stats.sac.az        # units: deg
            baz = tr.stats.sac.baz      # units: deg
        except AttributeError:
            # If for whatever reason SAC headers dont contain this info already
            # Use ObsPy to get great circle distance, azimuth and backazimuth
            dist, az, baz = gps2dist_azimuth(lat1=tr.stats.sac.evla,
                                             lon1=tr.stats.sac.evlo,
                                             lat2=tr.stats.sac.stla,
                                             lon2=tr.stats.sac.stlo)
            dist *= 1E-3   # units: m -> km

        # Ensure that 0 <= (b)az <= 360
        az = az % 360
        baz = baz % 360

        # Make sure distance is in the correct units, default units 'km'
        if self.distance_units == "deg":
            dist = kilometers2degrees(dist)  # units: km -> deg
        elif self.distance_units == "km_utm":
            # Overwrite `dist`, could probably skip that calc above but
            # leaving for now as I don't think this option will be used heavily.
            dist_deg = np.sqrt(
                ((tr.stats.sac.stlo - tr.stats.sac.evlo) ** 2) /
                ((tr.stats.sac.stla - tr.stats.sac.evla) ** 2)
            )
            dist = kilometers2degrees(dist_deg)  # units: km

        return dist, az, baz

    def get_amplitude_scaling(self):
        """
        Scale the amplitudes of all the waveforms by producing a Stream
        dependent scale factor based on user choice. It is expected that the
        output array will be DIVIDED by the data arrays:

        i.e., st[i].data /= self.amplitude_scaling[i]

        .. note::
            Needs to be run AFTER preprocessing because filtering etc. will
            affect the final amplitudes of the waveforms

        :rtype: np.array
        :return: an array corresponding to the Stream indexes which provides
            a per-trace scaling coefficient
        """
        logger.info(f"determining amplitude scaling with: {self.scale_by}")

        # Don't scale by anything
        if self.scale_by is None:
            amp_scaling = np.ones(len(self.st))
        # Scale by the max amplitude of each trace
        elif self.scale_by == "normalize":
            amp_scaling = self.max_amplitudes
            # When using absolute distance scale, scale waveforms to minmax dist
            if "abs" in self.sort_by:
                if "distance" in self.sort_by:
                    logger.info("scaling amplitudes for absolute distance")
                    scale = np.mean(self.distances) / 10
                elif "backazimuth" in self.sort_by:
                    logger.info("scaling amplitudes for absolute backazimuth")
                    scale = self.backazimuths.max() - self.backazimuths.min()
                    scale /= 100
                elif "azimuth" in self.sort_by:
                    logger.info("scaling amplitudes for absolute azimuth")
                    scale = self.azimuths.max() - self.azimuths.min()
                    scale /= 50
                # Divide to make waveforms larger
                amp_scaling /= scale
        # Scale by the largest amplitude shown
        elif self.scale_by == "global_norm":
            global_max = self.max_amplitudes[self.sorted_idx].max()
            amp_scaling = np.ones(len(self.st)) * global_max
        # Scale by the theoretical geometrical spreading factor
        elif self.scale_by == "geometric_spreading":
            amp_scaling = self._calculate_geometric_spreading()

        # Apply manual scale factor if provided, default value is 1 so nothing
        # Divide because the amplitude scale divides the data array, which means
        # `amplitude_scale_factor` will be MULTIPLIED by the data array
        amp_scaling /= self.amplitude_scale_factor

        return amp_scaling

    def _calculate_geometric_spreading(self):
        """
        Stations with larger source-receiver distances will have their amplitude
        scaled by a larger value.

        For information on geometric spreading, see Stein and Wysession,
        Section 4.3.4, Eq 20 (for Rayleigh waves, geometrical spreading
        factor = 0.5). For our purposes, this will fold the attenuation factor
        into the same single factor that accounts for geometric spreading.

        .. note::
            This does not take into account the variation in amplitude.

        TODO Plot geometric spreading with best-fitting curve

        :type max_amplitudes: list
        :param max_amplitudes: list of maximum amplitudes for each trace in the
            stream. IF not given, will be calculated on the fly. Optional incase
            this has been calculated already and you want to save time
        :rtype: list
        :return: scale factor per trace in the stream based on theoretical
            geometrical spreading factor. This is meant to be MULTIPLIED by the
            data arrays
        """
        raise NotImplementedError("This function is currently work in progress")

        logger.info("calculating geometrical spreading for amplitude "
                    "normalization")

        # Create a sinusoidal function based on distances in degrees
        sin_del = np.sin(np.array(self.dist) / (180 / np.pi))

        # !!! TODO Stein and Wysession and figure out vector names
        if self.geometric_spreading_k_val is not None:
            k_vector = self.max_amplitudes * \
                       (sin_del ** self.geometric_spreading_factor)
            k_val = np.median(k_vector)
        else:
            k_val = self.geometric_spreading_k_val

        w_vector = k_val / (sin_del ** self.geometric_spreading_factor)

        return w_vector

    def get_sorted_idx(self):
        """
        Sort the source-receiver pairs by the chosen parameters. Except for in
        `default` sorting, always maintains component ordering defined by the
        user, i.e., for the same station, components will be ordered by the
        `component` parameter.

        :rtype: np.array
        :return: returns an indexing list which is sorted based on the user
            defined `sort_by` argument.
        """
        logger.info(f"determining sort order with parameter: {self.sort_by}")

        # Retain input ordering but run sorting to allow reversing
        if "default" in self.sort_by:
            sorted_idx = sorted(self.idx, reverse=self.stats.reverse_sort)
        # Sort alphabetically BUT allow component sorting
        elif "alphabetical" in self.sort_by:
            sorted_idx = self._sort_by_alphabetical()
        # Sort by dist, az or baz
        else:
            components = self._get_component_order()
            # Azimuthal sorts are allowed to start at a value other than 0
            # Backazimuth needs to come first because 'azimuth' can return both
            if "backazimuth" in self.sort_by:
                sort_list = (self.backazimuths - self.azimuth_start_deg) % 360
            elif "azimuth" in self.sort_by:
                sort_list = (self.azimuths - self.azimuth_start_deg) % 360
            elif "distance" in self.sort_by:
                sort_list = self.distances

            # Sort by the values, AND the components (e.g., Z->R->T or whatever)
            _, _, sorted_idx = \
                zip(*sorted(zip(sort_list, components, self.idx),
                            reverse=self.stats.reverse_sort)
                    )
        sorted_idx = np.array(list(sorted_idx))

        return sorted_idx

    def _sort_by_alphabetical(self):
        """
        Sort by full station name in order. That is, network is sorted, then
        within network, station is sorted, and so on. Components are sorted
        last and NOT in alphabetical order. They are ordered based on the
        user-input `components` parameter.

        :rtype: np.array
        :return: indexing of networks and stations sorted alphabetically
        """
        networks = [_.split(".")[0] for _ in self.station_ids]
        stations = [_.split(".")[1] for _ in self.station_ids]
        locations = [_.split(".")[2] for _ in self.station_ids]
        components = self._get_component_order()
        zipped = zip(networks, stations, locations, components, self.idx)
        arr_out = np.array(
            list(zip(*sorted(zipped, reverse=self.stats.reverse_sort)))
        )

        return arr_out[-1].astype(int)

    def _get_component_order(self):
        """
        When we are sorting, we want components to be sorted by the
        preferred order (i.e, ZRT, Z index is 0, T is 2), which means the
        components have to be reordered to adhere to the sorting algorithm.
        ALSO we need to ensure that we honor component order even if
        everything else is being sorted reverse alphabetically so the
        components entry needs to be the opposite of stats.reverse_sort

        :rtype: list
        :return: components converted to integer values which can be used for
            sorting algorithms.
        """
        channels = [_.split(".")[3] for _ in self.station_ids]
        # ASSUMING component is the last entry in the channel, SEED convention
        components = [_[-1] for _ in channels]
        if self.stats.reverse_sort:
            comp_list = self.components
        else:
            comp_list = self.components[::-1]

        # Can't list comp here incase `comp_list` does NOT contain one of the
        # available components, which will throw a ValueError. Use a random
        # number to catch these.
        numbered_components = []
        for comp in components:
            try:
                numbered_components.append(comp_list.index(comp))
            except ValueError:
                numbered_components.append(-999)

        return numbered_components

    def get_y_axis_positions(self):
        """
        Determine how seismograms are vertically separated on the Y-Axis when
        plotting. Allows for constant separation, or absolute separation based
        on distance or (back)azimuth separation.

        :rtype: np.array
        :return: an array of actual y-axis positions that can be passed directly
            to plt.plot() (or similar)
        """
        logger.info(f"determining y-axis positioning for sort: {self.sort_by}")

        # Default weights provides constant `y_axis_spacing` between seismos
        if self.sort_by == "default" or "abs_" not in self.sort_by:
            y_range = np.arange(0, len(self.sorted_idx), 1)
            y_axis = self.y_axis_spacing * y_range
        # Absolute y-axis (i.e., absolute distance or (back)azimuth) will just
        # plot the actual distance or azimuth value
        else:
            if "backazimuth" in self.sort_by:
                y_axis = self.backazimuths
            elif "azimuth" in self.sort_by:
                y_axis = self.azimuths
            elif "distance" in self.sort_by:
                y_axis = self.distances

        return y_axis

    def process_st(self):
        """
        Preprocess the Stream with optional filtering in place.

        .. note::
            Data in memory will be irretrievably altered by running preprocess.

        TODO Add feature to allow list-like periods to individually filter
            seismograms. At the moment we just apply a blanket filter.
        """
        if self.preprocess is None:
            logger.info("no preprocessing applied")
            return
        elif self.preprocess == "st":
            logger.info(f"preprocessing {len(self.st)} `st` waveforms")
            preprocess_list = [self.st]
        elif self.preprocess == "st_syn":
            logger.info(f"preprocessing {len(self.st_syn)} `st_syn` waveforms")
            preprocess_list = [self.st_syn]
        elif self.preprocess == "both":
            logger.info(f"preprocessing {len(self.st) + len(self.st_syn)} "
                  f"`st` and `st_syn` waveforms")
            preprocess_list = [self.st, self.st_syn]

        for st in preprocess_list:
            # Fill any data gaps with mean of the data, do it on a trace by 
            # trace basis to get individual mean values
            for tr in st:
                tr.trim(starttime=tr.stats.starttime, endtime=tr.stats.endtime,
                        pad=True, fill_value=tr.data.mean())
            st.detrend("demean")
            st.taper(max_percentage=0.05, type="cosine")

            # Zero pad start and end of data if requested by user
            if self.zero_pad_s:
                _start, _end = self.zero_pad_s
                logger.info(f"padding zeros to traces with {_start}s before "
                            f"and {_end}s after")
                for idx, tr in enumerate(st):
                    tr.trim(starttime=tr.stats.starttime - _start,
                            endtime=tr.stats.endtime + _end,
                            pad=True, fill_value=0)

            # Allow multiple filter options based on user input
            # Min period but no max period == low-pass
            if self.max_period_s is not None and self.min_period_s is None:
                logger.info(f"apply lowpass filter w/ cutoff "
                            f"{1/self.max_period_s}")
                st.filter("lowpass", freq=1/self.max_period_s, zerophase=True)
            # Max period but no min period == high-pass
            elif self.min_period_s is not None and self.max_period_s is None:
                logger.info(f"apply highpass filter w/ cutoff "
                            f"{1/self.min_period_s}")
                st.filter("highpass", freq=1/self.min_period_s, zerophase=True)
            # Both min and max period == band-pass
            elif self.min_period_s is not None and \
                    self.max_period_s is not None:
                logger.info(f"applying bandpass filter w/ "
                      f"[{1/self.max_period_s}, {self.min_period_s}]")
                st.filter("bandpass", freqmin=1/self.max_period_s,
                            freqmax=1/self.min_period_s, zerophase=True)
            else:
                logger.info("no filtering applied")

            # Integrate and differentiate N number of times specified by user
            st.detrend("simple")
            if self.integrate != 0:
                if self.integrate < 0:
                    func = "differentiate"
                elif self.integrate > 0:
                    func = "integrate"
            for i in range(np.abs(self.integrate)):
                logger.info(f"{func} all waveform data x{abs(self.integrate)}")
                getattr(st, func)()

    def plot(self, subset=None, page_num=None, **kwargs):
        """
        Plot record sections based on all the available information

        :type subset: list of int
        :param subset: subset of `sorted_idx` if there are too many waveforms to
            plot on one page (set by `max_traces_per_rs`). e.g., to get the
            first 10 entries, subset=[0,10]
        :type page_num: int
        :param page_num: If multiple pages are required due to exceeding 
            `max_traces_per_rs`, page_num will make adjustments to the figure
            name and title to differentiate different pages of the same record
            section
        """
        if subset is None:
            start, stop = 0, None  # None will allow full list traversal
            nwav = len(self.sorted_idx)
        else:
            start, stop = subset
            nwav = stop - start

        logger.info(f"plotting record section for {nwav} waveforms")

        # Do a text output of station information so the user can check
        # that the plot is doing things correctly
        logger.debug("plotting line check starting from bottom (y=0)")
        logger.debug("\nIDX\tY\t\tID\tDIST\tAZ\tBAZ\tTSHIFT\tYABSMAX")
        self.f, self.ax = plt.subplots(figsize=self.figsize)

        log_str = "\n"
        # Allow choosing observed or synthetic data, defaults to observed
        # Allow different waveform looks based on observed vs. synthetic
        for choice in ["st", "st_syn"]:
            if getattr(self, choice) is None:
                continue
            # Main plotting call. Indexes should already be sorted
            for y_idx, idx in enumerate(self.sorted_idx[start:stop]):
                # Absolute scaling requires plotting actual dist or az values
                # Relative scaling just requires a linear index to stack seismos
                if "abs_" in self.sort_by:
                    y_index = idx
                else:
                    y_index = y_idx + start
                log_str += self._plot_trace(idx=idx, y_index=y_index,
                                            choice=choice, **kwargs)

        logger.debug(log_str)
        # Change the aesthetic look of the figure, should be run before other
        # set functions as they may overwrite what is done here
        self._set_plot_aesthetic()

        # Partition the figure by user-specified azimuth bins 
        if self.sort_by and "azimuth" in self.sort_by:
            self._plot_azimuth_bins()

        # Finalize the figure accoutrements
        self._plot_title(nwav=nwav, page_num=page_num)
        self._plot_axes(start=start, stop=stop)

        if self.save:
            # Allow appending page numbers to figure names,
            # e.g., default.png -> default_01.png
            if page_num:
                fid, ext = os.path.splitext(self.save)
                save_fid = f"{fid}_{page_num:0>2}{ext}"
            else:
                save_fid = self.save
            logger.info(f"saving figure to {save_fid}")
            plt.savefig(save_fid)
        if self.show:
            plt.show()

    def _plot_trace(self, idx, y_index, choice="st"):
        """
        Plot a single trace on the record section, with amplitude scaling,
        time shifts, etc. Observed waveforms are black, synthetics are red.

        :type idx: int
        :param idx: index of the trace to plot and all trace-ordered values like
            amplitude scaling and time shifts
        :type y_index: int
        :param y_index: numerical order which this trace was plotted, as each
            trace is plotted on top of the previous one. y_index should be
            iterated linearly in the loop that calls this function.
        :type choice: str
        :param choice: choice of 'st' or 'st_syn' depending on whether you want
            to plot the observed or synthetic waveforms
        """
        linewidth = self.kwargs.get("linewidth", .25)

        # Used to differentiate the two types of streams for plotting diffs
        choices = ["st", "st_syn"]
        assert (choice in choices)
        c = choices.index(choice)
        tr = getattr(self, choice)[idx]  # i.e., tr = self.st[idx]

        # Plot actual data on with amplitude scaling, time shift, and y-offset
        tshift = self.time_shift_s[idx]
        
        # These are still the entire waveform. Make sure we honor zero padding
        # and any time shift applied
        x = tr.times() + tshift
        if self.zero_pad_s is not None:
            x -= self.zero_pad_s[0]  # index 0 is start, index 1 is end
        y = tr.data / self.amplitude_scaling[idx] + self.y_axis[y_index]

        # Truncate waveforms to get figure scaling correct. 
        start, stop = self.xlim[idx]
        x = x[start:stop]
        y = y[start:stop]

        self.ax.plot(x, y, c=["k", "r"][c], linewidth=linewidth, zorder=10)

        # Sanity check print station information to check against plot
        log_str = (f"{idx}"
                   f"\t{int(self.y_axis[y_index])}"
                   f"\t{tr.get_id():<6}"
                   f"\t{self.distances[idx]:6.2f}"
                   f"\t{self.azimuths[idx]:6.2f}"
                   f"\t{self.backazimuths[idx]:6.2f}"
                   f"\t{self.time_shift_s[idx]:4.2f}"
                   f"\t{self.max_amplitudes[idx]:.2E}\n"
                   )

        # Retain some stats for global plot args
        self.stats.xmin.append(x.min())
        self.stats.xmax.append(x.max())
        self.stats.ymin.append(y.min())
        self.stats.ymax.append(y.max())

        return log_str

    def _plot_azimuth_bins(self):
        """
        If plotting by azimuth, create visual bin separators so the user has
        a visual understanding of radiation patterns etc.
        """
        azimuth_binsize = self.kwargs.get("azimuth_binsize", 45)

        # Look of the azimuth bin lines
        c = self.kwargs.get("azimuth_bin_c", "r")
        lw = self.kwargs.get("azimuth_bin_lw", .75)
        ls = self.kwargs.get("azimuth_bin_ls", "-")
        z = self.kwargs.get("azimuth_bin_zorder", 5)

        azimuth_bins = np.arange(self.azimuth_start_deg,
                                 self.azimuth_start_deg + 360,
                                 azimuth_binsize)
        # Make sure that the bins go from 0 <= azimuth_bins <= 360
        azimuth_bins = azimuth_bins % 360

        # In an absolute plot, the y values are simply the azimuth bins
        if "abs" in self.sort_by:
            y_vals = azimuth_bins
            s_vals = [f"{azbin}{DEG}" for azbin in azimuth_bins]
        # In a relative plot, the y values need to be found between seismograms
        else:
            s_vals, y_vals = [], []
            # Brute force determine where these azimuth bins would fit into the 
            # actual plotted azimuths, draw a line between adjacent seismograms
            # i.e., iterating and looking if the bin value fits between
            # two adjacent azimuth values
            for azbin in azimuth_bins:
                # Edge case for 0 deg azimuth bin since azimuth wraps on 0
                # Look for the top minimum azimuth value, place line above that
                if azbin == 0:
                    dy = abs(self.y_axis[1] - self.y_axis[0])
                    azis = self.azimuths[self.sorted_idx]
                    i = max(np.where(azis == azis.min())[0])
                    y_vals.append(self.y_axis[i] + dy/2)
                else:
                    for i, idx in enumerate(self.sorted_idx[1:]):
                        j = i + 1
                        idx_minus_one = self.sorted_idx[i]
                        azi_low = self.azimuths[idx]
                        azi_high = self.azimuths[idx_minus_one]
                        # Break if bin is in between azi values for two seismos
                        if azi_low <= azbin <= azi_high:
                            break
                    else:
                        continue
                    # Mean gives the space in between two adjacent seismograms
                    y_vals.append(np.mean([self.y_axis[i],  self.y_axis[j]]))

                s_vals.append(f"{azbin}{DEG}")

        for y, (s_val, y_val) in enumerate(zip(s_vals, y_vals)):
            # Dealing with the case where two bins occupy the same space,
            # only plot the first one that occurs
            if y_val == y_vals[y-1]:
                continue
            plt.axhline(y=y_val, c=c, linewidth=lw, linestyle=ls, zorder=z)
            plt.text(x=max(self.stats.xmax), y=y_val, s=s_val, c=c, ha="right")

    def _plot_axes(self, start=0, stop=None):
        """
        Contains the logic in how to handle the x- and y-axis labels, ticks etc.

        Logic options for how to plot the y-axis:
        - Relative: Relative plotting means the yticklabels will get replaced
            by station information. This can either be on the right or left
            side but will be attached to the actual axis
        - Absolute: Absolute plotting means the y-axis actual represents
            something (e.g., distance, azimuth). That means labels can not
            replace the y-ticks and they have to be placed within the figure

        .. note::
            Relative plotting can also place labels OFF the y-axis, at which
            point the y-axis is turned off because it contains no usable
            information

        :type start: int
        :param start: optional starting index for creating text labels
        :type stop: int
        :param stop: optional stop index for creating text labels
        """
        # Sort out how the y-axis will be labeled with station information and
        # dist/azimuth if we're doing absolute sorting
        if "abs_" in self.sort_by:
            self._set_y_axis_absolute()
            if self.y_label_loc is not None:
                # By default, place absolute sorted labels at xmin
                if self.y_label_loc == "default":
                    loc = "x_min"
                else:
                    loc = self.y_label_loc
                self._set_y_axis_text_labels(start=start, stop=stop, loc=loc)
        elif self.y_label_loc is not None:
            # By default, relative plotting shows labels on the right side
            if self.y_label_loc == "default":
                loc = "y_axis"
            else:
                loc = self.y_label_loc
            self._set_y_axis_text_labels(start=start, stop=stop, loc=loc)
        logger.info(f"placing station labels on y-axis at: {loc}")

        # User requests that we turn off y-axis
        if self.y_label_loc is None:
            logger.info(f"user requests turning off y-axis w/ "
                        f"`y_label_loc`== None")
            self.ax.get_yaxis().set_visible(False)
        # OR EDGE CASE: Relative plotting but labels are not placed on the
        # y-axis, turn off the y-axis cause it doesn't mean anything anymore
        elif self.y_label_loc not in ["default", "y_axis", "y_axis_right"] and \
                "abs" not in self.sort_by:
            logger.info("turning off y-axis as it contains no information")
            self.ax.get_yaxis().set_visible(False)

        # Reverse the y-axis if we are doing absolute y-axis and reversing
        if "abs_" in self.sort_by and "_r" in self.sort_by:
            logger.info("user requests inverting y-axis with absolute "
                        "reverse sort")
            self.ax.invert_yaxis()

        # X-axis label is different if we time shift
        if self.time_shift_s.sum() == 0:
            plt.xlabel("Time [s]")
        else:
            plt.xlabel("Relative Time [s]")

        # Allow user defined x-axis limits
        if self.xlim_s is None:
            self.ax.set_xlim([min(self.stats.xmin), max(self.stats.xmax)])
        else:
            self.ax.set_xlim(self.xlim_s)

        plt.tight_layout()

    def _set_y_axis_absolute(self):
        """
        If 'abs_' in sort_by, then the Y-axis should be shown in absolute scale.
        That means we need to format the text labels, add some labelling etc.
        """
        # Reset tick label size to be larger to match absolute x-axis size
        ytick_fontsize = self.kwargs.get("ytick_fontsize", 12)
        self.ax.tick_params(axis="y", labelsize=ytick_fontsize)

        if "distance" in self.sort_by:
            if "km" in self.distance_units:
                ytick_minor = self.kwargs.get("ytick_minor", 25)
                ytick_major = self.kwargs.get("ytick_major", 100)
            elif "deg" in self.distance_units:
                ytick_minor = self.kwargs.get("ytick_minor", 0.25)
                ytick_major = self.kwargs.get("ytick_major", 1.0)
            ylabel = f"Distance [{self.distance_units}]"
        elif "azimuth" in self.sort_by:
            ytick_minor = self.kwargs.get("ytick_minor", 45)
            ytick_major = self.kwargs.get("ytick_major", 90)
            ylabel = f"Azimuth [{DEG}]"

        # Set ytick label major and minor which is either dist or az
        self.ax.yaxis.set_major_locator(MultipleLocator(ytick_major))
        self.ax.yaxis.set_minor_locator(MultipleLocator(ytick_minor))
        self.ax.set_ylabel(ylabel)

    def _set_y_axis_text_labels(self, start=0, stop=-1, loc="y_axis"):
        """
        Plot a text label next to each trace describing the station,
        azimuth and source-receiver distance. We need to do this all at once
        because we have to replace all ticks and tick labels at once.

        .. note::
            if using the 'subset' option in plot, need to also tell the y-axis
            plotter that we are only plotting a subset of data by using the
            `start` and `stop` parameters

        :type start: int
        :param start: starting index for plotting, default to start 0
        :type stop: int
        :param stop: stop index for plotting, default to end -1
        :type loc: str
        :param loc: location to place the y_axis text labels, available:
            - y_axis: Place labels along the y-axis (left side of the figure)
                Will replace the actual y-tick labels so not applicable for
                absolute sorting which requries the y-axis labels
            - y_axis_right: same as `y_axis` but set on the right side of figure
            - x_min: Place labels on the waveforms at the minimum x value
            - x_max: Place labels on the waveforms at the maximum x value
        """
        c = self.kwargs.get("y_label_c", "k")
        fontsize = self.kwargs.get("y_label_fontsize", 10)

        y_tick_labels = []
        for idx in self.sorted_idx[start:stop]:
            str_id = self.station_ids[idx]
            if self.sort_by is not None and "backazimuth" in self.sort_by:
                # This is named `str_az` but it's actually backazimuths
                str_az = f"{self.backazimuths[idx]:6.2f}{DEG}"
            else:
                str_az = f"{self.azimuths[idx]:6.2f}{DEG}"

            # Allow degree distance to use the unicode * symbol
            if self.distance_units == "deg":
                du = DEG
            else: 
                du = self.distance_units

            str_dist = f"{self.distances[idx]:5.2f}{du}"

            # Looks something like: NN.SSS.LL.CC|30*|250.03km
            label = \
                f"{str_id:>{self.stats.longest_id}}|{str_az:>8}|{str_dist:>8}"
            # Add time shift if we have shifted at all
            if self.time_shift_s[idx] != 0:
                label += f"|{self.time_shift_s[idx]:.2f}s"
            y_tick_labels.append(label)

        if "y_axis" in loc:
            # For relative plotting (not abs_), replace y_tick labels with
            # station information
            self.ax.set_yticks(self.y_axis[start:stop])
            self.ax.set_yticklabels(y_tick_labels)
        else:
            # Trying to figure out where the min or max X value is on the plot
            if loc == "x_min":
                ha = "left"
                func = min
                x_val = func(self.stats.xmin)
            elif loc == "x_max":
                ha = "right"
                func = max
                x_val = func(self.stats.xmax)
            if self.xlim_s is not None:
                x_val = func([func(self.xlim_s), x_val])

            # Plotting y-axis labels for absolute scales
            if len(self.y_axis) == len(self.st):
                for idx, s in zip(self.sorted_idx[start:stop], y_tick_labels):
                    plt.text(x=x_val, y=self.y_axis[idx], s=s, ha=ha, c=c,
                             fontsize=fontsize)
            # Plotting y-axis labels for relative scales
            elif len(self.y_axis) == len(y_tick_labels):
                for y, s in zip(self.y_axis, y_tick_labels):
                    plt.text(x=x_val, y=y, s=s, ha=ha, c=c, fontsize=fontsize)

        if loc == "y_axis_right":
            self.ax.yaxis.tick_right()
            self.ax.yaxis.set_label_position("right")

    def _plot_title(self, nwav=None, page_num=None):
        """
        Create the title of the plot based on event and station information
        Allow dynamic creation of title based on user input parameters

        TODO Can we make this two-column to save space?

        :type nwav: int
        :param nwav: if using subset, the title needs to know how many waveforms
            it's showing on the page. self.plot() should tell it
        :type page_num: int
        :param page_num: same as nwav, we need to know what page number were on
        """
        # Defines the number of waveforms plotted on a single page, allowing
        # for subsets per page
        if nwav is None:
            nwav = len(self.sorted_idx)

        # Allow appending page numbers to title
        title_top = "RECORD SECTION"
        if page_num:
            title_top += f" PAGE {page_num:0>2}"

        # The y-label will show baz or az depending on user choice, distinguish
        if self.sort_by is not None and "backazimuth" in self.sort_by:
            az_str = "BAZ"
        else:
            az_str = "AZ"

        # Get the unique components that have been plotted, only
        cmp = "".join(np.unique([self.st[i].stats.component
                                 for i in self.sorted_idx]))

        # Y_FMT will include time shift IF there are time shifts
        y_fmt = f"Y_FMT: NET.STA.LOC.CHA|{az_str}|DIST"
        if self.time_shift_s.sum() != 0:
            y_fmt += "|TSHIFT"

        title = "\n".join([
            title_top,
            f"{'/' * len(title_top*2)}",
            f"ORIGINTIME: {min([tr.stats.starttime for tr in self.st])}",
            f"{y_fmt}",
            f"NWAV: {nwav}; NEVT: {self.stats.nevents}; "
            f"NSTA: {self.stats.nstation}; COMP: {cmp}",
            f"SORT_BY: {self.sort_by}; "
            f"SCALE_BY: {self.scale_by} * {self.amplitude_scale_factor}",
            f"FILT: [{self.min_period_s}, {self.max_period_s}]s; "
            f"MOVE_OUT: {self.move_out or 0}{self.distance_units}/s",
        ])
        self.ax.set_title(title)

    def _set_plot_aesthetic(self):
        """
        Give a nice look to the output figure by creating thick borders on the
        axis, adjusting fontsize etc. All plot aesthetics should be placed here
        so it's easiest to find.

        .. note::
            This was copy-pasted from Pyatoa.visuals.insp_plot.default_axes()
        """
        ytick_fontsize = self.kwargs.get("ytick_fontsize", 8)
        xtick_fontsize = self.kwargs.get("xtick_fontsize", 12)
        tick_linewidth = self.kwargs.get("tick_linewidth", 1.5)
        tick_length = self.kwargs.get("tick_length", 5)
        tick_direction = self.kwargs.get("tick_direction", "in")
        label_fontsize = self.kwargs.get("label_fontsize", 10)
        axis_linewidth = self.kwargs.get("axis_linewidth", 2.)
        title_fontsize = self.kwargs.get("title_fontsize", 10)
        xtick_minor = self.kwargs.get("xtick_minor", 25)
        xtick_major = self.kwargs.get("xtick_major", 100)
        spine_zorder = self.kwargs.get("spine_zorder", 8)

        # Re-set font sizes for labels already created
        self.ax.title.set_fontsize(title_fontsize)
        self.ax.tick_params(axis="both", which="both", width=tick_linewidth,
                            direction=tick_direction, length=tick_length)
        self.ax.tick_params(axis="x", labelsize=xtick_fontsize)
        self.ax.tick_params(axis="y", labelsize=ytick_fontsize)
        self.ax.xaxis.label.set_size(label_fontsize)

        # Thicken up the bounding axis lines
        for axis in ["top", "bottom", "left", "right"]:
            self.ax.spines[axis].set_linewidth(axis_linewidth)

        # Set spines above azimuth bins
        for spine in self.ax.spines.values():
            spine.set_zorder(spine_zorder)

        # Set xtick label major and minor which is assumed to be a time series
        self.ax.xaxis.set_major_locator(MultipleLocator(xtick_major))
        self.ax.xaxis.set_minor_locator(MultipleLocator(xtick_minor))

        plt.grid(visible=True, which="major", axis="x", alpha=0.5, linewidth=1)
        plt.grid(visible=True, which="minor", axis="x", alpha=0.2, linewidth=.5)


def parse_args():
    """
    Parse command line arguments to set record section parameters dynamically
    This arg parser provides a simplified interface for working with plotw_rs
    BUT it limits the flexibility of the code because things like long lists 
    are prohibitively verbose and not included in the arguments.

    Kwargs can be passed in in the same format thanks to:
        https://stackoverflow.com/questions/37367331/is-it-possible-to-use-\
                argparse-to-capture-an-arbitrary-set-of-optional-arguments

    .. note::
        Not all parameters are set here, some are left as default values
        Also some parameters are set different than the class defaults, so that
        when the user runs record_section.py without arguments, they get a
        reasonable result

    .. note::
        Do NOT use the command line if you want to exploit the expanded
        capabilities of the record section plotter, rather script it or call
        from an interactive environment.
    """
    parser = argparse.ArgumentParser(
        description="Input basic record section params",
        formatter_class=argparse.RawTextHelpFormatter,
                                     )

    parser.add_argument("-p", "--pysep_path", default="./", type=str, nargs="?",
                        help="path to Pysep output, which is expected to "
                             "contain trace-wise SAC waveform files which will "
                             "be read")
    parser.add_argument("--syn_path", default=None, type=str, nargs="?",
                        help="path to SPECFEM generated synthetics. Also "
                             "requires --cmtsolution_file and --stations_file")
    parser.add_argument("--cmtsolution", default=None, type=str, nargs="?",
                        help="required for synthetics, path to the CMTSOLUTION "
                             "file used to generate SPECFEM synthetics")
    parser.add_argument("--stations", default=None, type=str, nargs="?",
                        help="required for synthetics, path to the STATIONS "
                             "file used to generate SPECFEM synthetics")
    parser.add_argument("--sort_by", default="distance", type=str, nargs="?",
                        help=textwrap.dedent("""
            How to sort the Y-axis of the record section
            - None: Not set, don't sort, just iterate directly through Stream
            - 'alphabetical': sort alphabetically
            - 'azimuth': sort by source-receiver azimuth (deg) with constant
                vertical spacing
            - 'backazimuth': sort by source-receiver backazimuth (deg) with
                constant vertical spacing. Requires `azimuth_start_deg`
            - 'distance': sort by source-receiver distance (km) with constant
                vertical spacing. Requires `azimuth_start_deg` AND
                `distance_units`
            - 'abs_distance': absolute vertical spacing of waveforms defined by
                source-receiver distance (km). Requires `distance_units`
            - 'abs_azimuth': absolute vertical spacing of waveforms defined
                by source-receiver azimuth (deg). Requires
                `azimuth_start_deg`
            - 'abs_backazimuth': absolute vertical spacing of waveforms by
                source-receiver backazimuth (deg).
            - '*_r': Add a '_r' to any of the values about to REVERSE the sort,
                e.g., alphabetical_r sort will go Z->A"
                """)
                        )
    parser.add_argument("--scale_by", default="normalize", type=str, nargs="?",
                        help=textwrap.dedent("""
            How to sort the Y-axis of the record section
            - None: Not set, no amplitude scaling, waveforms shown raw
            - 'normalize': scale each trace by the maximum amplitude,
                i.e., > a /= max(abs(a))  # where 'a' is time series amplitudes
            - 'global_norm': scale by the largest amplitude to be displayed on
                the screen. Will not consider waveforms which have been 
                excluded on other basis (e.g., wrong component)
            - 'geometric_spreading': scale amplitudes globally by predicting the
                expected geometric spreading amplitude reduction and correcting
                for this factor. Requires `geometric_spreading_factor`, optional
                `geometric_spreading_k_val`
                """)
                        )
    parser.add_argument("--time_shift_s", default=None, type=float, nargs="?",
                        help="Set a constant time shift in unit: seconds")
    parser.add_argument("--move_out", default=None, type=float, nargs="?",
                        help="Set a constant velocity-based move out in units:"
                             "`distance_units`/s")
    parser.add_argument("--min_period_s", default=None, type=float, nargs="?",
                        help="Minimum filter period in unit seconds.")
    parser.add_argument("--max_period_s", default=None, type=float, nargs="?",
                        help="Maximum filter period in unit seconds.")
    parser.add_argument("--max_traces_per_rs", default=None, nargs="?",
                        help="Max waveforms to show on one page. If the number "
                             "of waveforms to show exceeds this value, "
                             "multiple pages will be saved")
    parser.add_argument("--integrate", default=0, type=int, nargs="?",
                        help="Integrate (positive values) or differentiate "
                             "(negative values) `integrate` times. e.g., -2 "
                             "will differentiate all waveforms twice.")
    parser.add_argument("--xlim_s", default=None, type=int, nargs="+",
                        help="Min and max x limits in units seconds. Input as "
                             "a list, e.g., --xlim_s 10 200 ...")
    parser.add_argument("--zero_pad_s", default=None, type=float, nargs="+",
                        help="Zero pad the start and end of each trace. Input "
                             "as two values in units of seconds."
                             "i.e., --zero_pad_s `START` `END` or"
                             "e.g., --zero_pad_s 30 0 to pad 30s before start "
                             "0s at end of trace")
    parser.add_argument("--components", default="ZRTNE12", type=str, nargs="?",
                        help="Ordered component list specifying 1) which "
                             "components will be shown, and in what order. "
                             "e.g., 'ZR' will show only Z and R components "
                             "with Z sorted above R.")
    parser.add_argument("--y_label_loc", default="default", type=str, nargs="?",
                        help="Station information label location on the y-axis")
    parser.add_argument("--y_axis_spacing", default=1, type=float, nargs="?",
                        help="For relative sorting, the y-axis spacing between "
                             "adjacent seismograms. If waveforms are "
                             "normalized then a default value of 1 is usually "
                             "normalized then a default value of 1 is usually "
                             "fine")
    parser.add_argument("--amplitude_scale_factor", default=1, type=float,
                        nargs="?",
                        help="A user dial allowing waveform amplitudes to be"
                             "scaled by an arbitrary float value")
    parser.add_argument("--azimuth_start_deg", default=0, type=float,
                        nargs="?",
                        help="When sorting by azimuth, choose the default "
                             "starting value, with azimuth 0 <= az <= 360")
    parser.add_argument("--distance_units", default="km", type=str,
                        nargs="?", help="Set units when sorting by distance")
    parser.add_argument("--save", default="./record_section.png", type=str,
                        nargs="?",
                        help="Path to save the resulting record section fig")
    parser.add_argument("-o", "--overwrite", default=False, action="store_true",
                        help="overwrite existing figure if path exists")
    parser.add_argument("--log_level", default="DEBUG", type=str,
                        help="verbosity of logger: 'WARNING', 'INFO', 'DEBUG'")

    # Keyword arguments can be passed directly to the argparser in the same 
    # format as the above kwargs (e.g., --linewidth 2), but they will not have 
    # help messages or type checking
    parsed, unknown = parser.parse_known_args()
    for arg in unknown:
        if arg.startswith(("-", "--")):
            parser.add_argument(arg.split("=")[0])

    return parser.parse_args()


def plotw_rs(*args, **kwargs):
    """
    Main call function, replacing `plotw_rs`. Run the record section plotting
    functions in order. Contains the logic for breaking up figure into multiple
    pages.

    .. note::
        All arguments should be parsed into the argparser, *args and **kwargs
        are just there to keep the IDE happy
    """
    _start = datetime.now()
    logger.info(f"starting record section plotter")

    rs = RecordSection(*args, **kwargs)
    rs.process_st()
    rs.get_parameters()
    # Simple case where all waveforms will fit on one page
    if len(rs.sorted_idx) <= rs.max_traces_per_rs:
        rs.plot()
    # More complicated case where we need to split onto multiple pages
    else:
        for i, start in enumerate(np.arange(0, len(rs.st),
                                            rs.max_traces_per_rs)):
            stop = start + rs.max_traces_per_rs
            # Case where the num waveforms is less than max_traces_per_rs
            if stop < rs.max_traces_per_rs:
                stop = len(rs.st)
            rs.plot(subset=[start, stop], page_num=i+1)

    _end = datetime.now()
    logger.info(f"finished record section in t={(_end - _start)}s")


def main():
    """
    Convenience 'main' function to play nice with entry scripts

    .. rubric::
        $ recsec -h
    """
    plotw_rs(**vars(parse_args()))


if __name__ == "__main__":
    main()

