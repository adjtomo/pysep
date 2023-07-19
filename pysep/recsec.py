#!/usr/bin/evn python3
"""
RECord SECtion plotting tool for seismic waveforms (observed and synthetic)

This is a refactor of Pysep's Python utility `plotw_rs`, a record section
plotting script. The intent of this script is to plot multiple time series'
based on source-receiver characteristics (i.e., src-rcv distance, backazimuth).

.. note:: 
    Code History:  
    - Written by Carl Tape (11/21/2011) and Yun Wang (11/2011) in Matlab  
    - Translated to Python by Nealy Sims (1/2021)  
    - Upgraded by Aakash Gupta (9/2021)  
    - Refactored by Bryant Chow (3/2022)  
    - Currently maintained by adjTomo Dev Team

.. requires::
    obspy >= 1.2 (expected to bring in numpy and matplotlib)

.. rubric:: 

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
        >>> from recsec import plotw_rs  # NOQA
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
from pysep.utils.cap_sac import origin_time_from_sac_header, SACDICT
from pysep.utils.io import read_sem, read_sem_cartesian
from pysep.utils.curtail import subset_streams
from pysep.utils.plot import plot_geometric_spreading, set_plot_aesthetic

# Unicode degree symbol for plot text
DEG = u"\N{DEGREE SIGN}"


class Dict(dict):
    """Simple dictionary overload for nicer get/set attribute characteristics"""
    def __setattr__(self, key, value):
        self[key] = value

    def __getattr__(self, key):
        return self[key]


class RecordSection:
    """
    Record section plotting tool which takes ObsPy streams and:

    1) preprocesses and filters waveforms,
    2) sorts source-receiver pairs based on User input,
    3) produces record section waveform figures.
    """
    def __init__(self, pysep_path=None, syn_path=None, stations=None,
                 cmtsolution=None, st=None, st_syn=None, sort_by="default",
                 scale_by=None, time_shift_s=None, zero_pad_s=None,
                 move_out=None, min_period_s=None, max_period_s=None,
                 preprocess=None, max_traces_per_rs=None, integrate=0,
                 xlim_s=None, components="ZRTNE12", y_label_loc="default",
                 y_axis_spacing=1, amplitude_scale_factor=1,
                 azimuth_start_deg=0., distance_units="km", tmarks=None,
                 geometric_spreading_factor=0.5, geometric_spreading_k_val=None,
                 geometric_spreading_exclude=None,
                 geometric_spreading_ymax=None, geometric_spreading_save=None,
                 figsize=(9, 11), show=True, save="./record_section.png",
                 overwrite=False, log_level="DEBUG", cartesian=False,
                 synsyn=False, **kwargs):
        """
         .. note::
            Used for reading in Pysep-generated waveforms

        :type pysep_path: str
        :param pysep_path: path to Pysep output, which is expected to contain
            trace-wise SAC waveform files which will be read in.

        .. note::
            Used for reading in SPECFEM-generated synthetic waveforms

        :type syn_path: str
        :param syn_path: full path to directory containing synthetic
            seismograms that have been outputted by SPECFEM. The synthetics
            are expected in the format: '*??.*.*.sem?*', which generally matches
            SPECFEM files like 'NZ.BFZ.BXE.semd'
        :type stations: str
        :param stations: full path to STATIONS file used to define the station
            coordinates. Format is dictated by SPECFEM
        :type cmtsolution: str
        :param cmtsolution: required for synthetics, full path to SPECFEM source
            file, which was used to generate SPECFEM synthetics. Example
            filenames are CMTSOLUTION, FORCESOLUTION, SOURCE.
        :type cartesian: bool
        :param cartesian: lets RecSec know that the domain is set in Cartesian
            coordinates (SPECFEM2D or SPECFEM3D_CARTESIAN in UTM), such that
            data/metadata reading will need to adjust acoordingly because the
            ObsPy tools used to read metadata will fail for domains defined in
            XYZ
        :type synsyn: bool
        :param synsyn: flag to let RecSec know that we are plotting two sets
            of synthetic seismograms. Such that both `pysep_path` and `syn_path`
            will be both attempt to read in synthetic data. Both sets of
            synthetics MUST share the same `cmtsolution` and `stations` metadata

        .. note::
            Used for defining user-input waveforms data

        :type st: obspy.core.stream.Stream
        :param st: Stream objects containing observed time series to be plotted
            on the record section. Can contain any number of traces
        :type st_syn: obspy.core.stream.Stream
        :param st_syn: Stream objects containing synthetic time series to be
            plotted on the record section. Must be same length as `st`

        .. note::
            Waveform plotting organization parameters

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
                vertical spacing. Smallest distances at the top of the figure.
                 Requires `distance_units`
            - 'abs_distance': absolute vertical spacing of waveforms defined by
                source-receiver distance. Smallest distance at the top of the
                figure. Requires `distance_units`
            - 'abs_azimuth': absolute vertical spacing of waveforms defined
                by source-receiver azimuth (deg).
            - 'abs_backazimuth': absolute vertical spacing of waveforms by
                source-receiver backazimuth (deg).
            - '*_r': Add a '_r' to any of the values about to REVERSE the sort,
                e.g., 'alphabetical_r' sort will go Z->A
        :type scale_by: str
        :param scale_by: scale amplitude of waveforms by available:

            - None: Not set, no amplitude scaling, waveforms shown raw
            - 'normalize': scale each trace by the maximum amplitude,
                i.e., > a /= max(abs(a))  # where 'a' is time series amplitudes
            - 'global_norm': scale by the largest amplitude to be displayed on
                the screen. Will not consider waveforms which have been
                excluded on other basis (e.g., wrong component)
                i.e., > st[i].max /= max([max(abs(tr.data)) for tr in st])
            - 'geometric_spreading': scale amplitudes by expected reduction
                through geometric spreading. Related parameters are:

                - `geometric_spreading_factor`  
                - `geometric_spreading_k_val`  
                - `geometric_spreading_exclude`  
                - `geometric_spreading_ymax`  A

                Equation is A(d) = k / sin(d) ** f

                Where A(d) is the amplitude reduction factor as a function of
                distnace, d. 'k' is the `geometric_spreading_k_val` and 'f' is
                the `geometric_spreading_factor`.
                'k' is calculated automatically if not given.
        :type time_shift_s: float OR list of float OR str
        :param time_shift_s: apply static time shift to waveforms, two options:

            1. float (e.g., -10.2), will shift ALL waveforms by
                that number (i.e., -10.2 second time shift applied)
            2. list (e.g., [5., -2., ... 11.2]), will apply individual time
                shifts to EACH trace in the stream. The length of this list MUST
                match the number of traces in your input stream.
            3. str: apply time shift based on a theoretical TauP phase arrival
                if available in the SAC header. These should have been appended
                by PySEP during data download. If no value is available in the
                SAC header, defaults to 0. This may have unintended consequences
                so you should manually check that all values are okay.
                Available options are:
                - 'first_arrival_time': shift based on earliest phase arrival
                - 'p_arrival_time': shift based on earliest P phase arrival
                - 's_arrival_time': shift based on earliest S phase arrival
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
        :type azimuth_start_deg: float
        :param azimuth_start_deg: If sorting by azimuth, this defines the
            azimuthal angle for the waveform at the top of the figure.
            Set to 0 for default behavior
        :type distance_units: str
        :param distance_units: Y-axis units when sorting by epicentral distance
            'km': kilometers on the sphere
            'deg': degrees on the sphere
            'km_utm': kilometers on flat plane, UTM coordinate system

        .. note::
            Data processing parameters

        :type preprocess: str
        :param preprocess: choose which data to preprocess, options are:

            - None: do not run preprocessing step, including filter (Default)
            - 'st': process waveforms in `st`
            - 'st_syn': process waveforms in `st_syn`. st still must be given
            - 'both': process waveforms in both `st` and `st_syn`
        :type min_period_s: float
        :param min_period_s: minimum filter period in seconds
        :type max_period_s: float
        :param max_period_s: maximum filter period in seconds
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

        .. note::
            Geometric spreading parameters, used for amplitude scaling

        :type geometric_spreading_factor: float
        :param geometric_spreading_factor: factor to scale amplitudes by
            predicting the expected geometric spreading amplitude reduction and
            correcting for this factor.
            Related optional parameter: `geometric_spreading_k_val`.
            For Rayleigh waves, `geometric_spreading_factor` == 0.5 (default)
        :type geometric_spreading_k_val: float
        :param geometric_spreading_k_val: Optional constant scaling value used
            to scale the geometric spreading factor equation. If not given,
            calculated automatically using max amplitudes
            Value should be between 0.5 and 1.0 for regional surface waves.
        :type geometric_spreading_exclude: list
        :param geometric_spreading_exclude: a list of station names that
            should match the input stations. Used to exclude stations from the
            automatic caluclation of the geometric
        :type geometric_spreading_ymax: float
        :param geometric_spreading_ymax: Optional value for geometric spreading
            plot. Sets the max y-value on the plot. If not set, defaults to
            whatever the peak y-value plotted is.
        :type geometric_spreading_save: str
        :param geometric_spreading_save: file id to save separate geometric
            spreading scatter plot iff `scale_by`=='geometric_spreading'. If
            NoneType, will not save. By default, turned OFF

        .. note::
            Figure generation control parameters

        :type max_traces_per_rs: int
        :param max_traces_per_rs: maximum number of traces to show on a single
            record section plot. Defaults to all traces in the Stream
        :type xlim_s: list of float
        :param xlim_s: [start, stop] in units of time, seconds, to set the
            xlimits of the figure
        :type y_axis_spacing: float
        :param y_axis_spacing: spacing between adjacent seismograms applied to
            Y-axis on relative (not absolute) scales. Defaults to 1.
        :type y_label_loc: str
        :param y_label_loc: Location to place waveform labels on the y-axis

            - 'default': auto choose the best location based on `sort_by`
            - 'y_axis': Replace tick labels on the y-axis (left side of figure),
                This won't work if using absolute sorting and will be over-
                written by 'default'
            - 'y_axis_abs': For absolute y-axis only. waveform labels plotted
               on the left side outside border, with y-axis labels overlapping
               (showing distance or azimuth)
            - 'y_axis_abs_right': For absolute y-axis only. waveform labels
               plotted on the right side outside border, with y-axis labels
               on the left side of the figure (showing distance or azimuth)
            - 'y_axis_right': Replace tick labels on the right side of the
                y-axis. This option won't work with absolute sorting
            - 'x_min': Place labels on top of the waveforms at the global min
                x-value on the figure
            - 'x_max': Place labels on top of the waveforms at the global max
                x-value on the figure
            - None: Don't plot any text labels
        :type tmarks: list of float
        :param tmarks: place vertical lines at given reference times. Used for 
            marking reference times such as the event origin, or phase arrival.
            Input as a list of times in units of seconds (where T=0 is the 
            event origin time). For example `tmarks`=[0, 100, 200] would set 
            vertical lines at 0s, 100s and 200s
        :type figsize: tuple of float
        :param figsize: size the of the figure, passed into plt.subplots()
        :type show: bool
        :param show: show the figure as a graphical output
        :type save: str
        :param save: path to save output figure, will create the parent
            directory if it doesn't exist. If None, will not save (default).

        .. note::
            Internal RecSec parameters

        :type overwrite: bool
        :param overwrite: if the path defined by `save` exists, will overwrite
            the existing figure
        :type log_level: str
        :param log_level: level of the internal logger. In order of ascending
            verbosity: 'CRITICAL', 'WARNING', 'INFO', 'DEBUG'.

        :raises AssertionError: if any parameters are set incorrectly
        """
        # Set the logger level before running anything
        logger.setLevel(log_level)

        # Check file path existence before we do any reading 
        if pysep_path is not None:
            assert(os.path.exists(pysep_path)), \
                    f"`pysep_path` given but does not exist: '{pysep_path}'"
        if syn_path is not None: 
            assert(os.path.exists(syn_path)), \
                    f"`syn_path` given but does not exist: '{syn_path}'"

        # Read files from path if provided
        if pysep_path is not None:
            if not synsyn:
                # Expecting to find SAC files labelled as such
                fids = glob(os.path.join(pysep_path, "*.sac"))
                if not fids:
                    # Or if legacy naming schema, assume that the sac files have
                    # a given file format: event_tag.NN.SSS..LL.CC.c
                    fids = glob(os.path.join(pysep_path, "*.*.*.*.*.?"))
                if fids:
                    logger.info(f"Reading {len(fids)} files from: {pysep_path}")
                    # Overwrite stream, so reading takes precedence
                    st = Stream()
                    for fid in fids:
                        st += read(fid)
            else:
                # User has told us that we want to read 'data' as synthetics,
                # useful for comparing e.g., a synthetic-synthetic inversion
                logger.debug("reading `pysep_path` expecting SPECFEM-generated "
                             "synthetic data")
                st = self._generate_synthetic_stream(syn_path=pysep_path,
                                                     source=cmtsolution,
                                                     stations=stations,
                                                     cartesian=cartesian)

        # Read in SPECFEM generated synthetics and generate SAC headed streams
        if syn_path is not None: 
            assert(os.path.exists(syn_path)), \
                    f"`syn_path` given but does not exist: '{syn_path}'"
            st_syn = self._generate_synthetic_stream(syn_path=syn_path,
                                                     source=cmtsolution,
                                                     stations=stations,
                                                     cartesian=cartesian)

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
        self.geometric_spreading_exclude = geometric_spreading_exclude or []
        self.geometric_spreading_ymax = geometric_spreading_ymax
        self.geometric_spreading_save = geometric_spreading_save

        # Time shift parameters
        self.move_out = move_out
        # float check incase command line input provides as a string, and also
        # to ensure int -> float
        try:
            self.time_shift_s = float(time_shift_s)
        except (TypeError, ValueError):
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
        self.tmarks = tmarks
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
        self.max_amplitudes_syn = []
        self.amplitude_scaling = []
        self.amplitude_scaling_syn = []
        self.y_axis = []
        self.xlim = []  # unit: samples
        self.xlim_syn = []
        self.sorted_idx = []

    def _generate_synthetic_stream(self, syn_path, source, stations,
                                   cartesian=False, fmt="*.*.*.sem?*"):
        """
        Convenience fucntion to read in synthetic seismograms from SPECFEM2D,
        SPECFEM3D or SPECFEM3D_GLOBE. Can be used to read in both `st` and
        `st_syn`

        :type syn_path: str
        :param syn_path: full path to directory containing synthetic
            seismograms
        :type source: str
        :param source: path to source file which defined the source that
            generated the synthetics. Acceptable values are CMTSOLUTION (from
            SPECFEM3D/GLOBE), and SOURCE (from SPECFEM2D)
        :type stations: str
        :param stations: full path to STATIONS file used to define the station
            coordinates. Format is dictated by SPECFEM
        :type fmt: str
        :param fmt: the expected filename format of the sythetics. Based on
            ASCII style files generated by SPECFEM. Defaults to '*??.*.*.sem?*'.
            Expected filename looks something like: 'NN.SSS.BXC.semd'
        :rtype: obspy.core.stream.Stream
        :return: Stream object with synthetic waveforms
        """
        assert(source is not None and os.path.exists(source)), (
            f"If `syn_path` is given, RecSec also requires `cmtsolution`"
        )
        assert(stations is not None and os.path.exists(stations)), (
            f"If `syn_path` is given, RecSec also requires `stations`"
        )
        fids = glob(os.path.join(syn_path, fmt))
        if fids:
            logger.info(f"Reading {len(fids)} synthetics from: {syn_path}")
            st_syn = Stream()
            for fid in fids:
                logger.debug(fid)
                if not cartesian:
                    st_syn += read_sem(fid=fid, source=source,
                                       stations=stations)
                else:
                    # If we are using SPECFEM2D synthetics, trying to read
                    # the SOURCE file will
                    st_syn += read_sem_cartesian(fid=fid, source=source,
                                                 stations=stations)
        return st_syn

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

        # Make sure stream and synthetic stream have the same length. If they
        # don't, subset so that they do.
        if self.st_syn is not None:
            self.st, self.st_syn = subset_streams(self.st, self.st_syn)

            if len(self.st) != len(self.st_syn):
                err.st_syn = f"length must match `st` (which is {len(self.st)})"

        # Check the `sort_by` sorting parameter options
        acceptable_sort_by = ["default", "azimuth", "backazimuth",
                              "distance", "alphabetical", "abs_azimuth",
                              "abs_distance", "abs_backazimuth"]
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
            acceptable_time_shift_s = [float, list, str]
            if type(self.time_shift_s) not in acceptable_time_shift_s:
                err.time_shift_s = f"must be in {acceptable_time_shift_s}"
            if isinstance(self.time_shift_s, list) and \
                    len(self.time_shift_s) != len(self.st):
                err.time_shift_s = f"must be list of length {len(self.st)}"
            elif isinstance(self.time_shift_s, str):
                if self.time_shift_s not in SACDICT:
                    err.time_shift_s = f"must be {SACDICT} ending w/ '_time'"

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

        if self.min_period_s is not None or self.max_period_s is not None:
            assert(self.preprocess is not None), \
                f"Setting filter bounds requires `preprocess` flag to be set"

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
                                  "x_max", "y_axis_abs", "y_axis_abs_right",
                                  None]
        if self.y_label_loc not in acceptable_y_label_loc:
            err.y_label_loc = f"must be in {acceptable_y_label_loc}"
        if "abs" in self.sort_by and "y_axis" in self.sort_by:
            err.y_label_loc = (f"absolute sorting means 'y_axis' label loc is" 
                               f"not available")

        if len(self.figsize) != 2:
            err.figsize = "must be tuple defining (horizontal, vertical) extent"

        if self.save:
            if os.path.exists(self.save) and not self.overwrite:
                err.save = (f"path {self.save} already exists. Use "
                            f"'--overwrite' flag to save over existing figures"
                            )
            else:
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
        logger.info(f"determining if any stations/channels should be skipped")
        skip_idx = []
        for idx in self.idx:
            tr = self.st[idx]
            # Component-wise removal
            if tr.stats.component not in self.components:
                logger.debug(f"skip '{tr.get_id()}' for non-matching component")
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
        # Needs to be run before getting xlims so that we know the time offset
        self.get_time_offsets()
        self.xlim = self.get_xlims()
        # Get xlims synthetic waveforms which may have different start/end times
        if self.st_syn is not None:
            self.xlim_syn = self.get_xlims(st=self.st_syn)

        # Max amplitudes should be RELATIVE to what were showing (e.g., if
        # zoomed in on the P-wave, max amplitude should NOT be the surface wave)
        for tr, xlim in zip(self.st, self.xlim):
            start, stop = xlim
            self.max_amplitudes = np.append(self.max_amplitudes,
                                            max(abs(tr.data[start:stop])))
        self.max_amplitudes = np.array(self.max_amplitudes)

        # Max amplitudes will be DIFFERENT for synthetics, which affects 
        # normalizations and thus requires its own array
        if self.st_syn is not None:
            for tr, xlim in zip(self.st_syn, self.xlim):
                start, stop = xlim
                self.max_amplitudes_syn = np.append(
                        self.max_amplitudes_syn, max(abs(tr.data[start:stop]))
                        )
        self.max_amplitudes_syn = np.array(self.max_amplitudes_syn)

        # Figure out which indices we'll be plotting
        sorted_idx = self.get_sorted_idx()
        skip_idx = self.get_skip_idx()
        # Remove skip indexes from sorted index to get the final ordered list
        self.sorted_idx = np.array([_ for _ in sorted_idx if _ not in skip_idx])

        # Figure out how to manipulate each of the traces in the Stream
        self.y_axis = self.get_y_axis_positions()
        self.amplitude_scaling = self.get_amplitude_scaling(_choice="st")
        if self.st_syn:
            self.amplitude_scaling_syn = \
                    self.get_amplitude_scaling(_choice="st_syn")

    def get_xlims(self, st=None):
        """
        The x-limits of each trace depend on the overall time shift (either 
        static or applied through move out), as well as the sampling rate of
        each trace (which can vary). Retrieve an index-dependent list of
        x-limits which can be used to truncate the time series during plotting.

        .. note::
            Requires that get_time_shifts() has already been run

        :type st: obspy.core.stream.Stream
        :param st: stream object to get xlims for. By default this is the 'data'
            stored in `st` but it can also be given `st_syn` to get synthetic
            x limits which may differ
        :rtype: np.array
        :return: an array of tuples defining the start and stop indices for EACH
            trace to be used during plotting. Already includes time shift 
            information so xlim can be applied DIRECTLY to the time shifted data
        """
        if st is None:
            st = self.st

        xlim = []
        if self.xlim_s is None:
            # None's will index the entire trace
            xlim = np.array([(None, None) for _ in range(len(st))])
        else:
            # Looping to allow for delta varying among traces,
            # AND apply the time shift so that indices can be used directly in
            # the plotting function
            for tr, tshift in zip(st, self.time_shift_s):
                start, stop = [int(_/tr.stats.delta) for _ in self.xlim_s]
                sshift = int(tshift / tr.stats.delta)  # unit: samples
                # These indices define the index of the user-chosen timestamp
                start_index = start - sshift
                end_index = stop - sshift
                # Shift by the total amount that the zero pad adjusted starttime
                if self.zero_pad_s:
                    zero_pad_index = int(self.zero_pad_s[0]/tr.stats.delta)
                    start_index += zero_pad_index
                    end_index += zero_pad_index

                # Address time shift introduced by traces which do not start
                # at the origin time
                if hasattr(tr.stats, "time_offset"):
                    start_index += abs(int(tr.stats.time_offset/tr.stats.delta))
                    end_index += abs(int(tr.stats.time_offset/tr.stats.delta))

                # When setting xlimits, cannot acces out of bounds (negative 
                # indices or greater than max). This might happen when User 
                # asks for `xlim_s` that is outside data bounds. In this case
                # we just plot up to the end of the data
                if start_index < 0:
                    logger.warning(f"trying to access negative time axis index "
                                   f"for xlimit of {tr.get_id()}, setting to 0")
                    start_index = 0
                if end_index > len(tr.data):
                    logger.warning(f"trying to access past time axis index "
                                   f"for xlimit of {tr.get_id()}, setting to "
                                   f"length of data trace")
                    end_index = len(tr.data)

                xlim.append((start_index, end_index))

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

    def get_time_offsets(self):
        """
        Find time shift to data constituting difference between trace start 
        time and event origin time. Both synthetics and data should have the 
        correct timing. Time offsets will be used during plotting to make 
        sure all traces have the same T=0

        .. note::
            Appends time offset directly to the stats header, overwriting any
            value that might have been there already
        """
        event_origin_time = origin_time_from_sac_header(self.st[0].stats.sac)

        logger.info(f"calculating starttime offsets from event origin time "
                    f"{event_origin_time}")

        for tr in self.st:
            tr.stats.time_offset = tr.stats.starttime - event_origin_time
        if self.st_syn is not None:
            for tr in self.st_syn:
                tr.stats.time_offset = tr.stats.starttime - event_origin_time

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
            if isinstance(self.time_shift_s, float):
                logger.info(f"apply constant time shift {self.time_shift_s}s")
                time_shift_arr += self.time_shift_s
            # User input an array which should have already been checked for len
            elif isinstance(self.time_shift_s, list):
                logger.info(f"apply user-defined array time shift values")
                time_shift_arr = self.time_shift_s
            # Allow shifting by a given phase arrival in the SAC header
            elif isinstance(self.time_shift_s, str):
                sac_key = SACDICT[self.time_shift_s]
                logger.info(f"apply time shift by {self.time_shift_s}")
                time_shift_arr = [-1 * tr.stats.sac[sac_key] for tr in self.st]
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
                ((tr.stats.sac.stlo - tr.stats.sac.evlo) ** 2) +
                ((tr.stats.sac.stla - tr.stats.sac.evla) ** 2)
            )
            dist = kilometers2degrees(dist_deg)  # units: km

        return dist, az, baz

    def get_amplitude_scaling(self, _choice="st"):
        """
        Scale the amplitudes of all the waveforms by producing a Stream
        dependent scale factor based on user choice. It is expected that the
        output array will be DIVIDED by the data arrays:

        i.e., st[i].data /= self.amplitude_scaling[i]

        .. note::
            Needs to be run AFTER preprocessing because filtering etc. will
            affect the final amplitudes of the waveforms

        :type _choice: str
        :param _choice: Internal choice only required for 'normalize' scaling. 
            `st` for amplitude scaling w.r.t data array, `st_syn` 
            for amplitude scaling w.r.t synthetics
        :rtype: np.array
        :return: an array corresponding to the Stream indexes which provides
            a per-trace scaling coefficient
        """
        logger.info(f"determining amplitude scaling w.r.t {_choice} with: "
                    f"'{self.scale_by}'")

        # Don't scale by anything
        if self.scale_by is None:
            amp_scaling = np.ones(len(self.st))
        # Scale by the max amplitude of each trace
        elif self.scale_by == "normalize":
            if _choice == "st":
                amp_scaling = self.max_amplitudes
            elif _choice == "st_syn":
                amp_scaling = self.max_amplitudes_syn
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
            if self.max_amplitudes_syn.any():
                global_max = max(global_max, 
                                 self.max_amplitudes_syn[self.sorted_idx].max()
                                 )
            amp_scaling = np.ones(len(self.st)) * global_max
        # Scale by the theoretical geometrical spreading factor
        elif self.scale_by == "geometric_spreading":
            amp_scaling = self.calculate_geometric_spreading()

        # Apply manual scale factor if provided, default value is 1 so nothing
        # Divide because the amplitude scale divides the data array, which means
        # `amplitude_scale_factor` will be MULTIPLIED by the data array
        amp_scaling /= self.amplitude_scale_factor

        return amp_scaling

    def calculate_geometric_spreading(self):
        """
        Stations with larger source-receiver distances will have their amplitude
        scaled by a larger value.

        For information on geometric spreading, see Stein and Wysession,
        Section 4.3.4, Eq 20 (for Rayleigh waves, geometrical spreading
        factor = 0.5). For our purposes, this will fold the attenuation factor
        into the same single factor that accounts for geometric spreading.

        Equation is for amplitude reduction 'A' as a factor of distance 'd':

            A(d) = k / sin(d) ** f

        where 'k' is the `geometric_spreading_k_val` and 'f' is the
        `geometric_spreading_factor`. 'k' is calculated automatically if not
        given by the User. In the code, A(d) is represented by `w_vector`

        .. note::
            This does not take into account the variation in amplitude from
            focal mechanism regions

        TODO
            - Look in Stein and Wysession and figure out vector names
            - Allow ignoring specific components, currently only allowed to
                exclude on a per-station basis

        :rtype: list
        :return: scale factor per trace in the stream based on theoretical
            geometrical spreading factor. This is meant to be MULTIPLIED by the
            data arrays
        """
        logger.info("using geometrical spreading factor: "
                    f"{self.geometric_spreading_factor}")

        # Ignore any amplitudes with 0 max amplitude, which will throw off calc
        for i, max_amp in enumerate(self.max_amplitudes):
            if max_amp == 0:
                station_name = self.st[i].get_id().split(".")[1]
                if station_name not in self.geometric_spreading_exclude:
                    self.geometric_spreading_exclude.append(station_name)
                    logger.warning(f"{station_name} has 0 max amplitude which "
                                   f"will negatively affect geometric "
                                   f"spreading calculation. excluding from list"
                                   )

        # To selectively remove stations from this calculation, we will need
        # to gather a list of indexes to skip which can be applied to lists
        indices = self.sorted_idx  # already some stations have been removed
        for station in set(self.geometric_spreading_exclude):
            # Matching station names to get indices
            remove_idx = [i for i, id_ in enumerate(self.station_ids)
                          if station in id_]
            if remove_idx:
                logger.debug(f"remove '{station}' from geometric spreading eq.")
            indices = [i for i in indices if i not in remove_idx]
        logger.info(f"excluded {len(self.sorted_idx) - len(indices)} traces "
                    f"from geometric spreading calculation")

        # Make sure distances are in units degrees
        if "km" in self.distance_units:
            distances = kilometers2degrees(self.distances)
        else:
            distances = self.distances

        # Create a sinusoidal function based on distances in degrees
        sin_delta = np.sin(np.array(distances) / (180 / np.pi))

        # Use or calculate the K valuse constant scaling factor
        if self.geometric_spreading_k_val is None:
            k_vector = self.max_amplitudes[indices] * \
                       (sin_delta[indices] ** self.geometric_spreading_factor)
            # OVERWRITING internal parameter here, which will be used elsewhere
            self.geometric_spreading_k_val = np.median(k_vector)

        logger.info(f"geometric spreading `k` vector = "
                    f"{self.geometric_spreading_k_val}")

        # This is the amplitude scaling that we are after, defined for ALL stas
        w_vector = (self.geometric_spreading_k_val /
                    (sin_delta ** self.geometric_spreading_factor))

        # Plot the geometric spreading figure
        if self.geometric_spreading_save is not None:
            # Curate station names so that it is just 'STA.COMP'
            station_ids = [f"{_.split('.')[1]}.{_.split('.')[-1][-1]}" for _ in
                           self.station_ids]

            plot_geometric_spreading(
                distances=distances, max_amplitudes=self.max_amplitudes,
                geometric_spreading_factor=self.geometric_spreading_factor,
                geometric_k_value=self.geometric_spreading_k_val,
                station_ids=station_ids, units=self.distance_units,
                include=indices, ymax=self.geometric_spreading_ymax,
                save=self.geometric_spreading_save
            )

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

    def get_x_axis_tick_values(self):
        """
        Determine, based on the length of the plotted traces, how often tick
        marks should be applied so that they are spaced as aesthetically
        pleasing values.
        """
        # Get a rough estimate for total length of time of the trace in seconds
        if self.xlim_s:
            rs_len = self.xlim_s[1] - self.xlim_s[0]
        else:
            # Assuming that all waveforms are trimmed to the same length
            rs_len = self.st[0].stats.endtime - self.st[0].stats.starttime

        # Find order of magnitude of the length to give a rough estimate
        oom = np.floor(np.log10(rs_len))
        xtick_major = 10 ** oom
        xtick_minor = xtick_major / 2

        return xtick_minor, xtick_major

    def process_st(self):
        """
        Preprocess the Stream with optional filtering in place.

        .. note::
            Data in memory will be irretrievably altered by running preprocess.

        TODO Add feature to allow list-like periods to individually filter
            seismograms. At the moment we just apply a blanket filter.
        """
        if self.preprocess is None:
            logger.info("no preprocessing (including filtering) applied")
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
            # Max period only == high-pass
            if self.max_period_s is not None and self.min_period_s is None:
                logger.info(f"apply highpass filter w/ cutoff "
                            f"{1/self.max_period_s}")
                st.filter("highpass", freq=1/self.max_period_s, zerophase=True)
            # Min period only == low-pass
            elif self.min_period_s is not None and self.max_period_s is None:
                logger.info(f"apply lowpass filter w/ cutoff "
                            f"{1/self.min_period_s}")
                st.filter("lowpass", freq=1/self.min_period_s, zerophase=True)
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
        assert(nwav > 0), (
            f"no waveforms available for plotting. Check skip criteria if you "
            f"think this is an issue"
        )

        # Do a text output of station information so the user can check
        # that the plot is doing things correctly
        logger.debug("plotting line check starting from bottom (y=0)")
        logger.debug("\nIDX\tY\t\tID\tDIST\tAZ\tBAZ\tTSHIFT\tTOFFSET\tYABSMAX")
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

        # Plot vertical bars at given reference times
        if self.tmarks:
            self._plot_tmarks()

        # Change the aesthetic look of the figure, should be run before other
        # set functions as they may overwrite what is done here
        _xtick_minor, _xtick_major = self.get_x_axis_tick_values()
        # Use kwarg values to avoid double inputs of the same parameter
        for name, val in zip(["xtick_minor", "xtick_major"],
                             [_xtick_minor, _xtick_major]):
            if name not in self.kwargs:
                self.kwargs[name] = val
        self.ax = set_plot_aesthetic(ax=self.ax, **self.kwargs)

        # Partition the figure by user-specified azimuth bins for relative
        # (back)azimuth sorting only (not absolute)
        if self.sort_by and ("azimuth" in self.sort_by) and \
                ("abs_" not in self.sort_by):
            self._plot_azimuth_bins(start=start, stop=stop)

        # Finalize the figure accoutrements
        self._plot_title(nwav=nwav)
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

        # Synthetics will already have a time offset stat from the 
        # 'read_sem' function which grabs T0 value from ASCII
        # Data will have a time offset relative to event origin time 
        if hasattr(tr.stats, "time_offset"):
            x += tr.stats.time_offset  
            toffset_str = f"{tr.stats.time_offset:4.2f}"
        else:
            toffset_str = "  None"

        # Flip the sign of the y-axis if we are doing normal absolute
        # sorting because we are flipping the y-axis in _plot_axes()
        if "abs_" in self.sort_by and "_r" not in self.sort_by:
            sign = -1
        else:
            sign = 1

        # Ampplitude scaling may change between observed and synthetic if e.g.,
        # we are doing trace-wise normalization
        if choice == "st":
            amplitude_scaling = self.amplitude_scaling
            max_amplitudes = self.max_amplitudes
        elif choice == "st_syn":
            amplitude_scaling = self.amplitude_scaling_syn
            max_amplitudes = self.max_amplitudes_syn

        y = sign * tr.data / amplitude_scaling[idx] + self.y_axis[y_index]

        # Truncate waveforms to get figure scaling correct. Make the distinction
        # between data and synthetics which may have different start and end 
        # times natively
        if choice == "st":
            start, stop = self.xlim[idx]
        elif choice == "st_syn":
            start, stop = self.xlim_syn[idx]

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
                   f"\t{toffset_str}"
                   f"\t{max_amplitudes[idx]:.2E}\n"
                   )

        # Retain some stats for global plot args
        self.stats.xmin.append(x.min())
        self.stats.xmax.append(x.max())
        self.stats.ymin.append(y.min())
        self.stats.ymax.append(y.max())

        return log_str

    def _plot_tmarks(self):
        """
        Plot vertical lines at given reference times based on user input values
        """
        c = self.kwargs.get("tmark_c", "r")
        lw = self.kwargs.get("tmark_lw", 1.5)
        ls = self.kwargs.get("tmark_ls", "-")
        alpha = self.kwargs.get("tmark_alpha", 0.75)
        z = self.kwargs.get("tmark_zorder", 5)

        for tmark in self.tmarks:
            plt.axvline(x=tmark, c=c, linewidth=lw, linestyle=ls, zorder=z,
                        alpha=alpha)

    def _plot_azimuth_bins(self, start=None, stop=None):
        """
        If plotting by azimuth, create visual bin separators so the user has
        a visual understanding of radiation patterns etc.

        :type start: int
        :param start: optional starting index to determine the bounds of
            the azimuth bins
        :type stop: int
        :param stop: optional stop index to determine the bounds of
            the azimuth bins
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

        # Cut down the azimuth bins to the closest values if we are doing
        # multi-page record sections so that we don't end up plotting too many
        if start is not None and stop is not None:
            min_az = self.azimuths[self.sorted_idx[start:stop]].min()
            max_az = self.azimuths[self.sorted_idx[start:stop]].max()
            # Find the closest azimuth bins
            idx_start = np.argmin(np.abs(azimuth_bins - min_az))
            idx_end = np.argmin(np.abs(azimuth_bins - max_az))
            azimuth_bins = azimuth_bins[idx_start:idx_end]

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
                # By default, place absolute sorted labels on y-axis
                if self.y_label_loc == "default":
                    loc = "y_axis_abs_right"
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

        # Reverse the y-axis if we are doing absolute y-axis so that smaller
        # values appear at the top, which is opposite of how the y-axis works.
        # !!! This also requires flipping the sign of the plotted data in
        # !!! _plot_trace() to deal with the axis flip.
        if "abs_" in self.sort_by and "_r" not in self.sort_by:
            self.ax.invert_yaxis()
        else:
            logger.info("user requests inverting y-axis with absolute "
                        "reverse sort")

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
            # Dynamically determine y-axis ticks based on the the total span
            # of the y-axis
            ymin, ymax = self.ax.get_ylim()
            dist_span = ymax - ymin
            oom = np.floor(np.log10(dist_span))
            _ytick_major = 10 ** oom
            _ytick_minor = _ytick_major / 2
            ytick_minor = self.kwargs.get("ytick_minor", _ytick_minor)
            ytick_major = self.kwargs.get("ytick_major", _ytick_major)
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
            - y_axis_abs: for absolute plotting, place waveform labels to the
                left side of the figure (outside border), co-existing with
                y-axis labels
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

        # Generate a 'template' y-axis format to help Users decipher labels
        if self.sort_by is not None and "backazimuth" in self.sort_by:
            az_str = "BAZ"
        else:
            az_str = "AZ"

        # Y_FMT will include time shift IF there are time shifts
        y_fmt = f"NET.STA.LOC.CHA|{az_str}|DIST"
        if self.time_shift_s.sum() != 0:
            y_fmt += "|TSHIFT"

        # Option 1: Replacing y-axis tick labels with waveform labels
        # Set the y-axis labels to the left side of the left figure border
        if loc == "y_axis":
            # For relative plotting (not abs_), replace y_tyick labels with
            # station information
            self.ax.set_yticks(self.y_axis[start:stop])
            self.ax.set_yticklabels(y_tick_labels)
            plt.text(-.01, .99, y_fmt, ha="right", va="top",
                     transform=self.ax.transAxes, fontsize=fontsize)
        # Set the y-axis labels to the right side of the right figure border
        elif loc == "y_axis_right":
            self.ax.set_yticks(self.y_axis[start:stop])
            self.ax.set_yticklabels(y_tick_labels)
            self.ax.yaxis.tick_right()
            self.ax.yaxis.set_label_position("right")
            plt.text(1.01, .99, y_fmt, ha="left", va="top",
                     transform=self.ax.transAxes, fontsize=fontsize)
        # Option 2: Plotting labels as text objects, separate from y-axis labels
        else:
            # 2a: Set the y-axis labels inside the figure border (xmin or xmax)
            # Trying to figure out where the min or max X value is on the plot
            if loc == "x_min":
                ha = "left"
                va = "top"
                func = min
                x_val = func(self.stats.xmin)
                plt.text(.01, .99, y_fmt, ha=ha, va=va,
                         transform=self.ax.transAxes, fontsize=fontsize)
            elif loc == "x_max":
                ha = "right"
                va = "top"
                func = max
                x_val = func(self.stats.xmax)
                plt.text(.99, .99, y_fmt, ha=ha, va=va,
                         transform=self.ax.transAxes, fontsize=fontsize)
            # 2b: Set the y-axis labels outside figure border, co-existing with
            # y-labels showing distance or azimuth
            elif loc == "y_axis_abs":
                ha = "right"
                va = "center"
                func = min
                x_val = func(self.stats.xmin)
                plt.text(0, .99, y_fmt, ha=ha, va=va,
                         transform=self.ax.transAxes, fontsize=fontsize)
            elif loc == "y_axis_abs_right":
                ha = "left"
                va = "center"
                func = max
                x_val = func(self.stats.xmax)
                plt.text(1., .99, y_fmt, ha=ha, va=va,
                         transform=self.ax.transAxes, fontsize=fontsize)

            if self.xlim_s is not None:
                x_val = func([func(self.xlim_s), x_val])

            # Plotting y-axis labels for absolute scales
            if len(self.y_axis) == len(self.st):
                for idx, s in zip(self.sorted_idx[start:stop], y_tick_labels):
                    plt.text(x=x_val, y=self.y_axis[idx], s=s, ha=ha, va=va,
                             c=c, fontsize=fontsize)
            # Plotting y-axis labels for relative scales
            elif len(self.y_axis) == len(y_tick_labels):
                for y, s in zip(self.y_axis, y_tick_labels):
                    plt.text(x=x_val, y=y, s=s, ha=ha, va=va, c=c,
                             fontsize=fontsize)

    def _plot_title(self, nwav=None):
        """
        Create the title of the plot based on event and station information
        Allow dynamic creation of title based on user input parameters

        :type nwav: int
        :param nwav: if using subset, the title needs to know how many waveforms
            it's showing on the page. self.plot() should tell it
        """
        title = self.kwargs.get("title", None)
        if title is not None:
            self.ax.set_title(title)
            return

        # Defines the number of waveforms plotted on a single page, allowing
        # for subsets per page
        if nwav is None:
            nwav = len(self.sorted_idx)

        # T0: Zero pad is either a list of length or None
        if self.zero_pad_s:
            _start, _end = self.zero_pad_s
        else:
            _start, _end = 0, 0
        # Origintime needs to be shifted by the zero pad offset value
        origintime = str(min([tr.stats.starttime + _start for tr in self.st]))
        # YYYY-MM-DDTHH:MM:SS.ssssssZ - > YYYY-MM-DD HH:MM:SS.ss
        origintime = origintime.replace("T", " ")[:-5]

        # HYPO: If only one event, show hypocenter information
        if self.stats.nevents == 1:
            sac = self.st[0].stats.sac
            # Accomodate unknown depth and magnitude
            if "evdp" in sac:
                evdp = f"{sac['evdp']:.2f}km"
            else:
                evdp = "None"
            if "mag" in sac:
                mag = f"M{sac['mag']:.2f}"
            else:
                mag = "None"

            # HYPO: lon, lat, depth, mag
            hypo = (
                f"\nHYPO: ({sac['evlo']:.2f}{DEG}, {sac['evla']:.2f}{DEG}), " 
                f"Z={evdp}, Mag={mag}"
            )
        else:
            hypo = ""

        # CMP: Get the unique components that have been plotted only
        cmp = "".join(np.unique([self.st[i].stats.component
                                 for i in self.sorted_idx]))

        # GEOSPREAD: if geometric spreading scaling, post equation variables
        if self.scale_by == "geometric_spreading":
            geospread = (f" (k={self.geometric_spreading_k_val:.2E}, "
                         f"f={self.geometric_spreading_factor:.2f})")
        else:
            geospread = ""

        # FILT: If preprocessing is turned OFF, no filtering was applied
        if self.preprocess is not None:
            filt = f"\nFILT: [{self.min_period_s}, {self.max_period_s}]s"
        else:
            filt = ""

        # MOVE_OUT: if move out was applied, post
        if self.move_out:
            move_out = f"\nMOVE_OUT: {self.move_out}{self.distance_units}/s"
        else:
            move_out = ""

        title = (
            f"T0: {origintime}{hypo}\n"
            f"NWAV: {nwav}; NEVT: {self.stats.nevents}; "
            f"NSTA: {self.stats.nstation}; COMP: {cmp}\n"
            f"SORT_BY: {self.sort_by}; "
            f"SCALE_BY: {self.scale_by}{geospread}"
            # The following title parts are optional depending on application
            f"{filt}{move_out}"
        )
        self.ax.set_title(title)

    def run(self):
        """
        Convenience run function to run the RecordSection workflow in order.
        No internal logic for breaking the figure up into multiple record
        sections. For that see main function `plotw_rs`.
        """
        self.process_st()
        self.get_parameters()
        self.plot()


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
        # formatter_class=argparse.RawTextHelpFormatter,
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
                vertical spacing. Smallest distances at top of figure. 
                Requires `distance_units`
            - 'abs_distance': absolute vertical spacing of waveforms defined by
                source-receiver distance (km). Smallest distances at top of 
                figure. Requires `distance_units`
            - 'abs_azimuth': absolute vertical spacing of waveforms defined
                by source-receiver azimuth (deg). Requires
                `azimuth_start_deg`
            - 'abs_backazimuth': absolute vertical spacing of waveforms by
                source-receiver backazimuth (deg).
            - '*_r': Add a '_r' to any of the values about to REVERSE the sort,
                e.g., alphabetical_r sort will go Z->A"
                """)
                        )
    parser.add_argument("--scale_by", default=None, type=str, nargs="?",
                        help=textwrap.dedent("""
            How to sort the Y-axis of the record section
            - None: Not set, no amplitude scaling, waveforms shown raw
            - 'normalize': scale each trace by the maximum amplitude,
                i.e., > a /= max(abs(a))  # where 'a' is time series amplitudes
            - 'global_norm': scale by the largest amplitude to be displayed on
                the screen. Will not consider waveforms which have been 
                excluded on other basis (e.g., wrong component)
                """)
                        )
    parser.add_argument("--geometric_spreading_factor", default=0.5, type=float,
                        help="Scale amplitudes by predicting the "
                             "expected geometric spreading amplitude reduction "
                             "and correcting for it. This defaults to 0.5. "
                             "Optional related parameter: "
                             "--geometric_spreading_k_val")
    parser.add_argument("--geometric_spreading_exclude", default=None,
                        nargs="+", type=str,
                        help="Comma separated list of stations to exclude from "
                             "the geometric spreading equation. e.g., "
                             "'STA1,STA2,STA3'")
    parser.add_argument("--geometric_spreading_ymax", default=None,
                        nargs="?", type=float,
                        help="Controls the y-max value on the geometric "
                             "spreading plot.")
    parser.add_argument("--time_shift_s", default=None, nargs="?",
                        help="Set a constant time shift in unit: seconds OR "
                             "shift by a given phase arrival in SAC header")
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
    parser.add_argument("--tmarks", default=None, type=float,
                        nargs="+", help="Plot vertical lines at given reference"
                                        "times [s]. Input as a list of values, "
                                        "e.g., --tmarks 0 100 200")
    parser.add_argument("--save", default="./record_section.png", type=str,
                        nargs="?",
                        help="Path to save the resulting record section fig")
    parser.add_argument("-o", "--overwrite", default=False, action="store_true",
                        help="overwrite existing figure if path exists")
    parser.add_argument("--cartesian", default=False, action="store_true",
                        help="Let RecSec know that your domain is defined in"
                             "Cartesian coordinates. Used for SPECFEM2D "
                             "and SPECFEM3D_Cartesian domains defined in XYZ.")
    parser.add_argument("--synsyn", default=False, action="store_true",
                        help="Let RecSec know that both `pysep_path` and "
                             "`syn_path` should be read in as SPECFEM-"
                             "generated synthetics.")
    parser.add_argument("-L", "--log_level", default="DEBUG", type=str,
                        help="verbosity of logger: 'WARNING', 'INFO', 'DEBUG'")

    # Keyword arguments can be passed directly to the argparser in the same 
    # format as the above kwargs (e.g., --linewidth 2), but they will not have 
    # help messages or type checking
    parsed, unknown = parser.parse_known_args()
    for arg in unknown:
        if arg.startswith(("-", "--")):
            parser.add_argument(arg.split("=")[0])

    return parser


def plotw_rs(*args, **kwargs):
    """
    Plot Waveform Record Section (main call function for RecordSection)

    Instantiates the RecordSection class, runs processing and parameter getting,
    and then plots record section. Contains additional logic for breaking up
    figures into multiple pages if requested by the user, while keeping sort
    order and waveform spacing consistent across multiple reord sections.

    .. note::

        See RecordSection.__init__() for acceptable args and kwargs
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
        # Iterate backwards through the range so that the order of pages is
        # more natural, i.e., plotting the record section from the top down,
        # rather than the bottom up
        rnge = np.arange(len(rs.sorted_idx), 0, -1 * rs.max_traces_per_rs)

        # When `max_traces_per_rs` is not an integer divisor of the number of
        # streams, the range will not hit zero, so we need to ensure it does
        if rnge[-1] != 0:
            rnge = np.append(rnge, 0)

        for i, stop in enumerate(rnge[:-1]):
            j = i + 1
            start = rnge[j]
            rs.plot(subset=[start, stop], page_num=j)

    _end = datetime.now()
    logger.info(f"finished record section in t={(_end - _start)}s")


def main():
    """
    Convenience 'main' function to play nice with entry scripts

    .. rubric::
        $ recsec -h
    """
    parser = parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    plotw_rs(**vars(parser.parse_args()))


if __name__ == "__main__":
    main()

