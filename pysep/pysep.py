"""
Python Seismogram Extraction and Processing

Download, pre-process, and organize seismic waveforms, event and station
metadata

TODO
    * flesh out the template config file
    * save intermediate files? or some form of checkpointing
"""
import argparse
import os
import sys
import yaml
import warnings
import llnl_db_client

from glob import glob
from pathlib import Path
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNBadRequestException
from obspy.core.event import Event, Origin, Magnitude

from pysep import logger
from pysep.utils.cap_sac import (append_sac_headers, write_cap_weights_files,
                                 format_sac_header_w_taup_traveltimes,
                                 format_sac_headers_post_rotation)
from pysep.utils.curtail import (curtail_by_station_distance_azimuth,
                                 quality_check_waveforms_before_processing,
                                 quality_check_waveforms_after_processing)
from pysep.utils.fmt import format_event_tag, format_event_tag_legacy, get_codes
from pysep.utils.io import read_yaml, read_event_file, write_stations_file
from pysep.utils.llnl import scale_llnl_waveform_amplitudes
from pysep.utils.process import (merge_and_trim_start_end_times, resample_data,
                                 format_streams_for_rotation, rotate_to_uvw,
                                 estimate_prefilter_corners,
                                 append_back_azimuth_to_stats)
from pysep.utils.plot import plot_source_receiver_map
from pysep.recsec import plotw_rs


class Pysep:
    """
    Download, preprocess, and save waveform data using ObsPy
    """
    def __init__(self, config_file=None, event_selection="default",
                 client="IRIS", origin_time=None, reference_time=None,
                 networks="*", stations="*", locations="*", channels="*",
                 event_latitude=None, event_longitude=None, event_depth_km=None,
                 event_magnitude=None, remove_response=True,
                 remove_clipped=False, water_level=60, detrend=True,
                 demean=True, taper_percentage=0, rotate=None,
                 pre_filt="default", mindistance=0, maxdistance=20E3,
                 minazimuth=0, maxazimuth=360,
                 minlatitude=None, minlongitude=None, maxlatitude=None,
                 maxlongitude=None, resample_freq=50, scale_factor=1,
                 phase_list=None, seconds_before_event=20,
                 seconds_after_event=20, seconds_before_ref=100,
                 seconds_after_ref=300, taup_model="ak135", output_unit="VEL",
                 user=None, password=None, client_debug=False,
                 log_level="DEBUG", timeout=600, write_files="all",
                 plot_files="all", llnl_db_path=None, output_dir=None,
                 overwrite=False, legacy_naming=False, overwrite_event_tag=None,
                 **kwargs):
        """
        Define a default set of parameters

        TODO
            * removed resample_TF, was this used?
            * load config FIRST and then set defaults, make all defaults None?

        :type client: str
        :param client: ObsPy FDSN client to query data from, e.g., IRIS, LLNL,
            NCEDC or any FDSN clients accepted by ObsPy
        :type network: str
        :param network: name or names of networks to query for, if names plural,
            must be a comma-separated list, i.e., "AK,AT,AV", else "AK".
            Wildcard okay
        :type station: str
        :param station: station name or names to query for, wildcard okay
        :type location: str
        :param location: locations name or names to query for, wildcard okay
        :type channel: str
        :param channel: channel name or names to query for, wildcard okay
        :type remove_response: bool
        :param remove_response: Remove instrument response or not
        :type remove_clipped: remove any clipped stations from final output
        :param remove_clipped: Check if waveforms are clipped, remove if so
        :type water_level: float
        :param water_level: water_level to apply during filtering for small vals
        :type detrend: bool
        :param detrend: apply simple detrend to data during preprocessing
        :type demean: bool
        :param demean: apply demeaning to data during preprocessing
        :type taper: float
        :param taper: apply a taper to the waveform with ObsPy taper, fraction
            between 0 and 1 as the percentage of the waveform to be tapered
            Applied generally used when data is noisy, e.g., HutchisonGhosh2016
            To get the same results as the default taper in SAC,
            use max_percentage=0.05 and leave type as hann.
            Tapering also happens while resampling (see util_write_cap.py)
        :type rotate: list of str
        :param rotate: choose how to rotate the waveform data. Can include
            the following options (order insensitive):
            * ZNE: Rotate from arbitrary components to North, East, Up
            * RTZ: Rotate from ZNE to Radial, Transverse, Up (requires ZNE)
            * UVW: Rotate from ZNE to orthogonal UVW orientation
        :type prefilt: str, tuple or NoneType
        :param prefilt: apply a pre-filter to the waveforms before deconvolving
            instrument response. Options are:
            * 'default': automatically calculate (f0, f1, f2, f3) based on the
                length of the waveform (dictating longest allowable period) and
                the sampling rate (dictating shortest allowable period). This is
                the default behavior.
            * NoneType: do not apply any pre-filtering
            * tuple of float: (f0, f1, f2, f3) define the corners of your pre
                filter in units of frequency (Hz)
        :type mindistance: float
        :param mindistance: get waveforms from stations starting from a minimum
            distance away  from event (km)
        :type maxdistance: float
        :param maxdistance: get waveforms from stations up to a maximum distance
            from event (km)
        :type minazimuth: float
        :param minazimuth: get waveforms from a station starting at a given
            minimum azimuth (deg) out to a final maximum azimuth (deg) where
            0 degrees points North from the event
        :type maxazimuth: float
        :param maxazimuth: get waveforms from a station starting at a given
            minimum azimuth (deg) out to a final maximum azimuth (deg) where
            0 degrees points North from the event
        :type minlatitude: float
        :param minlatitude: minimum latitude for bounding box to search for
            stations and events
        :type maxlatitude: float
        :param maxlatitude: maximum latitude for bounding box to search for
            stations and events
        :type minlongitude: float
        :param minlongitude: minimum longitude for bounding box to search for
            stations and events
        :type maxlongitude: float
        :param maxlongitude: maximum longitude for bounding box to search for
            stations and events
        :type resample_freq: float
        :param resample_freq: frequency to resample data at
        :type scale_factor: float
        :param scale_factor: scale data by a constant, for CAP use 10**2
            (to convert m/s to cm/s)
        :type phase_list: list of str
        :param phase_list: phase names to get ray information from TauP with.
            Defaults to direct arrivals 'P' and 'S'
        :type reference_time: str
        :param reference_time: Waveform origin time. Allows for a static time
            shift from the event origin time, e.g., if there are timing errors
            with relation to o time. If not given, defaults to event origin time
        :type seconds_before_ref: float
        :param seconds_before_ref: time before origintime to fetch waveform data
        :type seconds_after_ref: float
        :param seconds_after_ref: time after origintime to fetch waveform data
        :type taup_model: str
        :param taup_model: name of TauP model to use to calculate phase arrivals
        :type output_unit: str
        :param output_unit: the output format of the waveforms, something like
            'DISP', 'VEL', 'ACC'
        :type user: str
        :param user: User ID if IRIS embargoes data behind passwords
        :type password: str
        :param password: Password if IRIS embargoes data behind passwords
        :type overwrite: bool
        :param overwrite: overwrite an existing PySEP event directory if
            the current process is trying to write to it
        :type legacy_naming: bool
        :param legacy_naming: revert to old PySEP naming schema for event tag,
            which names the directory and output SAC files
        :type overwrite_event_tag: str
        :param overwrite_event_tag: option to allow the user to set their own
            event tag, rather than the automatically generated one
        """
        # Internal attribute but define first so that it sits at the top of
        # written config files
        self.event_tag = None
        self.config_file = config_file

        # ObsPy client-related parameters
        self.client = client
        self.client_debug = bool(client_debug)
        self.timeout = timeout
        self._user = user
        self._password = password
        self.taup_model = taup_model

        # Parameters related to event selection
        self.event_selection = event_selection
        try:
            self.origin_time = UTCDateTime(origin_time)
        except TypeError:
            self.origin_time = None
        self.seconds_before_event = seconds_before_event
        self.seconds_after_event = seconds_after_event

        # Optional: if User wants to define an event on their own.
        # `event_depth_km` and `event_magnitude` are also used for client query
        self.event_latitude = event_latitude
        self.event_longitude = event_longitude
        self.event_depth_km = event_depth_km
        self.event_magnitude = event_magnitude

        # Waveform and StationXML gathering parameters
        self.networks = networks
        self.stations = stations
        self.channels = channels
        self.locations = locations

        # Waveform collection parameters
        self.reference_time = reference_time or self.origin_time
        self.seconds_before_ref = seconds_before_ref
        self.seconds_after_ref = seconds_after_ref
        self.phase_list = phase_list

        # NOTE: This default is a UAF LUNGS system-specific database path
        self.llnl_db_path = (
                llnl_db_path or
                "/store/raw/LLNL/UCRL-MI-222502/westernus.wfdisc"
        )

        # Event and station search and curtailing criteria
        self.mindistance = mindistance
        self.maxdistance = maxdistance
        self.minazimuth = minazimuth
        self.maxazimuth = maxazimuth
        self.minlatitude = minlatitude
        self.maxlatitude = maxlatitude
        self.minlongitude = minlongitude
        self.maxlongitude = maxlongitude

        # Preprocessing flags
        self.demean = bool(demean)
        self.detrend = bool(detrend)
        self.taper_percentage = taper_percentage
        self.rotate = rotate or ["ZNE", "RTZ"]
        self.remove_response = bool(remove_response)
        self.output_unit = output_unit
        self.water_level = water_level
        self.pre_filt = pre_filt
        self.scale_factor = scale_factor
        self.resample_freq = resample_freq
        self.remove_clipped = bool(remove_clipped)

        # Program related parameters
        self.output_dir = output_dir or os.getcwd()
        self.write_files = write_files
        self.plot_files = plot_files
        self.log_level = log_level
        self._legacy_naming = legacy_naming
        self._overwrite = overwrite
        self._overwrite_event_tag = overwrite_event_tag

        # Internally filled attributes
        self.c = None
        self.st = None
        self.inv = None
        self.event = None

        # Allow the user to manipulate the logger during __init__
        if log_level is not None:
            logger.debug(f"`log_level` set to {log_level}")
            logger.setLevel(log_level)
        else:
            logger.disabled = True

    def check(self):
        """
        Check input parameter validity against expected Pysep behavior
        """
        self.origin_time = UTCDateTime(self.origin_time)
        if self.reference_time is None:
            self.reference_time = self.origin_time
        self.reference_time = UTCDateTime(self.reference_time)

        if self.event_selection == "search" or self.client == "LLNL":
            for par in [self.seconds_before_event, self.seconds_after_event]:
                assert(par is not None), (
                    "Event selection requires the following parameters: "
                    "`seconds_before_event` and `seconds_after_event`")
        elif self.event_selection == "default":
            for par in [self.event_latitude, self.event_longitude,
                        self.event_depth_km, self.event_magnitude]:
                assert(par is not None), (
                    "`event_selection`=='default' requires "
                    "`event_latitude`, `event_longitude`, `event_depth_km` "
                    "and `event_magnitude`")
        else:
            raise ValueError("`event_selection` must be one of the following: "
                             "'search' or 'default'")

        if self.client.upper() == "LLNL":
            assert(os.path.exists(self.llnl_db_path)), (
                f"`llnl_db_path` {self.llnl_db_path} does not exist but must "
                f"if `client`=='LLNL'. Please check path and try again"
            )

        if self.client.upper() == "NCEDC":
            assert("-" not in self.stations), (
                "`client`=='NCEDC' does not allow for '-' in station codes"
            )

        if self.networks == "*" and self.client.upper() == "IRIS":
            logger.warning("`networks`=='*' will search ALL networks, which "
                           "may take a long time depending on data availability"
                           )

        # Check preprocessing flags
        if self.remove_response:
            for par in [self.output_unit, self.water_level]:
                assert (par is not None), (
                    "`remove_response` requires parameters:"
                    "`output_unit`, `water_level` "
                )

        if not (0 <= self.minazimuth <= 360):
            _old_val = self.minazimuth
            self.minazimuth = self.minazimuth % 360
            logger.warning(f"0 <= `minazimuth` <= 360; "
                           f"{_old_val} -> {self.minazimuth}")

        if not (0 <= self.maxazimuth <= 360):
            _old_val = self.maxazimuth
            self.maxazimuth = self.maxazimuth % 360
            logger.warning(f"0 <= `maxazimuth` <= 360; "
                           f"{_old_val} -> {self.maxazimuth}")

        if self.rotate is not None:
            acceptable_rotations = {"RTZ", "UVW", "ENZ"}
            self.rotate = ["".join(sorted(val.upper())) for val in self.rotate]
            assert(set(self.rotate).issubset(acceptable_rotations)), (
                f"`rotate` must be a subset of: {acceptable_rotations}"
            )
            if "RTZ" in self.rotate:
                assert("UVW" not in self.rotate), (
                    f"rotate can only have one of the following: 'UVW', 'RTZ'"
                )
        else:
            self.rotate = []

        acceptable_units = ["DISP", "VEL", "ACC", "DEF"]
        self.output_unit = self.output_unit.upper()
        assert(self.output_unit in acceptable_units), (
            f"unnacceptable `output_unit` {self.output_unit}, must be in "
            f"{acceptable_units}")

        if self.pre_filt not in [None, "default"]:
            assert(len(self.pre_filt) == 4), (
                f"`pre_filt` must be a tuple of length 4, representing four "
                f"corner frequencies for a bandpass filter (f1, f2, f3, f4)"
            )

        # Enforce that `write_files` and `plot_files` are sets
        if self.write_files is None:
            self.write_files = {}
        else:
            try:
                self.write_files = {self.write_files}
            # TypeError thrown if we're trying to do {{*}}
            except TypeError:
                pass
            acceptable_write_files = self.write(_return_filenames=True)
            assert(self.write_files.issubset(acceptable_write_files)), (
                f"`write_files` must be a list of some or all of: "
                f"{acceptable_write_files}"
            )

        if self.plot_files is None:
            self.plot_files = {}
        else:
            try:
                self.plot_files = {self.plot_files}
            except TypeError:
                pass
            acceptable_plot_files = {"map", "record_section", "all"}
            assert(self.plot_files.issubset(acceptable_plot_files)), (
                f"`plot_files` must be a list of some or all of: "
                f"{acceptable_plot_files}"
            )

    def get_client(self):
        """
        Options to choose different Clients based on attribute `client` which
        will be used to gather waveforms and metadata

        :rtype: obspy.clients.fdsn.client.Client
        :return: Client used to gather waveforms and metadata
        """
        if self.client is None:
            return None

        # Lawrence Livermore Natinoal Lab internal waveform database
        if self.client.upper() == "LLNL":
            c = llnl_db_client.LLNLDBClient(self.llnl_db_path)
        # IRIS DMC PH5WS Station Web Service
        elif self.client.upper() == "PH5":
            c = Client(
                "http://service.iris.edu", debug=self.client_debug,
                service_mappings={
                    "station": "http://service.iris.edu/ph5ws/station/1",
                    "dataselect": "http://service.iris.edu/ph5ws/dataselect/1"
                }, user=self._user, password=self._password
            )
        # Default ObsPy FDSN webservice client
        else:
            c = Client(self.client, user=self._user, password=self._password,
                       debug=self.client_debug, timeout=self.timeout)

        return c

    def load(self, config_file=None, overwrite_event=True):
        """
        Overwrite default parameters using a YAML config file

        :type config_file: str
        :param config_file: YAML configuration file to load from
        :type overwrite_event: bool
        :param overwrite_event: overwrite event search parameters (origin time,
            lat, lon etc.) from the YAML config file. Defaults to True
        """
        if config_file is None:
            config_file = self.config_file
        if not overwrite_event:
            logger.info("will NOT overwrite event search parameters with "
                        "config file")

        if config_file is not None:
            logger.info(f"overwriting default parameters with config file: "
                        f"'{config_file}'")
            config = read_yaml(config_file)
            for key, val in config.items():
                if hasattr(self, key):
                    if not overwrite_event and \
                            key.startswith(("event_", "origin_time")):
                        continue
                    old_val = getattr(self, key)
                    if val != old_val:
                        logger.debug(f"{key}: {old_val} -> {val}")
                        setattr(self, key, val)

    def get_event(self):
        """
        Exposed API for grabbing event metadata depending on the
        `event_selection` choice.

        Options for `event_selection` are:
            'search': query FDSN with event parameters
            'default': create an event from scratch using user parameters
                or if 'client'=='LLNL', grab event from internal database

        :rtype: obspy.core.event.Event
        :return: Matching event given event criteria. If multiple events are
            returned with the query, returns the first in the catalog
        """
        if self.event_selection == "search":
            event = self._query_event_from_client()
        elif self.event_selection == "default":
            if self.client.upper() == "LLNL":
                event = self._get_event_from_llnl_catalog()
            else:
                event = self._create_event_from_scratch()
        else:
            event = None

        # This is how we will access the event info so print out and check
        lat = event.preferred_origin().latitude
        lon = event.preferred_origin().longitude
        depth_km = event.preferred_origin().depth * 1E-3
        otime = event.preferred_origin().time
        mag = event.preferred_magnitude().mag
        logger.info(f"event info summary - origin time: {otime}; "
                    f"lat={lat:.2f}; lon={lon:.2f}; depth[km]={depth_km:.2f}; "
                    f"magnitude={mag:.2f}")

        return event

    def _query_event_from_client(self, magnitude_buffer=0.1,
                                 depth_buffer_km=1.):
        """
        Retrieve an event catalog using ObsPy Client.get_events().
        Searches Client for a given origin time, location, depth (optional)
        and magnitude (optional).

        To use this, set attribute `event_selection`=='search'

        :type magnitude_buffer: float
        :param magnitude_buffer: if attribute `event_magnitude` is given,
            will search events for events with magnitude:
            event_magnitude +/- magnitude_buffer
        :type depth_buffer_km: float
        :param depth_buffer_km: if attribute `event_depth_km` is given,
            will search events for events with depth:
            event_depth_km +/- depth_buffer_km
        :rtype: obspy.core.event.Event
        :return: Matching event given event criteria. If multiple events are
            returned with the query, returns the first in the catalog
        """
        logger.info(f"getting event information with client {self.client}")

        # Allow user to guess magnitude of the event we're searching for
        if self.event_magnitude is not None:
            minmagnitude = self.event_magnitude - magnitude_buffer
            maxmagnitude = self.event_magnitude + magnitude_buffer
        else:
            logger.info("no `magnitude` specified, this will make the search "
                        "criteria more broad. If the returned catalog is too "
                        "large, consider setting a guess value for `magnitude`")
            minmagnitude = None
            maxmagnitude = None

        # Allow user to guess depth of the event we're searching for
        if self.event_depth_km is not None:
            mindepth = self.event_depth_km - depth_buffer_km
            maxdepth = self.event_depth_km + depth_buffer_km
        else:
            logger.info("no `depth` specified, this will make the search "
                        "criteria more broad. If the returned catalog is too "
                        "large, consider setting a guess value for `depth`")
            mindepth = None
            maxdepth = None

        cat = self.c.get_events(
            starttime=self.origin_time - self.seconds_before_event,
            endtime=self.origin_time + self.seconds_after_event,
            minlatitude=self.minlatitude, maxlatitude=self.maxlatitude,
            minlongitude=self.minlongitude, maxlongitude=self.maxlongitude,
            minmagnitude=minmagnitude, maxmagnitude=maxmagnitude,
            mindepth=mindepth, maxdepth=maxdepth, debug=self.client_debug,
        )

        logger.info(f"{len(cat)} matching events found for "
                    f"{self.seconds_before_event}s < {self.origin_time} "
                    f"< {self.seconds_after_event}, choosing first")

        event = cat[0]

        return event

    def _create_event_from_scratch(self):
        """
        Make a barebones event object based on user-defined parameters which
        will then be used to query for waveforms and StationXML data

        :rtype: obspy.core.event.Event
        :return: Event object with origin and magnitude information appended
        """
        logger.info("creating event metadata with user parameters")

        origin = Origin(latitude=self.event_latitude,
                        longitude=self.event_longitude,
                        depth=self.event_depth_km * 1E3,  # units: m
                        time=self.origin_time
                        )
        magnitude = Magnitude(mag=self.event_magnitude, magnitude_type="Mw")

        event = Event(origins=[origin], magnitudes=[magnitude])
        event.preferred_origin_id = origin.resource_id.id
        event.preferred_magnitude_id = magnitude.resource_id.id

        return event

    def _get_event_from_llnl_catalog(self):
        """
        Special getter function for Lawrence Livermore National Lab data
        LLNL database has a special client

        TODO Do we need more filtering in the catalog?

        :rtype: obspy.core.event.Event
        :return: Event information queried from LLNL database
        """
        logger.info("getting event information from LLNL database")

        cat = self.c.get_catalog()
        mintime_str = f"time > {self.origin_time - self.seconds_before_event}"
        maxtime_str = f"time < {self.origin_time + self.seconds_after_event}"

        cat = cat.filter(mintime_str, maxtime_str)
        logger.debug(f"{len(cat)} events in LLNL catalog, picking zeroth index")
        event = cat[0]

        # Apply preferred magnitude and origin if not already
        if not event.preferred_magnitude():
            event.preferred_magnitude_id = event.magnitudes[0].resource_id.id
        if not event.preferred_origin():
            event.preferred_origin_id = event.origins[0].resource_id.id

        return event

    def get_stations(self):
        """
        Exposed API for grabbing station metadata from client.
        Download station metadata using ObsPy get_stations() with a user-defined
        bounding box and for user-defined networks, stations etc.

        :rtype: obspy.core.inventory.Inventory
        :return: Station metadata queried from Client
        """
        logger.info(f"querying {self.client.upper()} for station metadata")
        if self.client.upper() == "LLNL":
            # LLNL DB client behaves differently than ObsPy clients
            inv = self.c.get_inventory()
            inv = inv.select(
                    network=self.networks, location=self.locations,
                    station=self.stations, channel=self.channels,
                    starttime=self.origin_time - self.seconds_before_ref,
                    endtime=self.origin_time + self.seconds_after_ref,
                    minlatitude=self.minlatitude, maxlatitude=self.maxlatitude,
                    minlongitude=self.minlongitude,
                    maxlongitude=self.maxlongitude
                    )
        else:
            # Standard behavior is to simply wrap get_stations()
            inv = self.c.get_stations(
                network=self.networks, location=self.locations,
                station=self.stations, channel=self.channels,
                starttime=self.origin_time - self.seconds_before_ref,
                endtime=self.origin_time + self.seconds_after_ref,
                minlatitude=self.minlatitude, maxlatitude=self.maxlatitude,
                minlongitude=self.minlongitude, maxlongitude=self.maxlongitude,
                level="response"
            )

        nnet = len(inv)
        ncha = len(inv.get_contents()["channels"])
        nsta = len(inv.get_contents()["stations"])

        logger.info(f"collected {ncha} channels from {nsta} stations in "
                    f"{nnet} networks from {self.client}")

        return inv

    def get_waveforms(self):
        """
        Exposed API for grabbing waveforms from client. Internal logic
        determines how waveforms are queried, but mainly it is controlled
        by the internal `inv` attribute detailing station information,
        and reference times for start and end times.

        .. note::
            We do not use the `minimumlength` variable so that we can figure
            out which stations have data gaps

        :rtype: obspy.core.stream.Stream
        :return: Stream of channel-separated waveforms
        """
        logger.info(f"querying client '{self.client.upper()}' for waveforms")

        # PH5 should not be queried in bulk as it is a small request
        if self.client.upper() == "PH5":
            st = self.c.get_waveforms(
                network=self.networks, location=self.locations,
                station=self.stations, channel=self.channels,
                starttime=self.origin_time - self.seconds_before_ref,
                endtime=self.origin_time + self.seconds_after_ref,
            )
        # LLNL DB has custom interface to client
        elif self.client.upper() == "LLNL":
            event_id = int(self.event.event_descriptions[0].text)
            st = self.c.get_waveforms_for_event(event_id)  # looks for integers
        else:
            st = self._bulk_query_waveforms_from_client()

        logger.info(f"{len(st)} waveforms returned after query")

        return st

    def _bulk_query_waveforms_from_client(self):
        """
        Make a bulk request query to the Client based on the internal `inv`
        attribute defining the available station metadata.

        TODO removed 'make_bulk_list_from_stalist' which used TauP to gather
            data between P and S arrivals. Do we need this?

        :rtype: obspy.core.stream.Stream
        :return: Stream of channel-separated waveforms
        """
        bulk = []
        t1 = self.reference_time - self.seconds_before_ref
        t2 = self.reference_time + self.seconds_after_ref
        for net in self.inv:
            for sta in net:
                # net sta loc cha t1 t2
                bulk.append((net.code, sta.code, "*", self.channels, t1, t2))

        try:
            logger.info(f"querying {len(bulk)} lines in bulk client request...")
            st = self.c.get_waveforms_bulk(bulk=bulk)
        except FDSNBadRequestException:
            logger.warning(f"client {self.client} returned no waveforms, "
                           f"please check your event and station "
                           f"parameters again and re-submit")
            sys.exit(-1)

        return st

    def curtail_stations(self):
        """
        Remove stations from `inv` based on station distance, azimuth, etc.

        .. note::
            One-function function currently, but we can expand curtailing here
            if need by

        :rtype: obspy.core.inventory.Inventory
        :return: station metadata that has been curtailed based on acceptable
            paramaters
        """
        inv = self.inv.copy()

        inv = curtail_by_station_distance_azimuth(
            event=self.event, inv=inv, mindistance=self.mindistance,
            maxdistance=self.maxdistance, minazimuth=self.minazimuth,
            maxazimuth=self.maxazimuth
        )

        return inv

    def preprocess(self):
        """
        Very simple preprocessing to remove response and apply a prefilter
        scale waveforms (if necessary) and clean up waveform time series

        :rtype: obspy.core.stream.Stream
        :return: a preprocessed stream with response removed, amplitude scaled
            (optional), and time series standardized
        """
        st_out = self.st.copy()
        if self.detrend:
            logger.info(f"applying linear detrend to all data")
            st_out.detrend("linear")
        if self.remove_response:
            logger.info(f"removing response, output units in: "
                        f"{self.output_unit}")
            if self.pre_filt is not None:
                logger.info(f"will apply pre-filter: {self.pre_filt}")
            if self.client.upper() == "LLNL":
                st_out = self._remove_response_llnl(st_out)
            else:
                for code in get_codes(st=self.st, choice="channel", suffix="?",
                                      up_to=True):
                    net, sta, loc, cha = code.split(".")
                    st_sta = st_out.select(network=net, station=sta,
                                           location=loc, channel=cha
                                           )
                    for tr in st_sta:
                        # Get trace-dependent pre-filtering if desired
                        if self.pre_filt == "default":
                            _pre_filt = estimate_prefilter_corners(tr)
                        else:
                            _pre_filt = self.pre_filt
                        try:
                            tr.remove_response(
                                inventory=self.inv,
                                water_level=self.water_level,
                                pre_filt=_pre_filt,
                                taper=bool(self.taper_percentage),
                                taper_fraction=self.taper_percentage,
                                zero_mean=self.demean,
                                output=self.output_unit
                                )
                        except ValueError as e:
                            logger.warning(f"can't remove response {code}: {e}"
                                           f"removing trace from stream")
                            st_out.remove(tr)
        if self.scale_factor:
            logger.info(f"applying amplitude scale factor: {self.scale_factor}")
            for tr in st_out:
                tr.data = tr.data * self.scale_factor
                tr.stats.sac["scale"] = self.scale_factor
        if self.client == "LLNL":
            # This won't do anything if we don't have any 'LL' network codes
            st_out = scale_llnl_waveform_amplitudes(st_out)

        # Apply pre-resample lowpass, resample waveforms, make contiguous
        st_out = merge_and_trim_start_end_times(st_out)
        st_out = resample_data(st_out, resample_freq=self.resample_freq)

        return st_out

    def _remove_response_llnl(self, st):
        """
        Remove response information from LLNL stations. This requires using
        the custom LLNL DB client. There are also some internal checks that
        need to be bypassed else they cause the program to crash
        """
        st_out = st.copy()

        # Mimicing the client check statements to remove stations that
        # have no response or matching time span
        for tr in st_out[:]:
            net = tr.stats.network
            sta = tr.stats.station
            cha = tr.stats.channel
            # Check that response info is available
            if (sta, cha) not in self.c.sensors:
                logger.warning(f"no response for {net}.{sta}.{cha}")
                st_out.remove(tr)
                continue
            # Check epoch times against the stream midpoint
            time = tr.stats.starttime + \
                   (tr.stats.endtime - tr.stats.starttime) / 2.
            for epoch in self.c.sensors[(sta, cha)]:
                if epoch.starttime <= time <= epoch.endtime:
                    break
                else:
                    logger.warning(f"wrong epoch for {net}.{sta}.{cha}")
                    try:
                        st_out.remove(tr)
                    except ValueError:
                        pass
                    continue

        self.c.remove_response(st_out, water_level=self.water_level,
                               output=self.output_unit,
                               pre_filt=self.pre_filt)

        return st_out

    def rotate_streams(self):
        """
        Rotate arbitrary three-component seismograms to desired orientation

        TODO Original code had a really complicated method for rotating,
            was that necesssary? Why was st.rotate() not acceptable?

        :rtype: obspy.core.stream.Stream
        :return: a stream that has been rotated to desired coordinate system
            with SAC headers that have been adjusted for the rotation
        """
        st_out = self.st.copy()
        st_out = format_streams_for_rotation(st_out)
        st_out = append_back_azimuth_to_stats(st=st_out, inv=self.inv,
                                              event=self.event)
        if "ZNE" in self.rotate:
            logger.info("rotating to components ZNE")
            st_out.rotate(method="->ZNE", inventory=self.inv)
        if "UVW" in self.rotate:
            logger.info("rotating to components UVW")
            st_out = rotate_to_uvw(st_out)
        elif "RTZ" in self.rotate:
            logger.info("rotating to components RTZ")
            st_out.rotate("NE->RT", inventory=self.inv)

        st_out = format_sac_headers_post_rotation(st_out)

        return st_out

    def write(self, write_files=None, _return_filenames=False,
              config_fid=None, stations_fid=None, inv_fid=None, event_fid=None,
              stream_fid=None, **kwargs):
        """
        Write out various files specifying information about the collected
        stations and waveforms.

        Options are:
            * weights_dist: write out 'weights.dat' file sorted by distance
            * weights_az: write out 'weights.dat' file sorted by azimuth
            * weights_code: write out 'weights.dat' file sorted by sta code
            * station_list: write a text file with station information
            * inv: write the inventory as a StationXML file
            * event: write the event as a QuakeML file
            * stream: write the stream as a single MSEED file
            * sac: write the stream as individual (per-channel) SAC files with
                the appropriate SAC header
            * config_file: write the current configuration as a YAML file
            * all: ignores all other options, writes everything listed above

        :type write_files: list of str
        :param write_files: list of files that should be written out, must
            match the acceptable list defined in the function or here in the
            docstring
        :type _return_filenames: bool
        :param _return_filenames: internal flag to not actually write anything
            but just return a list of acceptable filenames. This keeps all the
            file naming definitions in one function. This is only required by
            the check() function.
        :type config_fid: str
        :param config_fid: optional name for the configuration file name
            defaults to 'pysep_config.yaml'
        :type stations_fid: str
        :param stations_fid: optional name for the stations list file name
            defaults to 'station_list.txt'
        :type inv_fid: str
        :param inv_fid: optional name for saved ObsPy inventory object,
            defaults to 'inv.xml'
        :type event_fid: optional name for saved ObsPy Event object, defaults to
            'event.xml'
        :type stream_fid: optional name for saved ObsPy Stream miniseed object,
            defaults to 'stream.ms'
        """
        # This is defined here so that all these filenames are in one place,
        # but really this set is just required by check(), not by write()
        _acceptable_files = {"weights_az", "weights_dist", "weights_code",
                             "station_list", "inv", "event", "stream",
                             "config_file", "sac", "all"}
        if _return_filenames:
            return _acceptable_files

        # Allow the user to call write() with their own set of filenames if this
        # wasn't defined by the config or this is being scripted and they only
        # want certain files out at intermediate steps
        if write_files is None:
            write_files = self.write_files
        else:
            write_files = set(write_files)
            assert(write_files.issubset(_acceptable_files)), (
                f"`write_files` must be a list of some or all of: "
                f"{_acceptable_files}"
            )

        for weights_fid in ["weights_dist", "weights_az", "weights_code"]:
            if weights_fid in write_files or "all" in write_files:
                order_by = weights_fid.split("_")[1]
                write_cap_weights_files(st=self.st, order_by=order_by,
                                        path_out=self.output_dir)

        if "config_file" in write_files or "all" in write_files:
            logger.info("writing config YAML file")
            fid = os.path.join(self.output_dir,
                               config_fid or f"pysep_config.yaml")
            self.write_config(fid=fid)

        if "station_list" in write_files or "all" in write_files:
            fid = os.path.join(self.output_dir,
                               stations_fid or "stations_list.txt")
            logger.info("writing stations file")
            logger.debug(fid)
            write_stations_file(self.inv, self.event, fid)

        if "inv" in write_files or "all" in write_files:
            fid = os.path.join(self.output_dir, inv_fid or f"inv.xml")
            logger.info("writing inventory as StationXML")
            logger.debug(fid)
            self.inv.write(fid, format="STATIONXML")

        if "event" in write_files or "all" in write_files:
            fid = os.path.join(self.output_dir, event_fid or f"event.xml")
            logger.info("writing event as QuakeML")
            logger.debug(fid)
            self.event.write(fid, format="QuakeML")

        if "stream" in write_files or "all" in write_files:
            fid = os.path.join(self.output_dir, stream_fid or f"stream.ms")
            logger.info("writing waveform stream in MiniSEED")
            logger.debug(fid)
            with warnings.catch_warnings():
                # ignore the encoding warning that comes from ObsPy
                warnings.simplefilter("ignore")
                self.st.write(fid, format="MSEED")

        if "sac" in write_files or "all" in write_files:
            logger.info("writing each waveform trace in SAC format")
            # Old PySEP dumped all files into output dir. New PySEP makes subdir
            if self._legacy_naming:
                _output_dir = self.output_dir
            else:
                _output_dir = os.path.join(self.output_dir, "SAC")
            if not os.path.exists(_output_dir):
                os.makedirs(_output_dir)
            for tr in self.st:
                if self._legacy_naming:
                    # Legacy: e.g., 20000101000000.NN.SSS.LL.CC.c
                    _trace_id = f"{tr.get_id()[:-1]}.{tr.get_id()[-1].lower()}"
                    tag = f"{self.event_tag}.{_trace_id}"
                else:
                    # New style: e.g., 2000-01-01T000000.NN.SSS.LL.CCC.sac
                    tag = f"{self.event_tag}.{tr.get_id()}.sac"

                fid = os.path.join(_output_dir, tag)
                logger.debug(os.path.basename(fid))
                tr.write(fid, format="SAC")

    def write_config(self, fid=None):
        """
        Write a YAML config file based on the internal `Pysep` attributes.
        Remove a few internal attributes (those containing data) before writing
        and also change types on a few to keep the output file simple but
        also re-usable for repeat queries.

        :type fid: str
        :param fid: name of the file to write. defaults to config.yaml
        """
        if fid is None:
            fid = f"pysep_config.yaml"
        logger.debug(fid)
        dict_out = vars(self)

        # Drop hidden variables
        dict_out = {key: val for key, val in dict_out.items()
                    if not key.startswith("_")}
        # Internal attributes that don't need to go into the written config
        attr_remove_list = ["st", "event", "inv", "c", "write_files",
                            "plot_files", "output_dir"]

        if self.client.upper() != "LLNL":
            attr_remove_list.append("llnl_db_path")  # not important unless LLNL
        for key in attr_remove_list:
            del(dict_out[key])
        # Write times in as strings not UTCDateTime
        for key in ["origin_time", "reference_time"]:
            val = dict_out[key]
            if isinstance(val, UTCDateTime):
                dict_out[key] = str(val)

        with open(fid, "w") as f:
            yaml.dump(dict_out, f, default_flow_style=False, sort_keys=False)

    def plot(self):
        """
        Plot map and record section if requested

        TODO improve source receiver map plotting
        """
        if "map" in self.plot_files or "all" in self.plot_files:
            logger.info("plotting source receiver map")
            fid = os.path.join(self.output_dir, f"station_map.png")
            plot_source_receiver_map(self.inv, self.event, save=fid)

        if "record_section" in self.plot_files or "all" in self.plot_files:
            fid = os.path.join(self.output_dir, f"record_section.png")
            # Default settings to create a general record section
            plotw_rs(st=self.st, sort_by="distance_r", scale_by="normalize",
                     overwrite=True, show=False, save=fid)

    def _event_tag_and_output_dir(self):
        """
        Convenience function to establish and naming schema for files and
        directories. Also takes care of making empty directories.

        :rtype: tuple of str
        :return: (unique event tag, path to output directory)
        """
        # Options for choosing how to name things. Legacy, manual or new-style
        if self._legacy_naming:
            logger.debug("reverting to legacy style file naming")
            event_tag = format_event_tag_legacy(self.event)
        else:
            event_tag = format_event_tag(self.event)

        if self._overwrite_event_tag:
            logger.debug("overwriting automatically generated event tag")
            event_tag = self._overwrite_event_tag

        logger.info(f"event tag is: {event_tag}")
        full_output_dir = os.path.join(self.output_dir, event_tag)
        if not os.path.exists(full_output_dir):
            os.makedirs(full_output_dir)
        elif not self._overwrite:
            logger.warning(f"output directory '{full_output_dir}' exists and "
                           f"overwrite flag (-o/--overwrite) not set, exiting")
            sys.exit(0)
        return event_tag, full_output_dir

    def run(self, event=None, inv=None, st=None, **kwargs):
        """
        Run PySEP: Seismogram Extraction and Processing. Steps in order are:

            1) Set default parameters or load from config file
            2) Check parameter validity, exit if unexpected values
            3) Get data and metadata (QuakeML, StationXML, waveforms)
            4) Remove unacceptable stations based on user-defined criteria
            5) Remove unacceptable waveforms based on user-defined criteria
            6) Generate some new metadata for tagging and output
            7) Pre-process waveforms and standardize for general use
            8) Generate output files and figures as end-product

        :type event: obspy.core.event.Event
        :param event: optional user-provided event object which will force a
            skip over QuakeML/event searching
        :type inv: obspy.core.inventory.Inventory
        :param inv: optional user-provided inventory object which will force a
            skip over StationXML/inventory searching
        :type st: obspy.core.stream.Stream
        :param st: optional user-provided strean object which will force a
            skip over waveform searching
        """
        # Overload default parameters with event input file and check validity
        self.load(**kwargs)
        self.check()
        self.c = self.get_client()

        # Get metadata (QuakeML, StationXML)
        if event is None:
            self.event = self.get_event()
        else:
            self.event = event

        self.event_tag, self.output_dir = self._event_tag_and_output_dir()

        if inv is None:
            self.inv = self.get_stations()
        else:
            self.inv = inv
        self.inv = self.curtail_stations()

        # Get waveforms, format and assess waveform quality
        if st is None:
            self.st = self.get_waveforms()
        else:
            self.st = st
        self.st = quality_check_waveforms_before_processing(self.st)
        self.st = append_sac_headers(self.st, self.event, self.inv)
        if self.taup_model is not None:
            self.st = format_sac_header_w_taup_traveltimes(self.st, 
                                                           self.taup_model)

        # Waveform preprocessing and standardization
        self.st = self.preprocess()
        if self.rotate is not None:
            self.st = self.rotate_streams()
        self.st = quality_check_waveforms_after_processing(self.st)

        # Generate outputs for user consumption
        self.write(**kwargs)
        self.plot()


def parse_args():
    """
    Define command line arguments. Allow user to set Pysep completely from
    the command line (although that would be extremely cumbersome, best to
    use Config files)
    """
    parser = argparse.ArgumentParser(
        description="PYSEP: Python Seismogram Extraction and Processing",
    )
    # Exposing some of the more useful parameters as public arguments
    parser.add_argument("-c", "--config", default="", type=str, nargs="?",
                        help="path to a YAML config file which defines "
                             "parameters used to control PySEP")
    parser.add_argument("-p", "--preset", default="", type=str, nargs="?",
                        help="Overwrites '-c/--config', use one of the default "
                             "in-repo config files which have been previously "
                             "designed. See 'pysep/pysep/configs' for options")
    parser.add_argument("-e", "--event", default="", type=str, nargs="?",
                        help="Required if using '-p/--preset'. Each preset "
                             "config may define >1 event. This flag determines"
                             "which event to run. If `event`=='all', will "
                             "gather ALL events in the preset.")
    parser.add_argument("-E", "--event_file", default="", type=str, nargs="?",
                        help="Allows using a single config file to gather"
                             "data for multiple events. Overwrites "
                             "'origin_time' parameter in the original config "
                             "file.")
    parser.add_argument("-U", "--user", default=None, type=str, nargs="?",
                        help="Username if required to access IRIS webservices")
    parser.add_argument("-P", "--password", default=None, type=str, nargs="?",
                        help="Password if required to access IRIS webservices")
    parser.add_argument("-W", "--write", default=False, action="store_true",
                        help="Write out a blank configuration file to be "
                             "filled in by the User")
    parser.add_argument("-l", "--list", default=False, action="store_true",
                        help="list out avaialable `preset` config options")
    parser.add_argument("-L", "--log_level", default="INFO", type=str,
                        nargs="?", help="verbosity of logging: 'WARNING', "
                                        "'INFO', 'DEBUG'")
    parser.add_argument("--legacy_naming", default=False, action="store_true",
                        help="use the file naming schema and directory "
                             "structure of the legacy version of PySEP.")
    parser.add_argument("--overwrite_event_tag", default="", type=str,
                        nargs="?",
                        help="Manually set the event tag used to name the "
                             "output directory and SAC files. If None, will "
                             "default to a tag consisting of event origin time "
                             "and region, or just origin time if using "
                             "`--legacy_naming`")
    parser.add_argument("-o", "--overwrite", default=False, action="store_true",
                        help="overwrite existing directory which matches "
                             "the unique event tag")

    # Keyword arguments can be passed directly to the argparser in the same
    # format as the above kwargs (e.g., --linewidth 2), but they will not have
    # help messages or type checking
    parsed, unknown = parser.parse_known_args()
    for arg in unknown:
        if arg.startswith(("-", "--")):
            parser.add_argument(arg.split("=")[0])

    return parser.parse_args()


def get_data(config_file=None, event=None, inv=None, st=None, write_files=None,
             plot_files=None, log_level=None, *args, **kwargs):
    """
    Interactive/scripting function to run PySep and return quality controlled,
    SAC-headed stream object which can then be used for other processes.

    .. note::
        By default turns file writing and plotting OFF so that this function
        acts solely as a data collection/processing call.

    .. note::
        args and kwargs are passed directly to Pysep.__init__() so you can
        define all your parameters in this call, or through a config file

    .. rubric::
        >>> from pysep import get_data
        >>> st = get_data(config_file='config.yaml')

    :type config_file: str
    :param config_file: path to YAML config file which will overload any default
        configs
    :type event: obspy.core.event.Event
    :param event: optional user-provided event object which will force a
        skip over QuakeML/event searching
    :type inv: obspy.core.inventory.Inventory
    :param inv: optional user-provided inventory object which will force a
        skip over StationXML/inventory searching
    :type st: obspy.core.stream.Stream
    :param st: optional user-provided strean object which will force a
        skip over waveform searching
    :type write_files: list or None
    :param write_files: list of files to write, acceptable options defined in
        write(). Defaults to None, no files will be written
    :type plot_files: list or None
    :param plot_files: list of files to plot, acceptable options defined in
        plot(). Defaults to None, no figures will be made
    :type log_level: str or None
    :param log_level: verbosity of logger. Defaults to no logging to mimic
        a standard function call rather than a standalone package
    :rtype: tuple of (obspy.core.event.Event, obspy.core.inventory.Inventory,
                      obspy.core.stream.Stream)
    :return: returns obspy objects defining data and metadata that have been
        collected by PySEP
    """
    sep = Pysep(config_file=config_file, write_files=write_files,
                plot_files=plot_files, log_level=log_level, *args, **kwargs)
    sep.run(event=event, inv=inv, st=st)
    return sep.event, sep.inv, sep.st


def _print_preset_configs():
    """
    Print out a list of Config files that can be used to run PySEP
    """
    pysep_dir = Path(__file__).parent.absolute()

    filelist = glob(os.path.join(pysep_dir, "configs", "*", "*yaml"))
    filelist = sorted([fid.split("/")[-2:] for fid in filelist])
    print(f"-p/--preset -e/--event")
    for i, fid in enumerate(filelist):
        preset, event = fid
        print(f"-p {preset} -e {event}")


def _return_matching_preset_configs(args):
    """
    Return a list of preset configs from the itnernal PySEP directory
    Searches through the internal pysep/configs/*/*.yaml files for matching
    user queries.

    :type args: argparser args
    :param args: command line argumments
    :rtype: list
    :return: matching config files that will be used to run PySEP
    """
    pysep_dir = Path(__file__).parent.absolute()

    preset_glob = os.path.join(pysep_dir, "configs", "*")
    # Ignore the template config files, only look at directories
    presets = [os.path.basename(_) for _ in glob(preset_glob)
               if not _.endswith(".yaml")]
    if not args.list:
        assert (args.preset in presets), (
            f"chosen preset (-p/--preset) {args.preset} not in available: "
            f"{presets}"
        )

    event_glob = os.path.join(pysep_dir, "configs", args.preset, "*")
    event_paths = sorted(glob(event_glob))
    event_names = [os.path.basename(_) for _ in event_paths]

    # 'event' arg can be an index (int) or a name (str)
    try:
        event_idx = int(args.event)
        config_files = [event_paths[event_idx]]
    except (ValueError, TypeError):
        if args.event.upper() == "ALL":
            config_files = event_paths
        else:
            assert (args.event in event_names), (
                f"chosen event (-e/--event) must be an integer "
                f"index < {len(event_names)} or one of: {event_names}"
            )
            config_files = [event_paths[event_names.index(args.event)]]

    return config_files


def main():
    """
    Command-line-tool function which parses command line arguments and runs
    data gathering and processing. This is accessed through the command line:

    .. rubric::
        $ pysep -c config.yaml
    """
    config_files, events = [], None
    args = parse_args()
    # Write out a blank configuration file to use as a template
    if args.write:
        sep = Pysep()
        sep.write_config()
        return
    # List out available configurations inside the repo
    elif args.list:
        _print_preset_configs()
        return
    # Choose one of the available preset configuration files to run
    elif args.preset:
        config_files = _return_matching_preset_configs(args)
    # If not preset, user should specify the config file themselves
    elif args.config:
        assert(os.path.exists(args.config)), (
            f"config file (-c/--config) {args.config} does not exist"
        )
        config_files = [args.config]
        if args.event_file:
            events = read_event_file(args.event_file)
    # Or allow user to manually set all the input parameters through CLI
    else:
        config_files = [None]

    if events is None:
        # Normal run: parse through config files and run PySEP
        for config_file in config_files:
            sep = Pysep(config_file=config_file, **vars(args))
            sep.run()
    else:
        # Event File run, use the same config by parse through events
        for event in events:
            sep = Pysep(config_file=config_files[0],
                        origin_time=event["origin_time"],
                        event_latitude=event["event_latitude"],
                        event_longitude=event["event_longitude"],
                        event_depth_km=event["event_depth_km"],
                        event_magnitude=event["event_magnitude"],
                        **vars(args))
            sep.run(overwrite_event=False)


if __name__ == "__main__":
    main()
