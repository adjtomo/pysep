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
from pysep.utils.fmt import format_event_tag, get_codes
from pysep.utils.io import read_yaml, write_stations_file
from pysep.utils.llnl import scale_llnl_waveform_amplitudes
from pysep.utils.process import (merge_and_trim_start_end_times, resample_data,
                                 format_streams_for_rotation, rotate_to_uvw)
from pysep.utils.plot import plot_source_receiver_map
from pysep.recsec import plotw_rs


class Pysep:
    """
    Download, preprocess, and save waveform data from IRIS using ObsPy
    """
    def __init__(self, config_file=None, event_selection="default",
                 client="IRIS", origin_time=None, reference_time=None,
                 networks="*", stations="*", locations="*", channels="*",
                 event_latitude=None, event_longitude=None, event_depth_km=None,
                 event_magnitude=None, remove_response=True,
                 remove_clipped=False, water_level=60, detrend=True,
                 demean=True, taper_percentage=0, rotate=None, pre_filt=None,
                 mindistance=0, maxdistance=20E3, minazimuth=0, maxazimuth=360,
                 minlatitude=None, minlongitude=None, maxlatitude=None,
                 maxlongitude=None, resample_freq=50, scale_factor=1,
                 phase_list=None,
                 seconds_before_event=20, seconds_after_event=20,
                 seconds_before_ref=100, seconds_after_ref=300,
                 taup_model="ak135", output_unit="VEL",  user=None,
                 password=None, client_debug=False, log_level="DEBUG",
                 timeout=600, write_files="all", plot_files="all",
                 llnl_db_path=None, output_dir=None, overwrite=False, **kwargs):
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
        self._overwrite = overwrite

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

        if self.client.upper() == "PH5":
            assert(self._user is not None and self._password is not None), (
                "`client`=='PH5' requires `-U/--user` and `-P/--password`"
            )

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

        if self.pre_filt is not None:
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

    def load(self, config_file=None):
        """
        Overwrite default parameters using a YAML config file
        """
        if config_file is None:
            config_file = self.config_file

        if config_file is not None:
            logger.info(f"overwriting default parameters with config file: "
                        f"'{config_file}'")
            config = read_yaml(config_file)
            for key, val in config.items():
                if hasattr(self, key):
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

        :rtype: obspy.core.stream.Stream
        :return: Stream of channel-separated waveforms
        """
        logger.info(f"querying client '{self.client.upper()}' for waveforms")

        minimumlength = self.seconds_before_ref + self.seconds_after_ref
        logger.debug(f"asserting that the minimum waveform length should be"
                     f"{minimumlength:.2f}")

        # PH5 should not be queried in bulk as it is a small request
        if self.client.upper() == "PH5":
            st = self.c.get_waveforms(
                network=self.networks, location=self.locations,
                station=self.stations, channel=self.channels,
                starttime=self.origin_time - self.seconds_before_ref,
                endtime=self.origin_time - self.seconds_after_ref,
                minimumlength=minimumlength
            )
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
            st = self.c.get_waveforms_bulk(bulk=bulk, minimumlength=t2-t1)
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
            One function function currently, but we can expand curtailing here
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
            for code in get_codes(st=self.st, choice="channel", suffix="?",
                                  up_to=True):
                net, sta, loc, cha = code.split(".")
                st_sta = st_out.select(network=net, station=sta, location=loc,
                                       channel=cha)
                try:
                    st_sta.remove_response(inventory=self.inv,
                                           water_level=self.water_level,
                                           pre_filt=self.pre_filt,
                                           taper=bool(self.taper_percentage),
                                           taper_fraction=self.taper_percentage,
                                           zero_mean=self.demean,
                                           output=self.output_unit)
                except ValueError as e:
                    for tr in st_sta:
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

    def write(self, write_files=None, _return_filenames=False):
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
                write_cap_weights_files(
                    st=self.st, event=self.event, order_by=order_by,
                    path_out=self.output_dir
                )

        if "config_file" in write_files or "all" in write_files:
            logger.info("writing config YAML file")
            fid = os.path.join(self.output_dir, f"pysep_config.yaml")
            self.write_config(fid=fid)

        if "station_list" in write_files or "all" in write_files:
            fid = os.path.join(self.output_dir, "stations_list.txt")
            logger.info("writing stations file")
            logger.debug(fid)
            write_stations_file(self.inv, self.event, fid)

        if "inv" in write_files or "all" in write_files:
            fid = os.path.join(self.output_dir, f"inv.xml")
            logger.info("writing inventory as StationXML")
            logger.debug(fid)
            self.inv.write(fid, format="STATIONXML")

        if "event" in write_files or "all" in write_files:
            fid = os.path.join(self.output_dir, f"event.xml")
            logger.info("writing event as QuakeML")
            logger.debug(fid)
            self.event.write(fid, format="QuakeML")

        if "stream" in write_files or "all" in write_files:
            fid = os.path.join(self.output_dir, f"stream.ms")
            logger.info("writing waveform stream in MiniSEED")
            logger.debug(fid)
            with warnings.catch_warnings():
                # ignore the encoding warning that comes from ObsPy
                warnings.simplefilter("ignore")
                self.st.write(fid, format="MSEED")

        if "sac" in write_files or "all" in write_files:
            logger.info("writing each waveform trace in SAC format")
            _output_dir = os.path.join(self.output_dir, "SAC")
            if not os.path.exists(_output_dir):
                os.makedirs(_output_dir)
            for tr in self.st:
                # e.g., 2000-01-01T000000.NN.SSS.LL.CCC.sac
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
            fid = os.path.join(self.output_dir, f"station_map.png")
            plot_source_receiver_map(self.inv, self.event, fid)

        if "record_section" in self.plot_files or "all" in self.plot_files:
            fid = os.path.join(self.output_dir, f"record_section.png")
            # Default settings to create a general record section
            plotw_rs(st=self.st, sort_by="distance_r",
                     scale_by="normalize", overwrite=True, save=fid)

    def _event_tag_and_output_dir(self):
        """
        Convenience function to establish and naming schema for files and
        directories. Also takes care of making empty directories.

        :rtype: tuple of str
        :return: (unique event tag, path to output directory)
        """
        event_tag = format_event_tag(self.event)
        logger.info(f"event tag is: {self.event_tag}")
        full_output_dir = os.path.join(self.output_dir, event_tag)
        if not os.path.exists(full_output_dir):
            os.makedirs(full_output_dir)
        elif not self._overwrite:
            logger.warning(f"output directory '{full_output_dir}' exists and "
                           f"overwrite flag (-o/--overwrite) not set, exiting")
            sys.exit(0)
        return event_tag, full_output_dir

    def run(self, event=None, inv=None, st=None):
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
        self.load()
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
        self.st = format_sac_header_w_taup_traveltimes(self.st, self.taup_model)

        # Waveform preprocessing and standardization
        self.st = self.preprocess()
        self.st = self.rotate_streams()
        self.st = quality_check_waveforms_after_processing(self.st)

        # Generate outputs for user consumption
        self.write()
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
    parser.add_argument("-U", "--user", default=None, type=str, nargs="?",
                        help="Username if required to access IRIS webservices")
    parser.add_argument("-P", "--password", default=None, type=str, nargs="?",
                        help="Password if required to access IRIS webservices")
    parser.add_argument("-l", "--list", default=False, action="store_true",
                        help="list out avaialable `preset` config options")
    parser.add_argument("-L", "--log_level", default="INFO", type=str,
                        nargs="?", help="verbosity of logging: 'WARNING', "
                                        "'INFO', 'DEBUG'")
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


def main():
    """
    Command-line-tool function which parses command line arguments and runs
    data gathering and processing. This is accessed through the command line:

    .. rubric::
        $ pysep -c config.yaml
    """
    args = parse_args()
    # Allow grabbing preset Config files from inside the repo
    if args.preset or args.list:
        pysep_dir = Path(__file__).parent.absolute()
        # If '-l/--list', just print out available config objects
        if args.list:
            filelist = glob(os.path.join(pysep_dir, "configs", "*", "*"))
            filelist = sorted([fid.split("/")[-2:] for fid in filelist])
            print(f"-p/--preset -e/--event")
            for i, fid in enumerate(filelist):
                preset, event = fid
                print(f"-p {preset} -e {event}")
            return

        preset_glob = os.path.join(pysep_dir, "configs", "*")
        # Ignore the template config files, only look at directories
        presets = [os.path.basename(_) for _ in glob(preset_glob)
                   if not _.endswith(".yaml")]
        if not args.list:
            assert(args.preset in presets), (
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
                assert(args.event in event_names), (
                    f"chosen event (-e/--event) must be an integer "
                    f"index < {len(event_names)} or one of: {event_names}"
                )
                config_files = [event_paths[event_names.index(args.event)]]
    # If not preset, user should specify the config file
    elif args.config:
        assert(os.path.exists(args.config)), (
            f"config file (-c/--config) {args.config} does not exist"
        )
        config_files = [args.config]
    # Finally, allow user to manually set all the input parameters through
    # the command line using '--' flags. Not recommended
    else:
        config_files = [None]

    for config_file in config_files:
        sep = Pysep(config_file=config_file, **vars(args))
        sep.run()


if __name__ == "__main__":
    main()
