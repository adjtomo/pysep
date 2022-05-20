"""
Python Seismogram Extraction and Processing

Download, pre-process, and organize seismic waveforms, event and station
metadata
"""
import os
import sys
import numpy as np
import yaml
from obspy import UTCDateTime, Stream
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNBadRequestException
from obspy.core.event import Event, Origin, Magnitude

from pysep import logger
from pysep.utils.io import read_yaml
from pysep.utils.cap_sac import (append_sac_headers, write_cap_weights_files,
                                 format_sac_headers_post_rotation)
from pysep.utils.curtail import (curtail_by_station_distance_azimuth,
                                 quality_check_waveforms)
from pysep.utils.fetch import fetch_bulk_station_list_taup
from pysep.utils.llnl import scale_llnl_waveform_amplitudes
from pysep.utils.process import (zerophase_chebychev_lowpass_filter, get_codes,
                                 format_streams_for_rotation, rotate_to_UVW)
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
                 min_dist=0, max_dist=20E3, min_az=0, max_az=360,
                 min_lat=None, min_lon=None, max_lat=None, max_lon=None,
                 resample_freq=50, scale_factor=1, phase_list=None,
                 seconds_before_event=20, seconds_after_event=20,
                 seconds_before_ref=100, seconds_after_ref=300,
                 taupmodel="ak135", output_unit="VEL",  user=None,
                 password=None, client_debug=False, log_level="DEBUG",
                 timeout=600, write_files=None, plot_files=None,
                 llnl_db_path=None, output_dir=None, **kwargs):
        """
        Define a default set of parameters

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
        :type filter_type: str
        :param filter_type: type of ObsPy filter to use: 'bandpass', 'highpass',
            or 'lowpass'
        :type f1: float
        :param f1: minimum frequency in hz
        :type f2: float
        :param f2: maximum frequency in hz
        :type zerophase: bool
        :param zerophase: apply a two-pass filter to keep time shift
            information, defaults True
        :type corners: int
        :param corners: filter corners to apply
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
        :type rotateRTZ: bool
        :param rotateRTZ: rotate components to Radial, Transverse, Vertical (Z)
        :type rotateUVW: bool
        :param rotateUVW: rotate into a UVW coordinate system
        :type rotateZNE: bool
        :param rotateZNE: rotate to East, North, Vertical (Z)
        :type min_dist: float
        :param min_dist: get waveforms from stations starting from a minimum
            distance away  from event (km)
        :type max_dist: float
        :param max_dist: get waveforms from stations up to a maximum distance
            from event (km)
        :type min_az: float
        :param min_az: get waveforms from a station starting at a given
            minimum azimuth (deg) out to a final maximum azimuth (deg) where
            0 degrees points North from the event
        :type max_az: float
        :param max_az: get waveforms from a station starting at a given
            minimum azimuth (deg) out to a final maximum azimuth (deg) where
            0 degrees points North from the event
        :type min_lat: float
        :param min_lat: minimum latitude for bounding box to search for stations
        :type max_lat: float
        :param max_lat: maximum latitude for bounding box to search for stations
        :type min_lon: float
        :param min_lon: minimum longitude for bounding box to search for stas
        :type max_lon: float
        :param max_lon: maximum longitude for bounding box to search for stas
        :type resample_TF: bool
        :param resample_TF if False then resample_freq is taken from SAC files
        :type resample_freq: float
        :param resample_freq: frequency to resample data at
        :type scale_factor: float
        :param scale_factor: scale data by a constant, for CAP use 10**2
            (to convert m/s to cm/s)
        :type overwrite_ddir: bool
        :param overwrite_ddir: if True, delete any existing data directories
        :type phase_window: bool
        :param phase_window: trim waveforms between two phases determined by
            TauP
        :type phases: list of str
        :param phases: phases to write to SAC files or grab data from, e.g.,
            ['P', 'PP']
        :type sec_before_after_event: float
        :param sec_before_after_event: seconds before and after an event
            origin time to search for the event, allowing some wiggle room if
            the given origin_time is not exactly the catalog time
            will search O - T  <= t <= O + T
        :type tbefore_sec: float
        :param tbefore_sec: time before origintime to fetch waveform data
        :type tafter_sec: float
        :param tafter_sec: time after origintime to fetch waveform data
        :type use_catalog: bool
        :param use_catalog: do not gather event catalog from IRIS,
            but provide your own event information manually
        :type ifph5: bool
        :param ifph5: Use PH5 data format from PASSCAL (Denali nodal data
            specific)
        :type write_sac_phase: bool
        :param write_sac_phase: include phase information in output SAC files
        :type taupmodel: str
        :param taupmodel: name of TauP model to use to calculate phase arrivals
        :type output_cap_weight_file: bool
        :param output_cap_weight_file: output a 'weights.dat' file formatted for
            use in the Cut-And-Paste moment tensor inversion code
        :type output_event_info: bool
        :param output_event_info: output an event information file
        :type output_unit: str
        :param output_unit: the output format of the waveforms, something like
            'DISP', 'VEL', 'ACC'
        :type ipre_filt: int
        :param ipre_filt: apply a prefilter to incoming data
            0 No pre_filter
            1 default pre_filter (see getwaveform_iris.py)
            2 user-defined pre_filter (use this if using bandpass filter)
            pre-filter for deconvolution
            https://ds.iris.edu/files/sac-manual/commands/transfer.html
            Pre-filter will not be applied if remove_response == False
        :type iplot_response: bool
        :param iplot_response: plot instrument response
        :type icreateNull: int
        :param icreateNull: create Null traces to allow for trace rotation
            with missing components, useful for missing vertical components
        :type isave_raw: bool
        :param isave_raw: save raw waveform data (before preprocessing
        :type ifFilter: bool
        :param ifFilter: flag to turn filtering on or off
        :type ifsave_sacpaz: Save SAC poles-and-zeros information, which is
            required as input for the MouseTrap module
        :type ifplot_spectrogram: bool
        :param ifplot_spectrogram: plot spectrograms
        :type ifsave_stationxml: bool
        :param ifsave_stationxml: save StationXML files, for adjoint tomography
            purposes usually
        :type ifverbose: bool
        :param ifverbose: print status during waveform getting
        :type ifsave_asdf: bool
        :param ifsave_asdf: save waveforms in ASDF format
        :type ifmass_downloader: bool
        :param ifmass_downloader: use ObsPy's mass downloader
        :type user: str
        :param user: if IRIS embargoes data behind passwords, user ID
        :type password: str
        :param password: if IRIS embargoes data behind passwords, password
        """
        self.config_file = config_file

        # ObsPy client-related parameters
        self.client = client
        self.client_debug = bool(client_debug)
        self.timeout = timeout
        self._user = user
        self._password = password

        # Parameters related to event selection
        self.event_selection = event_selection

        try:
            self.origin_time = UTCDateTime(origin_time)
        except TypeError:
            self.origin_time = None
        self.seconds_before_event = seconds_before_event
        self.seconds_after_event = seconds_after_event

        # Optional: if User wants to define an event on their own
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
        # Reference time allows throwing a static time shift from the origin
        # time, e.g., if there are timing errors with relation to o time
        self.reference_time = reference_time or self.origin_time
        self.seconds_before_ref = seconds_before_ref
        self.seconds_after_ref = seconds_after_ref
        self.phase_list = phase_list

        # This default is a UAF LUNGS system-specific database path
        self.llnl_db_path = llnl_db_path or \
                            "/store/raw/LLNL/UCRL-MI-222502/westernus.wfdisc"

        # Station curtailing criteria
        self.min_dist = min_dist
        self.max_dist = max_dist
        self.min_az = min_az
        self.max_az = max_az
        self.min_lat = min_lat
        self.max_lat = max_lat
        self.min_lon = min_lon
        self.max_lon = max_lon

        # Preprocessing flags
        self.demean = bool(demean)
        self.detrend = bool(detrend)
        self.taper_percentage = taper_percentage
        self.rotate = rotate
        self.remove_response = bool(remove_response)
        self.output_unit = output_unit
        self.water_level = water_level
        self.pre_filt = pre_filt
        self.scale_factor = scale_factor
        self.resample_freq = resample_freq
        self.remove_clipped = bool(remove_clipped)

        # Program related parameters
        self.output_dir = output_dir or os.getcwd()
        self.write_files = write_files or {"all"}
        self.plot_files = plot_files or {"all"}
        self.log_level = log_level

        # Internally filled attributes
        self._client = None
        self.event_tag = None
        self._azimuths = {}
        self._distances = {}


    def check(self):
        """
        Check input parameter validity against expected Pysep behavior
        """
        self.origin_time = UTCDateTime(self.origin_time)
        if self.reference_time is None:
            self.reference_time = self.origin_time
        self.reference_time = UTCDateTime(self.reference_time)

        if self.event_selection == "catalog" or self.client == "LLNL":
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
            raise ValueError("`event_selection` must seconds_before_evente in "
                             "'catalog' or 'default'")

        if self.client.upper() == "PH5":
            assert(self._user is not None and self._password is not None), (
                "`client`=='PH5' requires `user` and `password`"
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

        assert(0 <= self.min_az <= 360), f"0 <= `min_az` <= 360"
        assert(0 <= self.max_az <= 360), f"0 <= `max_az` <= 360"

        if self.rotate is not None:
            acceptable_rotations = {"RTZ", "UVW", "ZNE"}
            self.rotate = [val.upper() for val in self.rotate]
            assert(set(self.rotate).issubset(acceptable_rotations)), (
                f"`rotate` must be a subset of: {acceptable_rotations}"
            )

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
        acceptable_write_files = {"weights", "weights_az",
                                  "weights_dist", "weights_name",
                                  "config", "all"}
        assert(self.write_files.issubset(acceptable_write_files)), (
            f"`write_files` must be a list of some or all of: "
            f"{acceptable_write_files}"
        )
        self.write_files = set(self.write_files)

        acceptable_plot_files = {"map", "record_section", "all"}
        assert(self.plot_files.issubset(acceptable_plot_files)), (
            f"`plot_files` must be a list of some or all of: "
            f"{acceptable_plot_files}"
        )
        self.plot_files = set(self.plot_files)


    def get_client(self):
        """
        Options to choose different Clients based on client name. These clients
        are ONLY for waveform gathering. We assume that
        """
        if self.client is None:
            return None

        # Lawrence Livermore Natinoal Lab internal waveform database
        if self.client.upper() == "LLNL":
            try:
                import llnl_db_client
                c = llnl_db_client.LLNLDBClient(self.llnl_db_path)
            except ImportError as e:
                logger.warning("Cannot import `llnl_db_client`, please "
                               "download and install from: "
                               "https://github.com/krischer/llnl_db_client")
                raise ImportError from e
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

    def load_config(self, config_file=None):
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
                setattr(self, key, val)

    def get_event(self):
        """
        Download event metadata from IRIS using ObsPy functionality
        """
        if self.event_selection == "catalog":
            event = self._get_event_catalog()
        else:
            if self.client.upper() == "LLNL":
                event = self._get_event_llnl()
            else:
                event = self._get_event_default()

        return event

    def _get_event_catalog(self):
        """
        TODO Assert origin time and second before/after event set
        :return:
        """
        logger.info(f"getting event information with client {self.client}")

        cat = self._client.get_events(
            starttime=self.origin_time - self.seconds_before_event,
            endtime=self.origin_time + self.seconds_after_event,
            debug=self.client_debug
        )
        logger.info(f"{len(cat)} matching events found for "
                    f"{self.seconds_before_event}s < {self.origin_time} "
                    f"< {self.seconds_after_event}, choosing first")

        event = cat[0]

        return event

    def _get_event_default(self):
        """
        Make a barebones event object based on user-defined parameters which
        will then be used to query for waveforms and StationXML data
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

    def _get_event_llnl(self):
        """
        Special getter function for Lawrence Livermore National Lab data
        LLNL database has a special client
        :return:
        """
        logger.info("getting event information from LLNL database")

        cat = self._client.get_catalog()
        mintime_str = f"time > {self.origin_time - self.seconds_before_event}"
        maxtime_str = f"time < {self.origin_time + self.seconds_after_event}"

        cat = cat.filter(mintime_str, maxtime_str)
        logger.debug(f"{len(cat)} events in LLNL catalog, picking zeroth index")
        event = cat[0]

        return event

    def get_waveforms(self):
        """
        Download data from IRIS using ObsPy functionality

        TODO add check for no '-' in station for NCEDC stations
        TODO add warning for "*" may take long
        :return:
        """
        if self.client.upper() == "PH5":
            st = self._get_waveforms_ph5()
        else:
            # Allow setting start and endtimes based on seismic phases
            if self.phase_list:
                bulk = fetch_bulk_station_list_taup()
            # Make the bulk selection list manually
            else:
                bulk = self._make_bulk_station_list()
            try:
                st = self._client.get_waveforms_bulk(bulk=bulk)
            except FDSNBadRequestException:
                logger.warning(f"client {self.client} returned no waveforms, "
                               f"please check your event and station "
                               f"parameters again and re-submit")
                sys.exit(-1)

        logger.info(f"{len(st)} waveforms returned after query")

        return st

    def _get_waveforms_ph5(self):
        """
        Download waveforms from PH5 dataservice
        :return:
        """
        logger.debug("")
        st = self._client.get_waveforms(
            network=self.networks,  location=self.locations,
            station=self.stations, channel=self.channels,
            starttime=self.origin_time - self.seconds_before_ref,
            endtime=self.origin_time - self.seconds_after_ref
        )
        return st

    def _make_bulk_station_list(self):
        """
        ObsPy Client.get_waveforms_bulk() requires a list input
        """
        bulk = []
        t1 = self.reference_time - self.seconds_before_ref
        t2 = self.reference_time + self.seconds_after_ref
        for net in self.inv:
            for sta in net:
                # net sta loc cha t1 t2
                bulk.append((net.code, sta.code, "*", self.channels, t1, t2))
        return bulk


    def get_stations(self):
        """
        Download station metadata from IRIS using ObsPy functionality
        :return:
        """
        inv = self._client.get_stations(
            network=self.networks, location=self.locations,
            station=self.stations, channel=self.channels,
            starttime=self.origin_time - self.seconds_before_ref,
            endtime=self.origin_time + self.seconds_after_ref,
            minlatitude=self.min_lat, maxlatitude=self.max_lat,
            minlongitude=self.min_lon, maxlongitude=self.max_lon,
            level="response"
        )

        nnet = len(inv)
        ncha = len(inv.get_contents()["channels"])
        nsta = len(inv.get_contents()["stations"])

        logger.info(f"collected {ncha} channels from {nsta} stations in "
                    f"{nnet} networks from {self.client}")

        return inv

    def curtail_stations(self):
        """
        Remove stations based on station distance, azimuth, etc.
        :return:
        """
        inv = self.inv.copy()
        inv, distances, azimuths = curtail_by_station_distance_azimuth(
            event=self.event, inv=inv, min_dist=self.min_dist,
            max_dist=self.max_dist, min_az=self.min_az, max_az=self.max_az
        )
        return inv, distances, azimuths

    def preprocess(self):
        """
        Very simple preprocessing to remove response and apply a prefilter

        TODO do we need the pre-filt if we are doing it in resample?
        """
        st_out = self.st.copy()
        if self.detrend:
            logger.info(f"applying linear detrend to all data")
            st_out.detrend("linear")
        if self.remove_response:
            logger.info(f"removing response, output units in: "
                        f"{self.output_unit}")
            st_out.remove_response(inv=self.inv, water_level=self.water_level,
                                   pre_filt=self.pre_filt,
                                   taper=bool(self.taper_percentage),
                                   taper_fraction=self.taper_percentage,
                                   zero_mean=self.demean,
                                   output=self.output_unit
                                   )
        if self.scale_factor:
            logger.info(f"applying amplitude scale factor {self.scale_factor}")
            for tr in st_out:
                tr.data = tr.data * self.scale_factor
                tr.stats.sac["scale"] = self.scale_factor
        if self.client == "LLNL":
            # This won't do anything if we don't have any 'LL' network codes
            st_out = scale_llnl_waveform_amplitudes(st_out)

        return st_out

    def _prep_resample(self):
        """
        Resample data to desired sampling rate and throw in a decimation filter

        TODO ignoring `resample_cut` which doesn't seem to do anything but
            is included in old code
        :return:
        """
        logger.info(f"resampling data to sampling rate: {self.resample_freq}")
        st_out = self.st.copy()
        for i, tr in enumerate(st_out[:]):
            target_nyquist = 0.5 * self.resample_freq
            current_nyquist = 0.5 * tr.stats.sampling_rate
            if target_nyquist < current_nyquist:
                try:
                    # This function affects the trace in-place
                    zerophase_chebychev_lowpass_filter(
                        tr=tr, freqmax=target_nyquist
                    )
                except Exception as e:
                    logger.warning(f"exception in lowpass filtering "
                                   f"{tr.get_id()}, remove")
                    logger.debug(e)
                    st_out.remove(tr)
            else:
                tr.detrend("linear")
            tr.taper(max_percentage=0.01, type="hann")
            # Enforce that this data array is 'c'ontiguous
            tr.data = np.require(tr.data, requirements=["C"])

        return st_out

    def resample_and_trim_start_end_times(self):
        """
        Trim the maximum start and minumum end times for each station to get
        all waveforms to start and end at the same time
        :return:
        """
        logger.info("trimming start and end times on a per-station basis")
        st_edit = self._prep_resample()
        st_out = Stream()
        # e.g., NN.SSS.LL.CC?  i.e., only getting per-station codes
        codes = get_codes(st=st_edit, choice="channel", suffix="?")

        for code in codes:
            # Subset the stream based on the N number of components
            net, sta, loc, cha = code.split(".")
            st_edit_select = st_edit.select(network=net, station=sta,
                                            location=loc, channel=cha)

            max_start = max([tr.stats.starttime for tr in st_edit_select])
            min_end = min([tr.stats.endtime for tr in st_edit_select])
            npts = int((min_end - max_start) * self.resample_freq)

            # Lanczos is the preferred ObsPy interpolation method. NOT the same
            # used by SAC (weighted_average_slopes)
            # `a` is passed to lanczos_interpolation and defines the width of
            #   the window in samples on either side. Larger values of `a` are
            #   more expensive but better for interpolation
            st_edit_select.interpolate(sampling_rate=self.resample_freq,
                                       method="lanczos", starttime=max_start,
                                       npts=npts, a=8)

            st_out += st_edit_select
        return st_out

    def rotate_streams(self):
        """
        Rotate arbitrary three-component seismograms to desired orientation

        TODO Original code had a really complicated method for rotating,
            is that necesssary? Why is st.rotate() not acceptable
        """
        st_out = self.st.copy()

        st_out = format_streams_for_rotation(st_out)
        if "ZNE" in self.rotate:
            logger.info("rotating to components ZNE")
            st_out.rotate(method="->ZNE", inventory=self.inv)
        if "UVW" in self.rotate:
            logger.info("rotating to components UVW")
            st_out = rotate_to_UVW(st_out)
        if "RTZ" in self.rotate:
            logger.info("rotating to components RTZ")
            st_out.rotate("NE->RT", inventory=self.inv)

        st_out = format_sac_headers_post_rotation(st_out)

        return st_out

    def write(self, fmt, fid):
        """
        Write out various files specifying information about the collected
        stations and waveforms

        :param fmt:
        :return:
        """
        for weights_fid in ["weights_dist", "weights_az", "weights_name"]:
            if weights_fid in self.write_files or "all" in self.write_files:
                order_by = weights_fid.split("_")[1]
                write_cap_weights_files(
                    st=self.st, event=self.event, order_by=order_by,
                    path_out=os.path.join(self.output_dir, self.event_tag)
                )

        if "config_file" in self.write_files or "all" in self.write_files:
            self._write_config_file()

        if "station_list" in self.write_files or "all" in self.write_files:
            self._write_station_list()

        if "inv" in self.write_files or "all" in self.write_files:
            fid = os.path.join(self.output_dir, f"inv_{self.event_tag}.xml")
            self.inv.write(fid, format="STATIONXML")

        if "event" in self.write_files or "all" in self.write_files:
            fid = os.path.join(self.output_dir, f"event_{self.event_tag}.xml")
            self.event.write(fid, format="QuakeML")

        if "stream" in self.write_files or "all" in self.write_files:
            fid = os.path.join(self.output_dir, f"stream_{self.event_tag}.ms")
            self.st.write(fid, format="MSEED")

        if "sac" in self.write_files or "all" in self.write_files:
            for tr in self.st:
                # e.g., 2000-01-01T000000.NN.SSS.LL.CCC.sac
                tag = f"{self.event_tag}.{tr.get_id()}.sac"
                fid = os.path.join(self.output_dir, tag)
                tr.write(fid, format="SAC")


    def _write_config_file(self):
        """
        Write a YAML config file based on the parameters
        :return:
        """
        event_tag = self.event_tag
        if event_tag is None:
            tag = "pysep"
        fid = os.path.join(self.output_dir, f"config_{tag}.yaml")
        dict_out = vars(self)

        # Drop hidden variables
        dict_out = {key: val for key, val in dict_out.items()
                                                     if not key.startswith("_")}
        with open(fid, "w") as f:
            yaml.dump(dict_out, f, default_flow_style=False, sort_keys=False)


    def _write_station_list(self, fid=None):
        """
        Write a list of station codes, distances, etc.
        """
        if fid is None:
            fid = os.path.join(self.event_tag, f"station_list.txt")

        with open(fid, "w") as f:
            for net in self.inv:
                for sta in net:
                    netsta_code = f"{net.code}.{sta.code}"
                    f.write(f"{sta.code} {net.code} {sta.latitude} "
                            f"{sta.longitude} {self._distances[netsta_code]}"
                            f"{self._azimuths[netsta_code]}\n")

    def plot(self):
        """
        Plot map and record section if requested
        """
        if "map" in self.plot_files or "all" in self.plot_files:
            fid = os.path.join(self.output_dir, f"map_{self.event_tag}.png")
            plot_source_receiver_map(self.inv, self.event, fid)

        if "record_section" in self.plot_files or "all" in self.plot_files:
            fid = os.path.join(self.output_dir,
                               f"record_section_{self.event_tag}.png")
            plotw_rs(st=self.st, sort_by="distance_r", save=fid)

    def main(self):
        """
        Run Pysep seismogram extraction
        :return:
        """
        # Overload default parameters with event input file and check validity
        self.load_config()
        self.check()
        self._client = self.get_client()

        # Get waveforms and metadata (QuakeML, StationXML)
        self.event = self.get_event()
        self.inv = self.get_stations()
        self.inv, self._azimuths, self._distances = self.curtail_stations()
        self.st = self.get_waveforms()

        # Format and assess waveform quality
        self.st = append_sac_headers(self.st, self.inv, self.event)
        self.event_tag = self.st[0].stats.sac["kevnm"]  # event time
        self.output_dir = os.path.join(self.output_dir, self.event_tag)
        self.st = quality_check_waveforms(st=self.st)

        # Preprocessing and standardization
        self.st = self.preprocess()
        self.st = self.resample_and_trim_start_end_times()
        self.st = self.rotate_streams()

        self.write()
        self.plot()


# if __name__ == "__main__":
#     # !!! Parse arguments
#     sep = Pysep(config_file=args.config_file)
#     sep.main()
