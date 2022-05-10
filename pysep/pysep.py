"""
Python Seismogram Extraction and Processing

Download, pre-process, and organize seismic waveforms, event and station
metadata
"""
import os
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.core.event import Event, Origin, Magnitude, Catalog

from pysep import logger
from pysep.utils.read import load_yaml

class Pysep:
    """
    Download, preprocess, and save waveform data from IRIS using ObsPy
    """
    def __init__(self, config_file=None,  client_name="IRIS", origin_time=None,
                 network="*", station="*", location="*", channel="*",
                 remove_response=True, remove_clipped=False,
                 filter_type="bandpass", f1=1/40, f2=1/10, zerophase=True,
                 corners=4, water_level=60, detrend=True, demean=True,
                 taper=False, rotateRTZ=True, rotateUVW=False, rotateENZ=True,
                 min_dist=0, max_dist=20E3, min_az=0, max_az=360,
                 min_lat=None, min_lon=None, max_lat=None, max_lon=None,
                 resample_TF=False, resample_freq=50, scale_factor=1,
                 overwrite_ddir=True, phase_window=False, phases=None,
                 seconds_before_event=20, seconds_after_event=20,
                 tbefore_sec=100, tafter_sec=300, use_catalog=False, ifph5=False,
                 write_sac_phase=False, taupmodel="ak135",
                 output_cap_weight_file=True, output_event_info=True,
                 outformat="VEL",  ipre_filt=1, iplot_response=False,
                 icreateNull=False, isave_raw=False, isave_raw_processed=True,
                 isave_rotateENZ=True, isave_ENZ=True,
                 ifFilter=False, ifsave_sacpaz=False, ifplot_spectrogram=False,
                 ifsave_stationxml=True, ifverbose=True, ifsave_asdf=False,
                 ifmass_downloader=False, user=None, password=None, **kwargs):
        """
        Define a default set of parameters

        :type client_name: str
        :param client_name: ObsPy FDSN client to query data from, e.g., IRIS, LLNL,
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
        :type rotateENZ: bool
        :param rotateENZ: rotate to East, North, Vertical (Z)
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
            the given origintime is not exactly the catalog time
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
        :type outformat: str
        :param outformat: the output format of the waveforms, something like
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
        self.client_name = client_name
        self.client = self.get_client(user=user, password=password)


        # Catalog event selection: Parameters for gathering from Client
        self.event_selection = "default"  # catalog, default
        self.seconds_before_event = seconds_before_event
        self.seconds_after_event = seconds_after_event

        # Default event selection: User-provided event parameters
        self.origin_time = UTCDateTime(origin_time)
        self.event_latitude = None
        self.event_longitude = None
        self.event_depth_km = None
        self.event_magnitude = None

        # Waveform selection criteria
        self.reference_time = reference_time
        self.seconds_before_ref = tbefore_sec
        self.seconds_after_ref = tafter_sec

        # Station-related parameters
        self.network = network
        self.station = station
        self.station = station
        self.channel = channel
        self.location = location
        self.min_dist = min_dist
        self.max_dist = max_dist
        self.min_az = min_az
        self.max_az = max_az
        self.min_lat = min_lat
        self.max_lat = max_lat
        self.min_lon = min_lon
        self.max_lon = max_lon
        self.overwrite_ddir = overwrite_ddir
        self.icreateNull = icreateNull
        self.phases = phases
        self.write_sac_phase = write_sac_phase
        self.taupmodel= taupmodel

        # Filter-related parameters
        self.ifFilter = ifFilter
        self.filter_type = filter_type
        self.f1 = f1
        self.f2 = f2
        self.zerophase = zerophase
        self.corners = corners
        self.phase_window = phase_window

        # Pre-filter-related parameters
        self.ipre_filt = ipre_filt
        self.water_level = water_level
        self.f0 = 0.5 * self.f1
        self.f3 = 2.0 * self.f2
        self.pre_filt = (self.f0, self.f1, self.f2, self.f3)

        # For CAP (Cut-And-Paste moment tensor inversion code)
        self.resample_TF = resample_TF
        self.resample_freq = resample_freq
        self.scale_factor = scale_factor

        # Pre-processing (mainly for CAP)
        self.detrend = detrend
        self.demean = demean
        self.taper = taper
        self.outformat = outformat
        self.remove_response = remove_response

        # Options for rotation
        self.rotateRTZ = rotateRTZ
        self.rotateUVW = rotateUVW
        self.rotateENZ = rotateENZ

        # Output flags for how to save the raw or processed data
        self.ph5 = ifph5
        self.output_event_info = output_event_info
        self.output_cap_weight_file = output_cap_weight_file
        self.save_sacpaz = ifsave_sacpaz
        self.plot_spectrogram = ifplot_spectrogram
        self.plot_response = iplot_response
        self.save_stationxml = ifsave_stationxml
        self.save_asdf = ifsave_asdf
        self.save_raw = isave_raw
        self.save_raw_processed = isave_raw_processed
        self.save_ENZ = isave_ENZ
        self.save_rotateENZ = isave_rotateENZ

        self.remove_clipped = remove_clipped
        if self.remove_clipped:
            print("Removing clipped stations forces `isave_raw`==True")
            self.isave_raw = True

        # Flags to deal with miscellaneous items
        self.ifmass_downloader = ifmass_downloader


        # Event information to be filled
        self.ev = Event()
        self.ref_time_place = Event()

    def get_client(self, user=None, password=None):
        """
        Options to choose different Clients based on client name

        TODO add log level to parameters and use to set client debug
        :return:
        """
        if self.client_name is not None:
            if self.client_name.upper() == "LLNL":
                try:
                    import llnl_db_client
                except ImportError:
                    logger.warning("Cannot import llnl_db_client, please "
                                   "download from "
                                   "https://github.com/krischer/llnl_db_client")
                client = llnl_db_client.LLNLDBClient(
                    "/store/raw/LLNL/UCRL-MI-222502/westernus.wfdisc"
                )
            # IRIS DMC PH5WS Station Web Service
            elif self.client_name.upper() == "PH5":
                client = Client(
                    "http://service.iris.edu", debug=True,
                    service_mappings={
                        "station": "http://service.iris.edu/ph5ws/station/1"
                    }
                )
            else:
                client = Client(self.client_name, user=user, password=password,
                                debug=False, timeout=600)
        else:
            client = None

        return client

    def overwrite_parameters(self):
        """
        Overwrite default parameters using a config file
        :return:
        """
        if self.config_file is not None:
            logger.info(f"overwriting default parameters with config file: "
                        f"{self.config_file}")
            config = load_yaml(self.config_file)

    def check(self):
        """
        Check that parameters are set acceptably
        :return:
        """

    def get_event(self):
        """
        Download event metadata from IRIS using ObsPy functionality
        :return:
        """
        if self.event_selection == "catalog":
            event = self._get_event_catalog()
        else:
            if self.client_name.upper() == "LLNL":
                event = self._get_event_llnl()
            else:
                event = self._get_event_default()

        return event

    def _get_event_catalog(self):
        """
        TODO Assert origin time and second before/after event set
        :return:
        """
        logger.info(f"getting event information with client {self.client_name}")

        cat = self.Client.get_events(
            starttime=self.origintime - self.seconds_before_event,
            endtime=self.origintime + self.seconds_after_event
        )
        event = cat[0]

        return event

    def _get_event_default(self):
        """
        Make a barebones event object based on user-defined parameters

        TODO Assert that event parameters are not None if event selection
        """
        logger.info("getting event information with user parameters")

        origin = Origin(latitude=self.event_latitude,
                        longitude=self.event_longitude,
                        depth=self.event_depth_km, time=self.origin_time
                        )
        magnitude = Magnitude(mag=self.event_magnitude, magnitude_type="Mw")

        event = Event(origins=[origin], magnitudes=[magnitude])

        return event

    def _get_event_llnl(self):
        """
        Special getter function for Lawrence Livermore National Lab data
        LLNL database has a special client

        TODO Check import llnl_db_client
        :return:
        """
        logger.info("getting event information from LLNL database")

        cat = self.client.get_catalog()
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
        self.st_raw = st
        self.inv = inv

    def get_stations(self):
        """
        Download station metadata from IRIS using ObsPy functionality

        TODO
        :return:
        """
        logger.info(f"collecting station data from {self.client_name}")
        inv = self.client.get_stations(
            network=self.network, location=self.location, station=self.station,
            channel=self.channel,
            starttime=self.reference_time - self.seconds_before_ref,
            endtime=self.reference_time + self.seconds_after_ref,
            minlatitude=self.min_lat, maxlatitude=self.max_lat,
            minlongitude=self.min_lon, maxlongitude=self.max_lon,
            level="response"
        )
        logger.info(f"collected {len(inv)} stations from {self.client_name}")

        inv = self.curtail_stations(inv)

        return inv

    def _curtail_stations(self):
        """
        Remove stations based on station distance, azimuth, etc.
        :return:
        """


    def append_sac_headers(self):
        """
        Write SAC headers to Stream objects
        :return:
        """


    def main(self):
        """
        Run Pysep seismogram extraction
        :return:
        """
        # Overload default parameters with event input file and check validity
        self.overwrite_parameters()
        self.check()

        self.event = self.get_event()

        # Get waveforms and add additional information
        self.st = self.get_waveforms()
        self.inv = self.get_stations()

        self.st = self.append_sac_headers()
        self.st = self.quality_check_waveforms()

        # Pre-filtering and preprocessing
        self.st = self.preprocess()
        self.get_unique_station_codes()

        # Making sure time series are standardized
        self.st = self.resample()
        self.st = self.trim_start_end_times()

        # Write intermediate files
        self.write()

        # Rotate streams to desired orientation
        self.rotate()

        # Output files for other programs
        self.write_weights(format="cap")
        self.write(format="all")


if __name__ == "__main__":
    # !!! Parse arguments
    sep = Pysep(config_file=args.config_file)
    sep.main()
