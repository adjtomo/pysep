"""
Python Seismogram Extraction and Processing

Download, pre-process, and organize seismic waveforms, event and station
metadata
"""
import os
import numpy as np
from obspy import UTCDateTime, Stream
from obspy.clients.fdsn import Client
from obspy.core.event import Event, Origin, Magnitude, Catalog

from pysep import logger
from pysep.utils.read import load_yaml
from pysep.utils.capsac import append_sac_headers, write_cap_weights_file
from pysep.utils.curtail import (curtail_by_station_distance_azimuth,
                                 quality_check_waveforms)
from pysep.utils.fetch import fetch_bulk_station_list_taup
from pysep.utils.llnl import scale_llnl_waveform_amplitudes
from pysep.utils.process import zerophase_chebychev_lowpass_filter


class Pysep:
    """
    Download, preprocess, and save waveform data from IRIS using ObsPy
    """
    def __init__(self, config_file=None,  client="IRIS", origin_time=None,
                 networks="*", stations="*", locations="*", channels="*",
                 remove_response=True, remove_clipped=False,
                 filter_type="bandpass", f1=1/40, f2=1/10, zerophase=True,
                 corners=4, water_level=60, detrend=True, demean=True,
                 taper=False, rotateRTZ=True, rotateUVW=False, rotateENZ=True,
                 min_dist=0, max_dist=20E3, min_az=0, max_az=360,
                 min_lat=None, min_lon=None, max_lat=None, max_lon=None,
                 resample_freq=50, scale_factor=1, phase_list=None,
                 seconds_before_event=20, seconds_after_event=20,
                 seconds_before_ref=100, seconds_after_ref=300,
                 write_sac_phase=False, taupmodel="ak135",
                 output_cap_weight_file=True, output_event_info=True,
                 output_unit="VEL",  user=None, password=None, **kwargs):
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
        self.client = client
        self.c = self.get_client(user=user, password=password)  # actual Client

        # Parameters related to event selection
        self.event_selection = event_selection

        self.origin_time = UTCDateTime(origin_time)
        self.seconds_before_event = seconds_before_event
        self.seconds_after_event = seconds_after_event

        # Optional: if User wants to define an event on their own
        self.event_latitude = event_latitude
        self.event_longitude = event_longitude
        self.event_depth_km = event_depth_km
        self.event_magnitude = magnitude

        # Waveform and StationXML gathering parameters
        self.networks = networks
        self.stations = stations
        self.channels = channels
        self.locations = locations

        # Waveform collection parameters
        try:
            # Reference time allows throwing a static time shift from the origin
            # time, e.g., if there are timing errors with relation to o time
            self.reference_time = UTCDateTime(reference_time)
        except TypeError:
            self.reference_time = self.origin_time
        self.seconds_before_ref = seconds_before_ref
        self.seconds_after_ref = seconds_after_ref
        self.phase_list = phase_list

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
        self.remove_response = bool(remove_response)
        self.output_unit = output_unit
        self.water_level = water_level
        self.pre_filt = pre_filt
        self.scale_factor = scale_factor
        self.resample_freq = resample_freq

        # Client related parameters
        self.debug = bool(debug)
        self.timeout = timeout

        # Internally filled attributes
        self.output_dir = None


    def check(self):
        """
        TODO
            * Check output format
            * Check preprocessing flags
            * Check user and password if PH5 data gathering
        """
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


        # Check preprocessing flags
        if self.remove_response:
            for par in [self.output_unit, self.pre_filt, self.water_level]:
                pass
        assert(0 <= self.min_az <= 360), f"0 <= `min_az` <= 360"
        assert(0 <= self.max_az <= 360), f"0 <= `min_az` <= 360"



    def get_client(self, user=None, password=None):
        """
        Options to choose different Clients based on client name. These clients
        are ONLY for waveform gathering. We assume that

        TODO add log level to parameters and use to set client debug
        :return:
        """
        if self.client is not None:
            # Lawrence Livermore Natinoal Lab internal waveform database
            if self.client.upper() == "LLNL":
                try:
                    import llnl_db_client
                    c = llnl_db_client.LLNLDBClient(
                        "/store/raw/LLNL/UCRL-MI-222502/westernus.wfdisc"
                    )
                except ImportError as e:
                    logger.warning("Cannot import llnl_db_client, please "
                                   "download from "
                                   "https://github.com/krischer/llnl_db_client")
                    raise ImportError from e
            # IRIS DMC PH5WS Station Web Service
            elif self.client.upper() == "PH5":
                c = Client(
                    "http://service.iris.edu", debug=self.debug,
                    service_mappings={
                        "station": "http://service.iris.edu/ph5ws/station/1",
                        "dataselect":
                            "http://service.iris.edu/ph5ws/dataselect/1"
                    }, user=user, password=password
                )
            # Default ObsPy FDSN webservice client
            else:
                c = Client(self.client, user=user, password=password,
                           debug=self.debug, timeout=self.timeout)
        else:
            c = None

        return c

    def overwrite_parameters(self):
        """
        Overwrite default parameters using a YAML config file
        """
        if self.config_file is not None:
            logger.info(f"overwriting default parameters with config file: "
                        f"{self.config_file}")
            config = load_yaml(self.config_file)
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

        cat = self.c.get_events(
            starttime=self.origin_time - self.seconds_before_event,
            endtime=self.origin_time + self.seconds_after_event,
            debug=self.debug
        )
        logger.info(f"{len(cat)} matching events found for "
                    f"{self.seconds_before_event}s < {self.origin_time} "
                    f"< {self.seconds_after_event}, choosing first")

        event = cat[0]

        return event

    def _get_event_default(self):
        """
        Make a barebones event object based on user-defined parameters

        TODO Assert that event parameters are not None if event selection
        TODO resource ids and preferred ids for magnitude and origin
        """
        logger.info("creating event metadata with user parameters")

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

        cat = self.c.get_catalog()
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

            st = self.c.get_waveforms_bulk(bulk=bulk)

        logger.info(f"{len(st)} waveforms returned after query")

        return st

    def _get_waveforms_ph5(self):
        """
        Download waveforms from PH5 dataservice
        :return:
        """
        logger.debug("")
        st = self.c.get_waveforms(
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
        t2 = self.reference_time - self.seconds_after_ref
        for net in self.inv:
            for sta in net:
                # net sta loc cha t1 t2
                bulk.append((net.code, sta.code, "*", self.channels, t1, t2))
        return bulk

    def get_codes(self, st=None, choice=None):
        """
        Get station codes from the internal stream attribute, where station
        codes are formatted NN.SSS.LL.CCc where N=network, S=station,
        L=location, C=channel, and c=component

        :type choice: str
        :param choice: choice of the part of the code returned, available:
            * 'network': return unique network codes (e.g., NN)
            * 'station': return unique network + station codes (e.g., NN.SSS)
            * 'location': return up to location (e.g., NN.SSS.LL)
            * 'channel': return up to channel, no component (e.g., NN.SSS.LL.CC)
            * else: return full station code (e.g., NN.SSS.LL.CCc)
        :rtype: list
        :return: unique station codes filtered by choice
        """
        if st is None:
            st = self.st.copy()

        full_codes = [tr.get_id() for tr in st]
        if choice == "network":
            codes = [code.split(".")[0] for code in full_codes]
        elif choice == "station":
            codes = [".".join(code.split(".")[:1]) for code in full_codes]
        elif choice == "location":
            codes = [".".join(code.split(".")[:2]) for code in full_codes]
        elif choice == "channel":
            codes = [code[:-1] for code in full_codes]
        else:
            codes = full_codes

        return list(set(codes))


    def get_stations(self):
        """
        Download station metadata from IRIS using ObsPy functionality
        :return:
        """
        logger.info(f"collecting station data from {self.client}")

        inv = self.c.get_stations(
            network=self.networks, location=self.locations,
            station=self.stations, channel=self.channels,
            starttime=self.origin_time - self.seconds_before_ref,
            endtime=self.origin_time + self.seconds_after_ref,
            minlatitude=self.min_lat, maxlatitude=self.max_lat,
            minlongitude=self.min_lon, maxlongitude=self.max_lon,
            level="response", debug=self.debug
        )

        logger.info(f"collected {len(inv)} stations from {self.client}")

        inv = self.curtail_stations(inv)
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
        return inv

    def preprocess(self):
        """
        Very simple preprocessing to remove response and apply a prefilter

        TODO do we need the pre-filt if we are doing it in resample?
        """
        st_out = self.st.copy()
        if self.demean:
            logger.info(f"applying demean to all data")
            st_out.detrend("demean")
        if self.detrend
            logger.info(f"applying linear detrend to all data")
            st_out.detrend("linear")
        if self.remove_response:
            logger.info(f"removing response, output units in: "
                        f"{self.output_unit}")
            st_out.remove_response(inv=self.inv, water_level=self.water_level,
                                   pre_filt=self.pre_filt,
                                   output=self.output_unit)
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
        codes = self.get_codes(st=st_edit)
        resampled_codes = []

        for code in codes:
            # Brute force check if we've done this station, so we don't do again
            break_loop = False
            for resampled in resampled_codes:
                if code.startswith(resampled):
                    break_loop = True
            if break_loop:
                break

            # Subset the stream based on the N number of components
            net, sta, loc, cha = code.split(".")
            cha_wild = f"{cha[:-1]}*"  # replace component with wildcard
            st_edit_select = st_edit.select(network=net, station=sta,
                                            location=loc, channel=cha_wild)

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
            resampled_codes.append(f"{code[-1]}")  # drop component

        return st_out

    def write(self, fmt, fid):
        """
        Write out various files specifying information about the collected
        stations and waveforms

        :param fmt:
        :return:
        """

    def _write_station_list(self, fid):
        """
        Write a list of station codes, distances, etc.
        """
        with open(fid, "w") as f:
            for net in self.inv:
                for sta in net:
                    netsta_code = f"{net.code}.{sta.code}"
                    f.write(f"{sta.code} {net.code} {sta.latitude} "
                            f"{sta.longitude} {self.distances[netsta_code]}"
                            f"{self.azimuths[netsta_code]}")

    def main(self):
        """
        Run Pysep seismogram extraction
        :return:
        """
        # Overload default parameters with event input file and check validity
        self.overwrite_parameters()
        self.check()

        # Get waveforms and metadata (QuakeML, StationXML)
        self.event = self.get_event()
        self.st = self.get_waveforms()
        self.inv = self.get_stations()
        self.inv = self.curtail_stations()

        self.st = append_sac_headers(self.st, self.inv, self.event)
        self.output_dir = self.st[0].stats.sac["kevnm"]
        self.st = quality_check_waveforms()

        # Pre-filtering and preprocessing
        self.st = self.preprocess()
        self.st = self.resample_and_trim_start_end_times()
        self.rotate()

        # Write intermediate files
        self.write()

        # Rotate streams to desired orientation

        # Output files for other programs
        write_cap_weights_file(st=self.st, event=self.event,
                               path_out=self.output_dir, order_by="dist")
        self.write(format="all")
        # plot_station_map(self.inv, self.event)


if __name__ == "__main__":
    # !!! Parse arguments
    sep = Pysep(config_file=args.config_file)
    sep.main()
