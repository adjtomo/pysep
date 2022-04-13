#!/usr/bin/env python3
"""
PYSEP: PYthon Seismic Extraction Program - waveform fetching and prepping.
A tool for interfacing with IRIS waveform and meta-data using ObsPy with options
for SAC input and output.
"""
import pickle
import warnings

from copy import deepcopy
from obspy.clients.fdsn import Client
from obspy.geodetics import kilometer2degrees
from obspy.core.event import Event, Origin, Magnitude, Catalog
from obspy.clients.fdsn.mass_downloader import (CircularDomain,
                                                RectangularDomain, Restrictions,
                                                MassDownloader)

from util_write_cap import *
from clipping_handler import remove_clipped
from util_helpers import get_streams_from_dir, get_inventory_from_xml

class GetWaveform:
    """
    Download, preprocess, and save waveform data from IRIS using ObsPy

    Note: station names expected in SEED format, that is
    NN.SSS.LL.CCC where: N=network, S=Station, L=Location, C=Channel

    .. note:: Embargoed Data
        For data embargoed data, register with IRIS here:
        http://ds.iris.edu/ds/nodes/dmc/forms/restricted-data-registration/
        And then set `user` and `password`

    .. note:: Filtering
        f1=fmin - highpass will keep frequencies larger than fmin
        f2=fmax - lowpass will keep frequencies lower than fmax
        f1 should consider the requested length of the time series
        f2 should consider the sampling rate for the desired channels
        'bandpass' both f1 and f2 are used
        'lowpass' only f2 is used
        'highpass' only f1 is used
        zerophase = False (causal/one-pass), = True (acausal/two-pass)
        4 pole filter is more sharper at the edges than 2 pole

    .. note:: Example uses
                                      ifFilter  zerophase  remove_response  ipre_filt
    A. CAP-ready waveforms [DEFAULT]: False     NA         True             1
    B. plot-ready waveforms, acausal: True      True       True             2
    C. plot-ready, causal waveforms:  True      False      True             0
    D. possible sensor issues:        True      False      False            NA
    """
    def __init__(self, client_name="IRIS", network="*", station="*",
                 location="*", channel="*", remove_response=True,
                 remove_clipped=False, filter_type="bandpass", f1=1/40,
                 f2=1/10, zerophase=True, corners=4,
                 water_level=60, detrend=True, demean=True, taper=False,
                 rotateRTZ=True, rotateUVW=False, rotateENZ=True,
                 min_dist=0, max_dist=20E3, min_az=0, max_az=360,
                 min_lat=None, min_lon=None, max_lat=None, max_lon=None,
                 resample_TF=False, resample_freq=50, scale_factor=1,
                 overwrite_ddir=True, phase_window=False, phases=None,
                 sec_before_after_event=20, tbefore_sec=100,
                 tafter_sec=300, use_catalog=False, ifph5=False,
                 write_sac_phase=False, taupmodel="ak135",
                 output_cap_weight_file=True, output_event_info=True,
                 outformat="VEL",  ipre_filt=1, iplot_response=False,
                 icreateNull=False, isave_raw=False, isave_raw_processed=True,
                 isave_rotateENZ=True, isave_ENZ=True,
                 ifFilter=False, ifsave_sacpaz=False, ifplot_spectrogram=False,
                 ifsave_stationxml=True, ifverbose=True, ifsave_asdf=False,
                 ifmass_downloader=False,
                 user=None, password=None
                 ):
        """
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
        self.client_name = client_name

        # Event-related parameters
        self.use_catalog = use_catalog
        self.sec_before_after_event = sec_before_after_event
        self.tbefore_sec = tbefore_sec
        self.tafter_sec = tafter_sec

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

        # Username and password for embargoed IRIS data
        self.user = user
        self.password = password

        # Output flags for how to save the raw or processed data
        self.ifph5 = ifph5
        self.output_event_info = output_event_info
        self.output_cap_weight_file = output_cap_weight_file
        self.ifsave_sacpaz = ifsave_sacpaz
        self.ifplot_spectrogram = ifplot_spectrogram
        self.iplot_response = iplot_response
        self.ifsave_stationxml = ifsave_stationxml
        self.ifsave_asdf = ifsave_asdf
        self.isave_raw = isave_raw
        self.isave_raw_processed = isave_raw_processed
        self.isave_ENZ = isave_ENZ
        self.isave_rotateENZ = isave_rotateENZ

        self.remove_clipped = remove_clipped
        if self.remove_clipped:
            print("Removing clipped stations forces `isave_raw`==True")
            self.isave_raw = True

        # Flags to deal with miscellaneous items
        self.ifverbose = ifverbose
        self.ifmass_downloader = ifmass_downloader

        # Internally used Event parameters when use_catalog == 0
        self.client = None
        self.evname = None
        self.otime = None
        self.elat = None
        self.elon = None
        self.edep = None
        self.emag = None
        # Dummy values for reference origin
        self.rlat = None
        self.rlon = None
        self.rtime = None
        self.ev = Event()
        self.ref_time_place = Event()

    def copy(self):
        """
        Create of copy of this class
        """
        return deepcopy(self)

    def run_get_waveform(self):
        """
        Get SAC waveforms for an event

        c              -  client
        event          -  obspy Event object
        ref_time_place -  reference time and place (other than origin time and place - for station subsetting)
        """
        event = self.ev
        ref_time_place = self.ref_time_place

        evtime = event.origins[0].time
        reftime = ref_time_place.origins[0].time

        if self.ifmass_downloader:
            stream_raw, inventory = self._get_waveform_mass_download()
        elif self.client_name == "LLNL":
            stream_raw, inventory = self._get_waveform_llnl()
        else:
            stream_raw, inventory = self._get_waveform()

        stream = set_reftime(stream_raw, evtime)
        
        print("--> Adding SAC metadata...")
        if self.ifverbose:
            print(stream.__str__(extended=True))
        st2 = add_sac_metadata(stream, client_name=self.client_name, ev=event, 
                               stalist=inventory, taup_model= self.taupmodel, 
                               phases=self.phases,
                               phase_write=self.write_sac_phase)
        
        do_waveform_QA(st2, self.client_name, event, evtime,
                       self.tbefore_sec, self.tafter_sec)

        st2 = self.preprocess(st2, inventory)

        # Set the sac header KEVNM with event name. This applies to the events
        # from the LLNL database
        # !!! Note: this command is needed at the time of writing files, so it
        # has to be set early
        st2, evname_key = rename_if_LLNL_event(st2, evtime)
        self.evname = evname_key

        # Note: Plotted are stations in the inventory and NOT the ones with
        # traces. It could be possible that there might not be waveforms for
        # some of these stations.
        try:
            fig = inventory.plot(projection="local", resolution="i",
                                 label=False, show=False)
            Catalog([self.ev]).plot(fig=fig,
                                    outfile=os.path.join(self.evname,
                                                         "station_map.pdf"))
        except:
            print("There is a problem with creating the station map!")

        # Get list of unique stations + locations (example: 'KDAK.00')
        stalist = []
        for tr in stream.traces:
            if self.ifverbose:
                print(f"\t{tr}")
            s = tr.stats
            station_key = \
                f"{s.network}.{s.station}.{s.location}.{s.channel[:-1]}"
            stalist.append(station_key)

        # Crazy way of getting a unique list of stations
        stalist = list(set(stalist))

        if self.resample_TF == True:
            # NOTE !!! tell the user if BOTH commands are disabled NOTE !!!
            if self.client_name == "IRIS":
                resample(st2, freq=self.resample_freq)
            elif (self.client_name == "LLNL"):
                resample_cut(st2, self.resample_freq, evtime,
                             self.tbefore_sec, self.tafter_sec)
        else:
            warnings.warn("Will not resample data")

        # Match start and end points for all traces
        st2 = trim_maxstart_minend(stalist, st2, self.client_name, event, evtime, 
                self.resample_TF, self.resample_freq, 
                self.tbefore_sec, self.tafter_sec, self.ifverbose)
        if len(st2) == 0:
            raise ValueError("no waveforms left to process!")

        # Save raw waveforms in SAC format
        if self.isave_raw:
            path_to_waveforms = evname_key + "/RAW"
            write_stream_sac_raw(stream_raw, path_to_waveforms, 
                                 evname_key, self.client_name, event,
                                 stations=inventory)


        if self.taper:
            st2.taper(max_percentage=self.taper, type='hann',max_length=None,
                      side='both')

        # save processed waveforms in SAC format
        # evname_key/RAW_processed = traces after waveform_QA + demean + detrend +
        #                            resample + remove response + filtering +
        #                            resampling + scaling + tapering
        # NOTE: The orientation is same as that of extracted waveforms
        #       Waveforms are rotated to ENZ, in case they are not already orientated,
        #       in the next step (self.rotateRTZ)
        if self.isave_raw_processed:
            path_to_waveforms = os.path.join(evname_key, 'RAW_processed')
            write_stream_sac(st2, path_to_waveforms, evname_key)

        # Rotate to ENZ (save: optional)
        #if self.rotateENZ:
        #st2 = rotate2ENZ(st2, evname_key, self.isave_ENZ, self.icreateNull,
        # self.ifverbose)

        if self.rotateENZ:
            st2 = rotate2ENZ(st2, evname_key, self.isave_ENZ, self.icreateNull,
                             self.ifverbose)

        # rotate to UVW and save
        if self.rotateUVW:
            rotate2UVW(st2, evname_key) 

        # Rotate to RTZ and save
        if self.rotateRTZ:
            rotate2RTZ(st2, evname_key, self.ifverbose) 

        # save CAP weight files
        if self.output_cap_weight_file:
            write_cap_weights(st2, evname_key, self.client_name, event,
                              self.ifverbose)

        # save event info
        if self.output_event_info:
            write_ev_info(event, evname_key)

        # Plot spectrograms
        if self.ifplot_spectrogram:
            plot_spectrogram(st2, evname_key)

        # save pole zero file (Needed for MouseTrap)
        if self.ifsave_sacpaz:
            write_resp(inventory,evname_key)

        # save station inventory as XML file
        if self.ifsave_stationxml:
            xmlfilename = evname_key + "/stations.xml"
            try:
                inventory.write(xmlfilename, format="stationxml", validate=True)
            except:
                print('Could not create stationxml file')
        
        # Path to the asdf_converter script        
        if self.ifsave_asdf:
            # save RTZ
            asdf_filename = evname_key + "/" + evname_key + ".h5"
            os.system("../asdf_converters/asdf_converters/sac2asdf.py "
                      + evname_key + " " + asdf_filename + " observed")
            # save NEZ
            nez_dir = evname_key + "/ENZ/"
            nez_asdf_filename = nez_dir + evname_key + ".h5"
            os.system("../asdf_converters/asdf_converters/sac2asdf.py "
                      + nez_dir + " " + nez_asdf_filename + " observed")
            
        if self.remove_clipped:
            remove_clipped(evname_key)

    def _get_waveform(self):
        """
        Get waveform for any ObsPy client
        """
        c = self.client

        reftime = self.ref_time_place.origins[0].time
        print(f"DATABASE >>> Sending request to {self.client_name} client "
              f"for data")

        # Client NCEDC does not understand '-XXX' station code
        if self.client_name == "NCEDC":
            assert("-" not in self.station), ("NCEDC client does not take '-' "
                                              "in station code")
        elif self.client_name == "IRIS":
            if '*' in self.network:
                warnings.warn("You have chosen to search ALL networks at "
                              "IRIS. This could take long!")
        if self.ifph5:
            STATION = 'http://service.iris.edu/ph5ws/station/1'
            c = Client('http://service.iris.edu',
                       service_mappings={'station': STATION}, debug=True
                       )

        print("Downloading stations...")
        stations = c.get_stations(network=self.network, location=self.location,
                                  station=self.station, channel=self.channel,
                                  starttime=reftime - self.tbefore_sec,
                                  endtime=reftime + self.tafter_sec,
                                  minlatitude=self.min_lat,
                                  maxlatitude=self.max_lat,
                                  minlongitude=self.min_lon,
                                  maxlongitude=self.max_lon,
                                  level="response")
        inventory = stations  # so that llnl and iris scripts can be combined

        if self.ifverbose:
            for sta in stations:
                print(f"\t{sta}")

        sta_limit_distance(self.ref_time_place, stations,
                           min_dist=self.min_dist, max_dist=self.max_dist,
                           min_az=self.min_az, max_az=self.max_az,
                           ifverbose=self.ifverbose)
        # print("Printing stations NEW")
        # print(stations)
        # print("Done Printing stations...")

        # stations.plotprojection="local")
        # Find P and S arrival times
        t1s, t2s = get_phase_arrival_times(stations, self.ev, self.phases,
                                           self.phase_window, self.taupmodel,
                                           reftime, self.tbefore_sec,
                                           self.tafter_sec)

        print("Downloading waveforms...")
        if self.ifph5:
            DATASELECT = 'http://service.iris.edu/ph5ws/dataselect/1'
            c = Client('http://service.iris.edu',
                       service_mappings={'dataselect': DATASELECT},
                       user=self.user, password=self.password, debug=True
                       )
            stream_raw = c.get_waveforms(network=self.network,
                                         location=self.location,
                                         station=self.station,
                                         channel=self.channel,
                                         starttime=reftime - self.tbefore_sec,
                                         endtime=reftime + self.tafter_sec
                                         )
        else:
            bulk_list = make_bulk_list_from_stalist(stations, t1s, t2s,
                                                    channel=self.channel)
            stream_raw = c.get_waveforms_bulk(bulk_list)

        # Save self as a pickle  (why do we want to do this? -B)
        file_out = os.path.join(self.evname, f"{self.evname}_ev_info.obj")
        with open(file_out, "wb") as f:
            pickle.dump(self, f)

        return stream_raw, inventory

    def _get_waveform_mass_download(self):
        """
        Apply ObsPy mass downloader to grab events from a pre-determined region
        """
        # Determine origin time and reference time in the given region
        evtime = self.ev.origins[0].time
        reftime = self.ref_time_place.origins[0].time


        mdl = MassDownloader()
        domain = CircularDomain(latitude=self.elat, longitude=self.elon,
                                minradius=kilometer2degrees(self.min_dist),
                                maxradius=kilometer2degrees(self.max_dist)
                                )
        res = Restrictions(starttime=reftime - self.tbefore_sec,
                           endtime=reftime + self.tafter_sec,
                           station_starttime=None,
                           station_endtime=None,
                           chunklength_in_sec=None,
                           network=self.network, station=self.station,
                           location=self.location, channel=self.channel,
                           reject_channels_with_gaps=False, minimum_length=0.0,
                           sanitize=True, minimum_interstation_distance_in_m=0,
                           # !!! These were commented, why? - B
                           # exclude_networks = (), exclude_stations = (),
                           # limit_stations_to_inventory=None,
                           # channel_priorities=(),
                           # location_priorities=()
                           )

        # Define the path structure of output files
        outdir = os.path.join(os.getcwd(), self.evname)
        mseed_storage = os.path.join(outdir, "mass_downloader", "waveforms")
        stationxml_storage = os.path.join(outdir, "mass_downloader", "stations")


        mdl.download(domain, res, mseed_storage=mseed_storage,
                     stationxml_storage=stationxml_storage,
                     download_chunk_size_in_mb=20, threads_per_client=3,
                     print_report=True)

        inventory = get_inventory_from_xml(stationxml_storage)
        stream_raw = get_streams_from_dir(mseed_storage)

        phases = self.phases

        # !!! This isnt used? - B
        t1s, t2s = get_phase_arrival_times(inventory, self.ev, self.phases,
                                           self.phase_window, self.taupmodel,
                                           reftime, self.tbefore_sec,
                                           self.tafter_sec)

        return stream_raw, inventory

    def _get_waveform_llnl(self):
        """
        Get waveforms from Lawrence Livermore National Laboratory
        """
        print("Preparing request for LLNL ...")
        c = self.client

        # Get event an inventory from the LLNL DB.
        event_number = int(self.ev.event_descriptions[0].text)
        # event = llnl_db_client.get_obspy_event(event)
        inventory = c.get_inventory()

        nsta_llnl = len(inventory.get_contents()["stations"])

        print(f"\tTotal stations in LLNL DB: {nsta_llnl}")
        sta_limit_distance(self.ev, inventory,
                           min_dist=self.min_dist,
                           max_dist=self.max_dist,
                           min_az=self.min_az,
                           max_az=self.max_az)
        print(f"\tStations after filtering for distance: "
              f"{len(inventory.get_contents()['stations'])}")

        stations = set([sta.code for net in inventory for sta in net])

        _st = c.get_waveforms_for_event(event_number)
        stream_raw = obspy.Stream()
        for tr in _st:
            if tr.stats.station in stations:
                stream_raw.append(tr)

        return stream_raw, stations

    def preprocess(self, st, inv):
        """
        Preprocess waveforms and return preprocessed waveforms
        """
        st2 = st.copy()

        if self.demean:
            st2.detrend('demean')
        if self.detrend:
            st2.detrend('linear')

        if self.ifFilter:
            prefilter(st2, self.f1, self.f2,
                      self.zerophase, self.corners, self.filter_type)

        if self.remove_response:
            resp_plot_remove(st2, self.ipre_filt, self.pre_filt,
                             self.iplot_response, self.water_level,
                             self.scale_factor,
                             inv, self.outformat, self.ifverbose)
        else:
            # output RAW waveforms
            decon = False
            print("WARNING -- NOT correcting for instrument response")

        if self.scale_factor > 0:
            amp_rescale(st2, self.scale_factor)
            if self.client_name == "LLNL":
                amp_rescale_llnl(st2, self.scale_factor)

        return st2


    def reference_time_place(self):
        '''
        returns an event object with different origin time and location 
        (i.e. not centered around the earthquake). Stations will be subsetted
        based on reference origin time and location
        '''

        self.ref_time_place = self.ev.copy()
        self.ref_time_place.origins[0].latitude = self.rlat
        self.ref_time_place.origins[0].longitude = self.rlon
        self.ref_time_place.origins[0].time = self.rtime

    def get_event_object(self):
        '''
        update events otime,lat,lon and mag with IRIS (or any other clients) catalog
        '''
        
        # get parameters from the cataog
        if self.use_catalog == 1:
            print("WARNING using event data from the IRIS catalog")
            cat = self.client_name.get_events(
                starttime = self.otime - self.sec_before_after_event,
                endtime = self.otime + self.sec_before_after_event)
            self.ev = cat[0]
            
            # use catalog parameters
            self.otime = self.ev.origins[0].time
            self.elat = self.ev.origins[0].latitude
            self.elon = self.ev.origins[0].longitude
            self.edep = self.ev.origins[0].depth
            self.emag = self.ev.magnitudes[0].mag
            
        # use parameters from the input file
        else:
            print("WARNING using event data from user-defined catalog")
            #self.ev = Event()
            org = Origin()
            org.latitude = self.elat
            org.longitude = self.elon
            org.depth = self.edep
            org.time = self.otime
            mag = Magnitude()
            mag.mag = self.emag
            mag.magnitude_type = "Mw"
            self.ev.origins.append(org)
            self.ev.magnitudes.append(mag)
    
    def get_events_client(self):
        '''
        get absolute and reference event object
        '''

        # IRIS
        if self.client_name != "LLNL":
            # import functions to access waveforms
            if not self.user and not self.password:
                self.client = Client(self.client_name,debug=True,
                                          timeout=600)
            else:
                self.client = Client(self.client_name, user=self.user,
                                     password=self.password, debug=True,
                                          timeout=600)
                # will only work for events in the 'IRIS' catalog
                # (future: for Alaska events, read the AEC catalog)

            # get event object
            self.get_event_object()

            # use a different reference time and place for station subsetting
            if self.rlat is not None:
                self.reference_time_place()
            # or use reference same as the origin
            else:
                self.ref_time_place = self.ev
        
        # LLNL
        if self.client_name == "LLNL":
            import llnl_db_client
            self.client = llnl_db_client.LLNLDBClient(
                "/store/raw/LLNL/UCRL-MI-222502/westernus.wfdisc")

            # get event time and event ID
            cat = self.client.get_catalog()
            mintime_str = "time > %s" % (self.otime - 
                                         self.sec_before_after_event)
            maxtime_str = "time < %s" % (self.otime + 
                                         self.sec_before_after_event)
            print(mintime_str + "\n" + maxtime_str)

            self.ev = cat.filter(mintime_str, maxtime_str)
        
            if len(self.ev) > 0:
                self.ev = self.ev[0]
                # Nothing happens here.  We can change later
                self.ref_time_place = self.ev
                print(len(self.ev))
            else:
                print("WARNING. No events in the catalog for the given time period")
                #sys.exit()

        # print client and event info
        print(self.ev)

    def save_extraction_info(self):
        # track git commit
        os.system('git log | head -12 > ./' + self.evname + 
                  '/' + self.evname + '_last_2git_commits.txt')
        # save filenames in a file for checking
        fname = self.evname + '/' + self.evname + '_all_filenames'
        fcheck = open(fname,'w')

        os.system('ls -1 ' + self.evname + '/* > ' +  fname)
