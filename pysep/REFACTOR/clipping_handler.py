import os
import obspy



def remove_clipped(path):
    ddir = path+'/'
    rawdir = os.path.join(ddir, 'RAW/')
    # Load all sac files as one stream
    data = obspy.read(os.path.join(rawdir,'*.sac'))
    # Define clipping factor for 24 bits signal
    q = ((2**(24-1))**2)**0.5

    # Empty list to contain all the clipped traces ids
    clipped_stations = []
    # Loop in traces stream
    for trace in data:
        n = len(trace.data[(trace.data**2)**(0.5) > 0.8*q])
        if not n == 0:
            clipped_stations.append(trace.id[:-1])
            print('station', trace.id[:-1], 'has', n, 'clipped values!')

    # Replacing surface wave weight in .dat files
    for file in os.listdir(ddir):
        if file.endswith('.dat'):
            with open(os.path.join(ddir, file), 'r') as read_obj:
                with open(os.path.join(ddir, file+'mod'), 'w') as tfile:
                    for line in read_obj:
                            if any(station in line for station in clipped_stations):
                                line = line.replace('1 1 1', '0 0 0')
                                tfile.writelines(line)
                            else:
                                tfile.writelines(line)
            os.replace(os.path.join(ddir, file+'mod'), os.path.join(ddir, file))



# DELETE ME

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
        self.seconds_before_ref = tbefore_sec
        self.seconds_after_ref = tafter_sec

        # Station-related parameters
        self.networks = network
        self.stations = station
        self.stations = station
        self.channels = channel
        self.locations = location
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
        self.output_unit = output_unit
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
        self.event = None
        self.inv = None
        self.st = None
        self.distances = None
        self.azimuths = None