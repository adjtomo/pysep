#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tools for interfacing IRIS data, ObsPy, and SAC input/output.
"""
from __future__ import print_function

import os
from copy import deepcopy

import obspy
from obspy.clients.fdsn import Client
from scipy import signal

from util_write_cap import *

from obspy.core.event import Event, Origin, Magnitude

class getwaveform:
    def __init__(self):
        """
        ---------------- copied from old getwaveform.py ----------------------
        min_dist - minimum station distance (default = 20)
        max_dist - maximum station distance (default =300)
        before -   time window length before the event time (default= 100)
        after  -   time window length after the event time (default = 300)
        network -  network codes of potential stations (default=*)
        channel -  component(s) to get, accepts comma separated (default='BH*')
        ifresample_TF   - Boolean. Request resample or not. Default = False
        resample_freq   - sampling frequency to resample waveforms (default 20.0)
        ifrotate - Boolean, if true will output sac files rotated to baz
        unrotated sac files will also be written
        ifCapInp - Boolean, make weight files for CAP
        ifEvInfo - Boolean, output 'ev_info.dat' containg event info (True)
        ifRemoveResponse - Boolean, will remove response (True)
        ifDetrend - Boolean, will remove linear trend from data (True)
        ifDemean  - Boolean, will insult the data (True)
        scale_factor - scale all data by one value (10.0**2)
        This usually puts the data in the units required by CAP
        From m/s to cm/s
        pre_filt  - list, corner frequencies of filter to apply before deconv
        a good idea when deconvolving (ifRemoveResponse=True)
        """
        # DEFAULT SETTINGS (see getwaveform.py)
        self.idb = 1    # default: =1-IRIS; =2-AEC; =3-LLNL

        # event parameters
        self.use_catalog = 1              # use an existing catalog (=1) or specify your own event parameters (see iex=9)
        self.sec_before_after_event = 10  # time window to search for a target event in a catalog
        self.min_dist = 0 
        self.max_dist = 20000
        self.min_az = 0 
        self.max_az = 360
        self.tbefore_sec = 100
        self.tafter_sec = 300
        # These are used only if self.use_catalog = 0
        self.otime = None
        self.elat = None
        self.elon = None
        self.edep = None
        self.emag = None

        # station parameters
        self.network = '*'                # all networks
        self.station = '*,-PURD,-NV33,-GPO'  # all stations
        self.channel = '*'                # all channels     
        self.overwrite_ddir = 1           # 1 = delete data directory if it already exists
        self.icreateNull = 1              # create Null traces so that rotation can work (obsby stream.rotate require 3 traces)
        
        # Filter parameters
        self.ifFilter = False 
        #------Filter--------------
        # for 'bandpass' both f1 and f2 are used
        # for 'lowpass' only f2 is used
        # for 'highpass' only f1 is used
        #
        # EXAMPLES
        #                                   ifFilter  zerophase  remove_response  ipre_filt
        # A. CAP-ready waveforms [DEFAULT]: False     NA         True             1
        # B. plot-ready waveforms:          True      True       True             2
        # C. plot-ready, causal waveforms:  True      False      True             0
        # D. possible sensor issues:        True      False      False            NA
        #
        self.filter_type = 'bandpass'
        # f1 should consider the requested length of the time series
        # f2 should consider the sampling rate for the desired channels
        self.f1 = 1/40 # fmin - highpass will keep frequencies larger than fmin
        self.f2 = 1/10 # fmax - lowpass will keep frequencies lower than fmax
        self.zerophase = True             # = False (causal), = True (acausal)
        # 4 pole filter is more sharper at the edges than 2 pole
        self.corners = 4                  # Is corner in Obspy same as Pole in SAC?
        
        # Pre-filter parameters
        self.ipre_filt = 1                # =0 No pre_filter
                                          # =1 default pre_filter (see getwaveform_iris.py)
                                          # =2 user-defined pre_filter (use this if you are using bandpass filter)
        # For tapering down the pre-filter
        # Perhaps you want to set ipre_filt = 0 to prevent extra filtering
        # pre-filter for deconvolution
        # https://ds.iris.edu/files/sac-manual/commands/transfer.html
        # Pre-filter will not be applied if remove_response = False 
        self.f0 = 0.5*self.f1
        self.f3 = 2.0*self.f2
        self.pre_filt=(self.f0, self.f1, self.f2, self.f3)    # applies for ipre_filt = 2 only
        # self.pre_filt = (0.005, 0.006, 10.0, 15.0) # BH default

        # For CAP
        self.resample_TF = True           # if False then resample_freq is taken from SAC files
        self.resample_freq = 50           # 0 causes errors. Use resample_TF instead
        self.scale_factor = 10**2         # for CAP use 10**2  (to convert m/s to cm/s)

        # Pre-processing (manily for CAP)
        self.rotateRTZ = True
        self.rotateUVW = False            # This option works only if 'rotateRTZ = True'
        self.output_cap_weight_file = True
        self.detrend = True
        self.demean = True
        self.taper = False                # this could also be a fraction between 0 and 1 (fraction to be tapered from both sides)
        self.output_event_info = True
        self.outformat = 'VEL'            # Intrument response removed waveforms could be saved as 'VEL' 'DISP' 'ACC'
        self.ifsave_sacpaz = False        # save sac pole zero (needed as input for MouseTrap module)
        self.ifEvInfo=True
        self.remove_response = True       # remove instrument response 
        self.iplot_response = False       # plot response function
        self.ifplot_spectrogram = False   # plot spectrograms 

        # Refernce origin (dummy values)
        self.rlat = None
        self.rlon = None
        self.rtime = None

    def run_get_waveform(self,c, event, ref_time_place):
        """
        Get SAC waveforms for an event
        
        basic usage:
        run_get_waveform(event)
        
        c              -  client
        event          -  obspy Event object
        ref_time_place -  reference time and place (other than origin time and place - for station subsetting)
        """
        
        evtime = event.origins[0].time
        reftime = ref_time_place.origins[0].time
        
        if self.idb==1:
            print("Preparing request for IRIS ...")
            # BK network doesn't return data when using the IRIS client.
            # this option switches to NCEDC if BK is 
            if "BK" in self.network:
                client_name = "NCEDC"
                print("\nWARNING. Request for BK network. Switching to NCEDC client")
                c = Client("NCEDC")
            else:
                client_name = "IRIS" 
                
                print("Download stations...")
                stations = c.get_stations(network=self.network, station=self.station, 
                                          channel=self.channel,
                                          starttime=reftime - self.tbefore_sec, endtime=reftime + self.tafter_sec,
                                          level="response")
                inventory = stations    # so that llnl and iris scripts can be combined
                print("Printing stations")
                print(stations)
                print("Done Printing stations...")
                sta_limit_distance(ref_time_place, stations, min_dist=self.min_dist, max_dist=self.max_dist, min_az=self.min_az, max_az=self.max_az)
                
                print("Downloading waveforms...")
                bulk_list = make_bulk_list_from_stalist(
                    stations, reftime - self.tbefore_sec, reftime + self.tafter_sec, channel=self.channel)
                stream_raw = c.get_waveforms_bulk(bulk_list)
                
        elif self.idb==3:
            client_name = "LLNL"
            print("Preparing request for LLNL ...")
            
            # Get event an inventory from the LLNL DB.
            event_number = int(event.event_descriptions[0].text)
            # event = llnl_db_client.get_obspy_event(event)
            inventory = c.get_inventory()
            
            print("--> Total stations in LLNL DB: %i" % (
                    len(inventory.get_contents()["stations"])))
            sta_limit_distance(event, inventory, min_dist=self.min_dist, max_dist=self.max_dist, min_az=self.min_az, max_az=self.max_az)
            print("--> Stations after filtering for distance: %i" % (
                    len(inventory.get_contents()["stations"])))
            
            stations = set([sta.code for net in inventory for sta in net])
            
            _st = c.get_waveforms_for_event(event_number)
            stream_raw = obspy.Stream()
            for tr in _st:
                if tr.stats.station in stations:
                    stream_raw.append(tr)
    
        # set reftime
        stream = obspy.Stream()
        stream = set_reftime(stream_raw, evtime)
        
        print("--> Adding SAC metadata...")
        st2 = add_sac_metadata(stream,idb=self.idb, ev=event, stalist=inventory)
        
        # Do some waveform QA
        # - (disabled) Throw out traces with missing data
        # - log waveform lengths and discrepancies
        # - Fill-in missing data -- Carl request
        do_waveform_QA(st2, client_name, event, evtime, self.tbefore_sec, self.tafter_sec)
        
        if self.demean:
            st2.detrend('demean')
            
        if self.detrend:
            st2.detrend('linear')
            
        if self.ifFilter:
            prefilter(st2, self.f1, self.f2, self.zerophase, self.corners, self.filter_type)
            
        if self.remove_response:
            resp_plot_remove(st2, self.ipre_filt, self.pre_filt, self.iplot_response, self.scale_factor, stations, self.outformat)
        else:
            # output RAW waveforms
            decon=False
            print("WARNING -- NOT correcting for instrument response")

        if self.scale_factor > 0:
            amp_rescale(st2, self.scale_factor)
            if self.idb ==3:
                amp_rescale_llnl(st2, self.scale_factor)


        # Set the sac header KEVNM with event name
        # This applies to the events from the LLNL database
        # NOTE this command is needed at the time of writing files, so it has to
        # be set early
        st2, evname_key = rename_if_LLNL_event(st2, evtime)

        # Get list of unique stations + locaiton (example: 'KDAK.00')
        stalist = []
        for tr in st2.traces:
            # stalist.append(tr.stats.station)
            stalist.append(tr.stats.network + '.' + tr.stats.station +'.'+ tr.stats.location + '.'+ tr.stats.channel[:-1])

        # Crazy way of getting a unique list of stations
        stalist = list(set(stalist))

        # match start and end points for all traces
        st2 = trim_maxstart_minend(stalist, st2, client_name, event, evtime, self.resample_TF, self.resample_freq, self.tbefore_sec, self.tafter_sec)
        if len(st2) == 0:
            raise ValueError("no waveforms left to process!")

        if self.resample_TF == True:
            # NOTE !!! tell the user if BOTH commands are disabled NOTE !!!
            if (client_name == "IRIS"):
                resample(st2, freq=self.resample_freq)
            elif (client_name == "LLNL"):
                resample_cut(st2, self.resample_freq, evtime, self.tbefore_sec, self.tafter_sec)
        else:
            print("WARNING. Will not resample. Using original rate from the data")

        # save raw waveforms in SAC format
        path_to_waveforms = evname_key + "/RAW"
        write_stream_sac_raw(stream_raw, path_to_waveforms, evname_key, self.idb, event, stations=inventory)

        # Taper waveforms (optional; Generally used when data is noisy- example: HutchisonGhosh2016)
        # https://docs.obspy.org/master/packages/autogen/obspy.core.trace.Trace.taper.html
        # To get the same results as the default taper in SAC, use max_percentage=0.05 and leave type as hann.
        if self.taper:
            st2.taper(max_percentage=self.taper, type='hann',max_length=None, side='both')

        # save processed waveforms in SAC format
        path_to_waveforms = evname_key 
        write_stream_sac(st2, path_to_waveforms, evname_key)
        
        if self.rotateRTZ:
            rotate_and_write_stream(st2, evname_key, self.icreateNull, self.rotateUVW)

        if self.output_cap_weight_file:
            write_cap_weights(st2, evname_key, client_name, event)

        if self.ifEvInfo:
            write_ev_info(event, evname_key)

        if self.ifplot_spectrogram:
            plot_spectrogram(st2, evname_key)

        if self.ifsave_sacpaz:
            write_resp(inventory,evname_key)

        # save station inventory as XML file
        xmlfilename = evname_key + "/stations.xml"
        stations.write(xmlfilename, format="stationxml", validate=True)

    def copy(self):
        '''
        create of copy of itself
        '''
        return deepcopy(self)

    def reference_time_place(self,ev):
        '''
        returns an event object with different origin time and location (i.e. not centered around the earthquake). 
        Stations will be subsetted based on reference origin time and location
        '''

        ref_time_place = ev
        ref_time_place.origins[0].latitude = self.rlat
        ref_time_place.origins[0].longitude = self.rlon
        ref_time_place.origins[0].time = self.rtime

        return(ref_time_place)

    def get_event_object(self,c):
        '''
        update events otime,lat,lon and mag with IRIS (or any other clients) catalog
        '''
        
        # get parameters from the cataog
        if self.use_catalog == 1:
            print("WARNING using event data from the IRIS catalog")
            cat = c.get_events(starttime = self.otime - self.sec_before_after_event,\
                                        endtime = self.otime + self.sec_before_after_event)
            ev = cat[0]
            
            # use catalog parameters
            self.otime = ev.origins[0].time
            self.elat = ev.origins[0].latitude
            self.elon = ev.origins[0].longitude
            self.edep = ev.origins[0].depth
            self.emag = ev.magnitudes[0].mag
            
        # use parameters from the input file
        else:
            print("WARNING using event data from user-defined catalog")
            ev = Event()
            org = Origin()
            org.latitude = self.elat
            org.longitude = self.elon
            org.depth = self.edep
            org.time = self.otime
            mag = Magnitude()
            mag.mag = self.emag
            mag.magnitude_type = "Mw"
            ev.origins.append(org)
            ev.magnitudes.append(mag)

        return(ev)


                    
            
