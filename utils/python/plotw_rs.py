
import numpy as np
import numpy.matlib
from math import radians, degrees, sin, cos, asin, acos, sqrt
from geopy import distance
from matlab2datetime import matlab2datetime
import warnings
from numpy import arctan,arctan2,random,sin,cos,degrees
import math
import statistics
from scipy import signal
from obspy.geodetics import kilometers2degrees
import matplotlib.pyplot as plt
from matplotlib import dates
from obspy import read, Stream
from obspy.core import read, UTCDateTime
from datetime import datetime
from obspy.geodetics import gps2dist_azimuth
from scipy.signal import butter, filtfilt, sosfiltfilt


def plotw_rs(win,elat=[], elon=[], rssort=2, iabs=0, tshift=[], tmark=[], T1=[], T2=[], pmax=50, iintp=0, inorm=[1], tlims=[], nfac=1, azstart=[], iunit=1, imap=1, wsyn=[], bplotrs=True, displayfigs='on'):
    '''
    %PLOTW_RS processes waveform object and plots as a record section
    %
    % INPUT
    %   w       waveform object (WHAT ADDED FIELDS SHOULD WE REQUIRE?)
    % INPUT (OPTIONAL) -- these can be omitted or set as [] to get default settings
    %   rssort  =0 input sorting of records
    %           =1 azimuthal sorting of records
    %           =2 distance sorting of records
    %           =3 alphabetical sorting of records (by station name or event ID)
    %   iabs    =1 to plot by absolute distance or azimuth (this means unequal
    %              vertical spacing between the waveforms)
    %   tshift  time shift (in seconds) to apply to the waveforms ([] for none)
    %   tmark   absolute time markers (Matlab serial date) ([] for none)
    %           note: if tshift varies, then you can't use this option
    %   Tfilt   =[Tmin Inf] for high-pass filter
    %           =[0 Tmax] for low-pass filter
    %           =[Tmin Tmax] for band-pass filter
    %   T1      minimum period of filter: =[] for low-pass or no filter
    %   T2      maximum period of bandpass: =[] for high-pass or no filter
    %   pmax    maximum number of time series per record section
    %   iintp   =1 to integrate, =-1 to differentiate, =0 for nothing
    %   inorm   =0 for no amplitude scaling
    %           =1 for scaling by max(abs(d(t))) per trace
    %           =2 for no amplitude scaling except correcting for geometric spreading
    %              inorm can have up to 4 entries:
    %                   inorm(=2)
    %                   GEOFAC
    %                   plot_geometric_spreading
    %                   K
    %   tlims   time limits for x-axis ([] to use default)
    %           note that the 'tshift' field will shift each seismogram
    %   nfac    controls spacing between seismograms
    %   azstart   azimuthal angle for top record (applies with rssort=1 only)
    %   iunit   units for distance record section (applies with rssort=2 only)
    %           =1 for km sphere, =2 for deg sphere, =3 for km flat
    %   imap    plot a simple map showing the station and source locations
    %   wsyn    second waveform object to superimpose on w (e.g., synthetic seismograms)
    %   bplotrs =true to plot record section (can have two additional arguments odir and ftag -- see below)
    %           =false to not plot record section (but presumably to return processed waveforms)
    %
    % OUTPUT (OPTIONAL)
    %   ivec  index in which w(i) are ordered in the plot;
    %            w(ivec) gives the plotting order from top to bottom (if iabs=0)
    %   w     processed waveform (typically data)
    %   wsyn  processed waveform (typically synthetics)
    %   K     2nd parameter for geometric spreading
    %   fH	  figure handles for plots
    %   
    %
    % DETAILS ABOUT THE TIME AXIS
    %   The following variables are pertinent to the time axis:
    %   tshift, tmark, and tlims. tmark is intended to mark absolute times with
    %   vertical bars. If tshift varies from station to station, then the time
    %   axis is relative time, and tmark cannot be used. If there is no tshift
    %   specified, then t=0 will be the earliest time of all the seismograms,
    %   and tlims will be specified w.r.t. this t=0.
    %
    % See examples in run_getwaveform_short.m
    %
    % To print figures to files, set bprint_record_section=true
    %
    % FUTURE WORK:
    %   - if padding zeros, it should be done AFTER demean, taper, etc
    %   - allow pmax to represent the first X sorted seismograms, while not
    %     plotting the rest (which may have too-high SNR, for example)
    %   - combine station and network code to allow for a station to have two
    %     different network codes (e.g., MCK)
    %
    % See GISMO plotting in GISMO/@waveform/plot.m
    %
    % Carl Tape 11/21/2011
    % Yun Wang 11/2011
    %==========================================================================
    '''
    start = datetime.now()
    print('--> entering plotw_rs.m')
    
    narg0 = 18;         # number of input arguments
    spdy = 86400;       # seconds per day
    synplot = 'r';      # plotting a second set of seismograms
    deg = 180/np.pi;
    GEOFAC = 0.5;       # default geometric spreading factor
                        # =0.5 for surface waves (between 0.5 and 1.0 for regional surface waves)
                        # note: GEOFAC = inorm(2)
    bplot_geometric_speading = True;
    T1 = np.atleast_1d(T1)
    T2 = np.atleast_1d(T2)
    tshift = np.atleast_1d(tshift)
    tmark = np.atleast_1d(tmark)
    inorm = np.atleast_1d(inorm)
    tlims = np.atleast_1d(tlims)
    azstart=np.atleast_1d(azstart)
    # options for printing record sections (see also bplotrs)
    bprint_record_section = False;
    odir = './';
    otag = '';
    bprint_map = False;
    fhct = 1; # initialize number of figure handle counts
    #--------------------
    # check input arguments
    w=win.copy()
    print(len(w))
    #w.merge(method=1, fill_value=0)#,interpolation_samples=0)  ## merge any traces with duplicate sta/chans
    if len(w)==0:
        print('empty w') 
        return
    wtemp=Stream()
    for tr in w:
        if statistics.mean(tr.data)==0:
                print('nan trace')
        else:
            wtemp.append(tr)
    w=wtemp
    #print('%i/%i input variables:' % (nargin,narg0))
    # note: the variable will not be listed if it is not present
    #whos w rssort iabs tshift tmark T1 T2 pmax iintp inorm tlims nfac azstart iunit imap wsyn bplotrs
    
    if len(wsyn)==0:
        isyn=0
        wsyn=None
    else:
        print('second waveform object detected -- waveforms will be superimposed');
        isyn=1
    geoinorm = inorm
    # exit here if user enters impermissible values
    if rssort not in [0, 1, 2, 3]: 
        raise ValueError('input rssort = %f must be 0, 1, 2, 3' % (rssort))
    if iabs not in [0, 1]:
        raise ValueError('input iabs = %f must be 0 or 1' % (iabs))
    if len(inorm) >= 2:
        if len(inorm)==3: 
            bplot_geometric_speading = inorm[2]
        GEOFAC = inorm[1]
        inorm = inorm[0]
        print('seismogram normalization:')
        print('  inorm = %s' % (inorm))
        print('  bplot_geometric_speading = %s' % (bplot_geometric_speading))
        print('  GEOFAC = %.2f' % (GEOFAC))
    
    if iintp not in [-1, 0, 1]:
        raise ValueError('input iintp = %f must be -1 or 0 or 1' % (iintp))
    if bplotrs==False:
        raise ValueError('check input: you say you do not want a plot and do not want to return any variables')
    
    # if vertical axis is absolute (distance or azimuth), then plot all
    # waveforms on the same record section (pmax huge)
    if iabs==1:
        pmax = 1000
        if rssort not in [1, 2]:
            print(iabs, rssort)
            raise ValueError('if iabs=1, then rssort=1 (az) or rssort=2 (dist)')
            
    nw=len(w)
    
    ncomp = 1   # w could be nw x ncomp
    # warning: a 1 x 3 w could still all have the same component
    if ncomp not in [1,2,3]:
        #w = w[:]
        #w = np.concatenate(w)
        nw,ncomp = np.shape(w)[0],np.shape(w)[1];
    
    #w = np.concatenate(w)               # convert w to vector
    #tshift = tshift[:]
    ##if len(tshift) != 0:
        ##tshift = np.concatenate(tshift)
    nseis = len(w)
    if pmax==0:
        pmax=nseis
   
    # time shifts (allows for alignment of seismograms on a time that varies
    # from one trace to the next, say, a P arrival)
    if len(tshift) == 0:
        print('no time shift applied (default)')
        tshift = np.zeros((nseis,1))
        
    elif len(tshift)==1:
        print('time shift of %.2f s applied to all waveforms' % (tshift[0]))
        tshift = tshift*np.ones((nseis,1))
        
    elif len(tshift)==nw:
        print('input tshift has dimension %i x %i' % (np.shape(tshift)))
        tshift = tshift[:]
        tshift = np.matlib.repmat(tshift,1,ncomp);
        print('output tshift has dimension %i x %i' % (np.shape(tshift)));
        
    else:
        if nw != 1:
            raise ValueError('tshift (%i) must be same length as w (%i)' % (len(tshift),nseis))
                             
    # time markers (absolute times)
    if len(tmark) == 0:
        nmark = 0
        print('no time markers')
    else:
        if len(np.unique(tshift)) > 1:
            nmark = 0
            print('variable time shifts, so there can be no absolute time markers')
        else:
            nmark = len(tmark)
            print('%i time markers to plot' % (nmark))
        
    # relative time or absolute time
    if len(np.unique(tshift)) > 1:
        itrel = 1
    else:
        tshift0 = tshift[0]
        itrel = 0
    
    #-------------------

    starttime=[]
    endtime=[]
    netwk=[]
    chans=[]
    sta=[]
    rlat=[]
    rlon=[]
    edep=[]
    eid=[]
    mag=[]
    loc=[]
    tdata = []
    trtimes=[]
    evidst=[]
    print(len(w))
    for i, tr in enumerate(w):
        chans.append(tr.stats.channel)
        try:
            rlat.append(tr.stats.sac.stla)
            rlon.append(tr.stats.sac.stlo)
        except:
            rlat.append(tr.stats.coordinates[0])
            rlon.append(tr.stats.coordinates[1])
        starttime.append(tr.stats.starttime)
        endtime.append(tr.stats.endtime)
        evid=str(tr.stats.starttime).replace(":",'')
        evidst.append(evid.replace("-",''))
        netwk.append(tr.stats.network)
        sta.append(tr.stats.station)
        loc.append(tr.stats.location)
        edep.append('dep ')
        eid.append(evidst[0])
        mag.append('NaN')
        elat= elat*np.ones((len(w),1))
        elon= elon*np.ones((len(w),1))
        tdata.append(tr.data)
        trtimes.append(tr.times("timestamp"))
        
    if nseis ==1:
        stchan = chans
    else:
        #unique channels
        uchan = np.unique(chans)
        nuchan = len(uchan)
        stchan = []
        ii=0
        for ii in range(nuchan):
            stchan.append(uchan[ii])
    # FUTURE WORK
    # NEED TO EXIT IF ANY OF THE ABOVE FIELDS ARE EMPTY (isempty does not work)
    
    # we assume that there is either one event or one station in the set of waveforms
    
    if nseis==1:
        nsta = 1
        neve = 1
        irs = 1
        eid = [str(ee) for ee in eid]
        ii=0
        while ii< len(netwk):
            slabs.append(str(netwk[ii]) + '.' + str(sta[ii]) + '.' + str(loc[ii]))
            ii +=1
        sta = [str(ss) for ss in sta];

    else:
        #nsta = length(sta);     % temporary (MCK with two networks)
        nsta = len(np.unique(sta))
        neve = len(np.unique(eid))
        #[nsta,~] = size(unique([rlon(:) rlat(:)],'rows'));
        #[neve,~] = size(unique([elon(:) elat(:)],'rows'));
        if neve==1 and nsta>0:
            irs=1;  # 1 event, multiple stations
            #nsta = nseis;
            #slabs = sta;
            #slabs = np.matlib.repmat(str(''),nseis,1)
            slabs=np.empty(((nseis),1), dtype=object)
            ii=0
            while ii < nseis:
                if nuchan > 1:
                    slabs[ii]=(str(netwk[ii]) + '.' + str(sta[ii]) + '.' + str(loc[ii]) + '.' + str(chans[ii]))
                else:
                    slabs[ii]=(str(netwk[ii]) + '.' + str(sta[ii]) + '.' + str(loc[ii]))
                ii+=1
        elif nsta==1 and neve>0:
            irs=0;  # 1 station, multiple events
            slabs = eid
            #neve = nseis;  
        else:
            # consider TCOL and COLA -- we might want to compare these in a
            # record section for the same event, so here we check the station
            # locations and do NOT exit with an error if the stations are close
            # to each other
            
            irs=0  # 1 station, multiple events
            slabs = eid
            
            dmax = 1e6     # initialize to large number
            xdists=[]
            i=0
            if nsta > 1:
                #xdists = distance.distance(rlat[0]*np.ones(np.shape(rlon)), rlon[0]*np.ones(np.shape(rlon)), rlat, rlon).km
                while i < len(rlat):
                    xdists.append(distance.distance((rlat[0], rlon[0]), (rlat[i], rlon[i])).km)
                    i+=1
                dmax = max(xdists)
    
            if dmax > 0.5:
                print('nsta = %i, neve = %i -- error with waveform data' % (nsta,neve))
                print('check headers KEVNM and station:')
                print(eid)
                print(sta)
                print('check headers KEVNM and station')
                print('must have a single event or station common to all input waveforms')
                
    print('%i event, %i station' % (neve,nsta))
    # reference start time
    tstartmin = min(starttime)
    imin = starttime.index(tstartmin)      
    print('minimum start time of all waveforms is %s (%s)' % (tstartmin,sta[imin]))
    tref = dates.date2num(tstartmin)*spdy
    if itrel==1 and len(tmark)==1:
        tref = dates.date2num(tmark[0])*spdy
    
    stref = print('reference time is %s' % (dates.num2date(tref/spdy)));
    #print(stref);
    print('--> this will be subtracted from all time vectors')
    
    # compute distances and azimuths to stations (or events)
    if irs==1:
        lat1 = elat
        lon1 = elon
        lat2 = rlat
        lon2 = rlon
    else:
        lat1 = rlat[0]
        lon1 = rlon[0]
        lat2 = elat
        lon2 = elon

    
    ###[dist] = distance(lat1,lon1,lat2,lon2, 'degrees')
    def get_bearing(lat1, long1, lat2, long2):
        dLon = (long2 - long1)
        x = math.cos(math.radians(lat2)) * math.sin(math.radians(dLon))
        y = math.cos(math.radians(lat1)) * math.sin(math.radians(lat2)) - math.sin(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.cos(math.radians(dLon))
        brng = arctan2(x,y)
        brng = degrees(brng)

        return brng
    
    dist = []
    azi=[]
    i=0
    #while i < len(lat2):
    for i in range(len(lat2)):
        #dist.append(distance.distance((lat1[i], lon1[i]), (lat2[i], lon2[i])).km)
        dist.append((gps2dist_azimuth(lat1[i], lon1[i], lat2[i], lon2[i])[0])/1000)
        #azimuth=get_bearing(lat1[i], lon1[i], lat2[i], lon2[i])
        azimuth=gps2dist_azimuth(lat1[i], lon1[i], lat2[i], lon2[i])[1]
        
        if azimuth < 0:
            azi.append(azimuth + 360)
            #azi.append(get_bearing(lat1[i], lon1[i], lat2[i], lon2[i]))
        else:
            azi.append(azimuth)
        
        #i+=1
    
    dist_deg = dist
    if iunit not in [1,2,3]:
        warnings.warn('iunit not 1,2,3 -- setting iunit = 1');
        iunit = 1
    if iunit ==1:
        sunit = 'km'
        #dist = deg2km(dist)
    if iunit == 2:
        sunit = 'deg';
        for d in range(len(dist)):
            dist[d]=kilometers2degrees(dist[d])
    if iunit == 3:
        # input lon2 and lon1 are assumed to be utmx and utmy,
        # so dist (from above) is over-written
        sunit = 'km';
        i=0
        dist=[]
        while i < len(lat2):
            dist.append( 1e-3 * np.sqrt( (lon2[i]-lon1[i])^2 + (lat2[i]-lat1[i])^2 ))
            i+=1 
     
    dran = max(dist) - min(dist)
    aran = max(azi) - min(azi)
    print('distance range is %.1f - %.1f = %.1f %s' % (max(dist),min(dist),dran,sunit))
    print(' azimuth range is %.1f - %.1f = %.1f' % (max(azi),min(azi),aran))
    
    # sort
    # NEED TO IMPLEMENT A MULTI-SORT FOR THE CASE OF MULTIPLE SEISMOGRAMS AT
    # THE SAME SITE (if the distance or az are exactly the same, then sort
    # based on the net.sta.loc.chan string)
    sortlab = ['input','azimuth','distance','label']
    if rssort == 0:      # default is no sorting
        nsei=nseis+1
        #ivec = [1:nseis].T
        ivec =arange(1,nsei)[:, newaxis]
    if rssort == 1:      # azimuth
        if len(azstart)==0:
            print(len(azstart))
            azstart = 0
            warnings.warn('azstart not specified -- setting azstart = 0');
        
        ###ivec = np.argsort((azi-azstart) % 360);
        ivec = np.argsort(azi)
        #[~,ivec] = np.argsort(azi-azstar)
    if rssort == 2:      # distance
        ivec = np.argsort(dist);
    if rssort == 3:      # alphabetical
        ivec = np.argsort(slabs);
    
    # labels for record section -- UNSORTED
    rlabels=np.empty(((nseis),1), dtype=object)
    if iabs==0:
        jj=0
        for jj in range(nseis):
            # variable time shift
            if len(np.unique(tshift)) > 1:
                if max(tshift) < spdy:
                    if tshift[jj] < 100:
                        stshift = ('DT %.1f s' % (tshift[jj]))
                    else:
                        stshift = ('DT %.1f min' % (tshift[jj]/60))
                    
                else:
                    stshift = ''
                
                #stshift = sprintf('DT %.1f',tshift(jj)-min(tshift));
            else:
                stshift = ''
            
            if irs==1:   # 1 event, multiple stations
                rlabels[jj] = ('%s (%.0f, %.0f %s) %s' % 
                    (str(slabs[jj])[2:-2],math.floor(azi[jj]),dist[jj],sunit,stshift))
            else :       # 1 station, multiple events (list event depths)
                rlabels[jj] = ('%s %s (%.0f, %.0f %s) %.0f km %.1f %s' % 
                    (str(slabs[jj])[2:-2], chans[jj], math.floor(azi[jj]), dist[jj], sunit, edep[jj], mag[jj], stshift))
            
            
        
    else:
        rlabels = slabs
    
    
    #--FILTER WAVEFORMS--------------------------------------------------------
       
    RTAPER = 0.05
    
    # filter
    if len(T1)==0 and len(T2)==0:
        ifilter = 0
        print('NO FILTER WILL BE APPLIED')
        stfilt = '--'
        wtemp=Stream()
        for tr in w:
            fval=statistics.mean(tr.data)
            #tr.trim(min(starttime), max(endtime), pad=True, fill_value=fval)
            tr.trim((tr.stats.starttime),(tr.stats.endtime), pad=True, fill_value=fval)
            wtemp.append(tr)
        w=wtemp
        
    else:
        ifilter = 1;
        npoles = 2;
        
        # fill gaps with mean value
        #w = fillgaps(w,'meanAll');
        wtemp=Stream()
        for tr in w:
            fval=statistics.mean(tr.data)
            tr.trim((tr.stats.starttime),(tr.stats.endtime), pad=True, fill_value=fval)
            #tr.trim(min(starttime), max(endtime), pad=True, fill_value=fval)
            wtemp.append(tr)
        w=wtemp
        print(max(w[0]))
        # these operations might depend on whether the input is displacements
        # (which could have static offsets) or velocities
        print('pre-processing: detrend, demean, taper');
        
        w.detrend('demean')
        '''wtemp=Stream()
        for tr in w:
            fval=statistics.mean(tr.data)
            dmean=tr.data/abs(fval)
            tr.data=dmean
            wtemp.append(tr)
        w=wtemp'''
        #w = demean(w); 
        w.taper(max_percentage=RTAPER, type='cosine');
        wfilt=Stream()
        Tmax_for_mHz = 100     # list mHz, not Hz, for T2 >= Tmax_for_mHz
        if len(T1) != 0 and len(T2) == 0:
            print('%i-pole low-pass T > %.1f s (f < %.2f Hz)' % (npoles,T1[0],1/T1[0]))
            #f = filterobject('L',1/T1,npoles);
            stfilt = ('T > %.1f s (f < %.2f Hz)' % (T1[0],1/T1[0]))
            if T1[0] >= Tmax_for_mHz: 
                stfilt = ('T > %.1f s (f < %.1f mHz)' % (T1[0],1/T1[0]*1e3))
            try:
                wfilt=w.filter("lowpass", freq=1/T1[0],zerophase=True)
            except:
                print("filter didn't work")
                
        elif len(T1) == 0 and len(T2) != 0:
            print('%i-pole high-pass T < %.1f s (f > %.2f Hz)' % (npoles,T2[0],1/T2[0]))
            #f = filterobject('H',1/T2,npoles);        
            stfilt = ('T < %.1f s (f > %.2f Hz)' % (T2[0],1/T2[0]))
            if T2[0] >= Tmax_for_mHz:
                stfilt = ('T < %.1f s (f > %.1f mHz)' % (T2[0],1/T2[0]*1e3))
            try:
                wfilt=w.filter("highpass", freq=1/T2[0],zerophase=True)
            except:
                print("filter didn't work")
            
        elif len(T1) != 0 and len(T2) != 0:
            print('%i-pole band-pass filter between T = %.1f - %.1f s' % (npoles,T1[0],T2[0]))
            #f = filterobject('B',[1/T2 1/T1],npoles);
            stfilt = ('T = %.1f-%.1f s (%.2f-%.2f Hz)' % (T1[0],T2[0],1/T2[0],1/T1[0]))
            if T2[0] >= Tmax_for_mHz:
                stfilt = ('T = %.1f-%.1f s (%.1f-%.1f mHz)' % (T1[0],T2[0],1/T2[0]*1e3,1/T1[0]*1e3))
            try:
                
                for tr in w:
                    #sos = butter(5, [1/T2[0], 1/T1[0]], 'bandpass', output='sos');
                    lowcut=1/T2[0]
                    highcut=1/T1[0]
                    nyq = 0.5 * tr.stats.sampling_rate
                    low = lowcut / nyq
                    high = highcut / nyq
                    sos = butter(5, [low, high], analog=False, btype='band', output='sos')
                    #y = filtfilt(b, a, tr.data, padtype = 'odd', padlen=3*(max(len(b),len(a))-1))
                    y = sosfiltfilt(sos, tr.data.copy())
                    tr.data=y
                    wfilt.append(tr)
                    
                #wfilt=w.filter("bandpass", freqmin=1/T2[0], freqmax=1/T1[0],zerophase=True)
            except:
                print("filter didn't work")
        w=wfilt
        
        '''wtemp=Stream()
        for tr in w:
            fval=statistics.mean(tr.data)
            dmean=tr.data/abs(fval)
            tr.data=dmean
            wtemp.append(tr)
        w=wtemp'''
        
        # apply identical filtering to synthetics, if present
        
        if isyn==1:
            wsyn.detrend()
            #wsyn = demean(wsyn); 
            wsyn.taper(max_percentage=RTAPER)
            wsyn.filter("bandpass", freqmin=1/T2[0], freqmax=1/T1[0],corners=npoles,zerophase=True)
    ##for trr in w: 
        ##trr.data=trr.data/max(abs(trr.data))
    # integrate or differentiate
    # note: units of w will automatically change
    units = 'nm / sec'
    if iintp==1:
        if ifilter==0: 
            w.detrend()
        w=w.integrate()
        units = 'nm'
        if isyn==1:
            if ifilter==0:
                wsyn.detrend()
            wsyn.integrate()
        
    if iintp==-1:
        w.differentiate()        # help waveform/diff
        units = 'nm / sec^2'
        if isyn==1:
            wsyn.differentiate()
    
    # compute amplitude scaling for plots
    # NOTE: This will normalize by the full time series, not simply by the time
    #       specified within the window denoted by tlims.
    nlabs = ['none','max(abs(d_i))',('(sin D)^-%.2f' % (GEOFAC))];
    nlab = 'norm --> ' + str(nlabs[inorm[0]]);
    nvec = np.ones((nseis,1))                          # normalization vector for seismograms
    #w.merge(fill_value=0)
    wmaxvec = []
    #w.normalize(global_max=True)
    ftrs=[]
    fttimes=[]
    for tra in w:
        #print(max(abs(tra.data)))
        wmaxvec.append(max(abs(tra.data[100:-100])))
        ftrs.append(tra.data)
        fttimes.append(tra.times("timestamp"))
        #tra.plot()
    
    #wmaxvec = max(abs(w.merge(fill_value=0)))          # will handle empty records
    #print(wmaxvec)
    '''
    ### come back to deal with synthetic arguments ###
    if isyn==1
        wsynmaxvec = max(abs(fillgaps(wsyn,0)));
        disp('amplitude comparison between data and synthetics:');
        for ii=1:length(w)
            disp(sprintf('%5s %s ln(data/syn) = %6.2f',...
                sta{ii},chans{ii},log(wmaxvec(ii) / wsynmaxvec(ii))));
        end
        wmaxvec = max([wmaxvec' ; wsynmaxvec'])';
    end
    '''
    K = []
    fH = []
    if geoinorm[0] == 2:
    #if inorm(1)==2
        # For information on geometric spreading, see Stein and Wysession,
        # Section 4.3.4, Eq 20 (which is for Rayleigh waves. GEOFAC = 0.5).
        # For our purposes, this will fold the attenuation factor into the same
        # single factor that accounts for geometric spreading.
        # WARNING: This does not take into account the variation in amplitude. 
        print('correct for geometrical spreading using (sin x)^-%.2f' % (GEOFAC))
        sindel = np.sin(np.asarray(dist_deg) / deg)
        # note: stations that are father away get enhanced by a larger valuye
        Kvec = wmaxvec * (sindel**GEOFAC)
        # single parameter for fitting (GEOFAC fixed)
        if len(geoinorm)==4:
            K = geoinorm[3]               # user-specified K value
        else:
            K = statistics.median(Kvec)
 
        wvec = K/(sindel**GEOFAC)     # factor will depend on source-station distance
                                        # but not on source depth
        ii=0
        while ii < len(w):
            # THIS CHANGES THE AMPLITUDE OF w(ii)
            w[ii].data = w[ii].data / wvec[ii] 
            if isyn==1:
                wsyn[ii] = wsyn[ii] / wvec[ii] 
            print('%5s %8.3f deg, max = %8.2e, maxG = %5.2f, K/sin(del)^x = %16.2f' %
               (sta[ii],dist_deg[ii],wmaxvec[ii],max(abs(w[ii].data)),wvec[ii]))
            ii+=1
        # extra figure with best-fitting curve
        if bplot_geometric_speading:
            fsizeg = 14 
            msizeg = 20
            if len(w) > 10: 
                fsizeg = 10 
                msizeg = 12
            
            geofig=plt.figure()
            plt.plot(dist_deg,wmaxvec,'ko',markersize=msizeg,markerfacecolor='r')
            plt.text(dist_deg,wmaxvec,sta,fontsize=fsizeg)
            # best-fitting curve
            x = np.linspace(0,min([1.1*max(dist_deg), 180]))
            plt.plot(x,K/(sin(x/deg)**GEOFAC),'r--',linewidth=2)
            plt.ylim(0, 1.05*max(wmaxvec))
            plt.xlabel('\Delta, source-station arc distance, deg',fontsize=14);
            plt.ylabel('max( |v(t)| )',fontsize=14);
            plt.title('A(\\Delta) = %.2e / (sin \\Delta)^{%.2f}' % (K,GEOFAC),fontsize=16)
        
        # update wmaxvec
        #wmaxvec = max(abs(fillgaps(w,0)))
        '''
        if isyn==1:
            wsynmaxvec = max(abs(fillgaps(wsyn,0)))
            wmaxvec = max([wmaxvec' ; wsynmaxvec'])'
            '''
       
    wmax = max(wmaxvec)    
    
    #--PLOT RECORD SECTIONS----------------------------------------------------
    
    if bplotrs == False:
        return
    
    # get units, which may have changed
    units = 'nm / sec'
    if nseis==1:
        stunit = units;
    else:
        # unique units
        #uunit = np.unique(get(w,'units'));
        uunit = np.unique(units)
        nuunit = len(uunit)
        if nuunit > 1:
            warnings.warn('multiple units are present on waveforms')
        stunit = []
        ii=0
        while ii < nuunit: 
            stunit.append(uunit[ii])
            ii +=1
    
    nfig = np.ceil(nseis/pmax);
    fsize = 10
    kk = 0
    az1 = azstart
    azinc = 45
    #azbin = np.arange(azstart , azstart+360, azinc) % 360
    #print(azbin)
    try:
        azbin = np.arange(azstart , azstart+360, azinc) # %360
    except:
        azbin=[]
    azbin_fwid = 0.1    
    for a in range(len(azbin)):
        azbin[a]=azbin[a]%360
    # vertical separation between seismograms
    wsep = 1.5*statistics.mean(wmaxvec)
    yshift = wsep*nfac
    
    if inorm[0]==1:
        #nvec = max(abs(w)); 
        nvec = wmaxvec
        yshift = nfac
    
    if iabs==1:
        nvec=np.ones((len(wmaxvec),1))
        if inorm[0]==1:
            if rssort==1:
                yran=aran
            else:
                yran=dran
                if iunit==2:
                    yran = kilometers2degrees(dran)
            
            for ii in range(len(wmaxvec)):
                nvec[ii] = nfac * nseis/yran*wmaxvec[ii]
            
        else:        # inorm = 0 or 2
            wvec = wmax*np.ones((nseis,1))
            for jj in range(len(wvec)):
                nvec[jj] = nfac*wvec[jj]
        
        # sign flip is needed, since y-axis direction is flipped
        nvec = -nvec
    
    print('wmax = %.3e, wsep = %.3e, yshift = %.3e' % (wmax,wsep,yshift))
    
    # debugging for variable nvec:
    #disp(sprintf('summary of nvec (%i): min/mean/median/max = %.3e / %.3e / %.3e / %.3e',...
    #    length(nvec),min(abs(nvec)),mean(abs(nvec)),median(abs(nvec)),max(abs(nvec))));
    
    # get time limits for record section
    # NOTE: THIS SHOULD BE AVOIDED SINCE IT INVOLVES READING ALL THE WAVEFORMS
    if len(tlims) == 0:
        #tic
        tlim1 = np.zeros((nw,1))
        tlim2 = np.zeros((nw,1))
        ii=0
        while ii < nw:
            #ti3 = get(w(ii),'timevector')
            ti3 = trtimes[ii]
            # KEY: time vector for plotting
            tplot = (ti3 - tref) - tshift[ii];
            tlim1[ii] = min(tplot);
            tlim2[ii] = max(tplot);
            ii+=1
        #print('it took %.2f s to establish the time limits for plotting' % (toc))
        tlims = [min(tlim1), max(tlim2)]
        print('tlims',tlims)
    pp=0
    jj=0
    kk=0
    
    while pp < nfig:
        print('record section page %i/%i (max %i per page)' % (pp,nfig,pmax))
        
        #tempfh = figure('Visible',displayfigs); hold on;
        figname='fig'+str(pp)
        figname=plt.figure(figsize=(9,11))
        ax = plt.subplot(111)

        #fig=plt.figure()
    
        if nfig==1:
            jmax = nseis;
        else:
            jmax = pmax;
            
        # initialize arrays
        dy = np.zeros((jmax,1))
        rlabs=np.empty(((jmax),1), dtype=object)
        dplotmax = 0;
        
        wtemp=Stream()
        jj=0
        
        for jj in range(jmax):       # loop over seismograms
            if kk<nseis:
            
                ii = ivec[kk];  # key sorting
                rlabs[jj] = rlabels[ii]
                # get seismogram
                ti3 = fttimes[ii]
                di3 = ftrs[ii]
                
                # KEY: time vector for plotting
                tplot = (ti3 - tref) - tshift[ii]
                if iabs==0:
                    # plot from top to bottom
                    dy[jj] = (jmax + 1 - jj)*yshift
                else:
                    # use absolute scaling (e.g., plot records at their actual distance)
                    if rssort==1: 
                        dy[jj] = azi[ii] 
                    else: 
                        dy[jj] = dist[ii]
                #norm=np.linalg.norm(di3)
                
                dplot = di3/nvec[ii]      # key amplitude scaling
                dplotshift = dplot + dy[jj]
                
                # get the max value of the seismogram within the plotted time interval
                #btplot = and(tplot > tlims(1),tplot < tlims(2));
                #dplotmax = max([dplotmax max(abs(dplot(btplot)))]);
                #plt.plot(ti3,dplot,'b')
              
                ax.plot(tplot,dplotshift,'b', linewidth=0.5);
                # PLOT LABELS FOR EACH WAVEFORM
                txtplace=max(rlabels, key=len)
                txtplce=len(str(rlabels[ii])[2:-2])
                
                plt.text(tlims[0], statistics.mean(dplotshift), str(rlabels[ii])[2:-2], ha='right', fontsize=8)
                
                
                # specify the amplitude of the first seismogram plotted
                # (note that this is not the maximum over all seismograms plotted)
                
                if jj==0:
                    imx = np.argmax(abs(di3))
                    mx=di3.max()
                    stmx = ('%s max %.2e %s at t = %.1f s ' % (sta[ii],di3[imx],units,tplot[imx]))
                
                '''
                ### come back to deal with synthetic arguments ###
                if isyn==1
                # be careful about the reference time
                [tisyn,disyn,tstartsyn] = getm(wsyn(ii),'timevector','data','start');
                tplot = (tisyn - tstartsyn)*spdy - tshift(ii);
                dplotshift = disyn/nvec(ii) + dy(jj);   % key amplitude scaling
                plot(tplot,dplotshift,synplot); 
                end
                '''
                
                # partition if plotting by azimuth (see azbin above)
                if rssort==1 and iabs==0:
                    # azimuth of previous (az0) and current (az1) stations in the sorted list
                    az0 = az1
                    az1 = azi[ii]
                    
                    # (this boolean clause could probably be simplified) 
                    #if (az1 > az0 and (any(az1>azbin) and any(azbin>az0))) or(az1 < az0 and (any(az1>azbin)or any(azbin>az0))):
                    if az1 > az0:
                        for ii in range(len(azbin)):
                            if az1>azbin[ii] and azbin[ii]>az0:
                                
                                # note: it would nice if these bars extended to the LEFT,
                                # outside the plotting axes
                                tp1 = tlims[0];
                                tp2 = tlims[0] + azbin_fwid*(tlims[1]-tlims[0]);
                                #disp(sprintf('%i %i %.2f %.2f',pp,jj,tp1,tp2));  % testing
                                if wsyn != None: # and var != None:
                                    pc = 'r' 
                                else: 
                                    pc = 'k'
                                plt.plot([tp1, tp2],[dy[jj]+yshift/2 ,dy[jj]+yshift/2],pc,'linewidth')
                                break
                    elif az1<az0:
                        for ii in range(len(azbin)):
                            if az1>azbin[ii] or azbin[ii]>az0:
                                
                                # note: it would nice if these bars extended to the LEFT,
                                # outside the plotting axes
                                tp1 = tlims[0];
                                tp2 = tlims[0] + azbin_fwid*(tlims[1]-tlims[0]);
                                #disp(sprintf('%i %i %.2f %.2f',pp,jj,tp1,tp2));  % testing
                                if wsyn != None: # and var != None:
                                    pc = 'r' 
                                else: 
                                    pc = 'k'
                                plt.plot([tp1, tp2],[dy[jj]+yshift/2 ,dy[jj]+yshift/2],pc,'linewidth')
                                break
                # exit early for the last page of the multi-page record section
                if pp == nfig and jj == nseis % pmax:
                    print('ending early')
                    if iabs==0:
                        dy = (jmax + 1 - arange(jmax))*yshift
                    return
    
                
                kk = kk+1
            
        if iabs==0:
            
            # this may chop off a large-amplitude trace near the boundary, but
            # at least the gap will be the same for a sequence of plots

            ylims = [0 ,yshift*(jmax+2)]
            plt.yticks([])
            #plt.ylim(ylims)
            # this will provide more space if one of the traces has a large
            # amplitude, but the gap may vary for a sequence of plots
            #ylims = yshift*[1 ,jmax] + dplotmax*[-1 ,1]   
            
        else:
            # using actual values of distance or azimuth
            ax.yaxis.tick_right()
            if rssort==1:
                ylims = [min(azi)-5, max(azi)+5];
                if max(azi)-min(azi) > 300:
                    ylims = [-10, 370] 
                
            else:
                ylims = [min(dist)-0.1*dran, max(dist)+0.1*dran]
            
        print('tlims = %.2f to %.2f' % (tlims[0],tlims[1]))
        print('ylims = %.3e to %.3e (yshift = %.3e, jmax=%i)' % (ylims[0],ylims[1],yshift,jmax))
        
        # plot tshift marker
        #plt.plot([0, 0],ax0[2:3],'r','linewidth')
        plt.vlines(0,ylims[0],ylims[1], 'r')
        # plot absolute-time marker
        mm=0
        for mm in range(nmark):
            tm = (dates.date2num(tmark[mm]) - dates.date2num(tstartmin))*spdy - tshift0;
            #print((dates.date2num(tmark[mm]) - dates.date2num(tstartmin))*spdy)
            #plt.plot(tm*[1, 1],ax0[2:3],'r','linewidth')
            plt.vlines(tm,ylims[0],ylims[1], 'r')
        plt.xlim(tlims)
        plt.ylim(ylims[0],ylims[1])
        if iabs==1:
        #if rssort ==1:
            plt.gca().invert_yaxis()
        #plt.subplots_adjust(left=0.10)
        plt.xlabel('Time (s)', fontweight='bold')
        
        
        # title string explaining the time axis
        if itrel==0:
            
            t1a = dates.date2num(tstartmin) + (tshift0+tlims[0])/spdy
            t2a = t1a + np.diff(tlims)
            
            stdur = ('%s + %.2f s' % (dates.num2date(np.double(t1a)),t1a))
        else:
            # relative time shifts
            stdur = ('variable time shifts: %s' % (stref))
            #  if length(tmark)==1
            #      stdur = sprintf('variable time shifts: %s',stref);
             #else
             #    % just pick the iith record as an example -- this should be the
             #    % bottom time series on the record section
             #    ix = ii;
             #    t1a = tstartmin + (tshift(ix)+tlims(1))/spdy;
             #    t2a = t1a + diff(tlims);
             #    stdur = sprintf('variable time shifts: w(%i) (%s) is %s + %.2f s',...
             #        ix,sta{ix},datestr(double(t1a),31),t2a-t1a);
    
        # plot title
        # note: might want to label the channel, too
        #if ifilter==1, stfilt = sprintf('T = %.1f-%.1f s (%.2f-%.2f Hz)',T1,T2,1/T2,1/T1); else stfilt = '--'; end
        
        st0 = ('%s [%s, %s]' % (stchan,stunit,stfilt))
        if irs==1:
            
            stline1 = ('event %s (%s, M%s, %.1f, %.1f, z = %s km)' %
                (eid[0],starttime[0],mag[0],elon[0],elat[0],edep[0]))
            stline2 = ('%i / %i seismograms (%i stations) ordered by %s, %s' %
               (jmax,nseis,nsta,sortlab[rssort],nlab))
        else:
            stline1 = ('station %s (%.1f, %.1f)' % (sta[0],rlon[0],rlat[0]))
            stline2 = ('%i / %i seismograms (%i events) ordered by %s, %s' %
               (jmax,nseis,neve,sortlab[rssort],nlab))
        
        plt.title(str(stdur) + ';  ' + str(stmx)+ '\n' + str(st0)+ '\n' + str(stline1)+ '\n' + str(stline2), fontsize=8)
        #title({[stdur ';  ' stmx],sprintf('%s %s',st0),stline1,stline2},'fontsize',fsize,'interpreter','none');
                
        
        
        #plt.ion
        plt.tight_layout()
        #plt.show()
        pp+=1
    #return fig
    
    #--------------------------------------------------------------------------
    # plot map
    # future work would be to add some more options for alaska plots (alaska_basemap.m)
    if imap==1:
        if iunit==3: 
            fac = 1e-3    
        else:
            fac = 1
        iplotsrc = 1
        fsize = 10
        figsta=plt.figure(figsize=(8,10)) 
        plt.plot(rlon*fac,rlat*fac,'bv') 
        #if iplotsrc==1:
        plt.plot(elon*fac,elat*fac,'k*',markersize=15,markerfacecolor='r')
        # plot station labels (or eid labels)
        if irs==1:
            for i,lab in enumerate(sta):
                plt.text(rlon[i]*fac,rlat[i]*fac,str(lab),fontsize=fsize)
        else:
            plt.text(elon[i]*fac,elat[i]*fac,eid[i],fontsize=fsize)
        
        plt.title(str(stline1) +'\n'+ str(stchan),fontsize=8)
        plt.show()
    print('--> leaving plotw_rs.m')
    # Print download time
    d = datetime.now() - start
    print (d, " s to plot waveforms")

    #==========================================================================

    return w
    
    
    
    
