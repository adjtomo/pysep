function [w,fnames] = loadsac_all(idir,ftag)
%LOADSAC_ALL load a set of sac files
%
% INPUT
%   idir    directory containing waveforms
%   ftag    suffix for file names, e.g.
%               ftag = {'*.z','*.r','*.t'};
%               
%
% calls wupdate_sensor.m
%

if ~iscell(ftag), ftag = cellstr(ftag); end
ntag = length(ftag);

xx = 0;
for kk=1:ntag
    filetag = strcat(idir,ftag{kk})
    a = dir(filetag);
    n = length(a);
    if n==0
        filetag
    else
        files = {a.name};
        fnames = cellstr(repmat(' ',n,1));
        % initialize
        if xx==0, w(n,1) = waveform; end
        for ii=1:n
            fnames{ii} = strcat(idir,files{ii});
            xx = xx+1;
            % THERE IS NO CHECK FOR MISSING FILES -- IF fnames INCLUDES A
            % FILE NAME THAT DOES NOT EXIST, THIS WILL CRASH WITHOUT ERROR
            % MESSAGE.
            w(xx) = loadsac(waveform,fnames{ii});
        end
    end
end

% create RESPFILE field from headers on sac files generated from run_getwaveform.py
w = wupdate_sensor(w);
% fill units field based on KUSER0 and KUSER1
w = wupdate_units(w);
% add last digit for event id based on origin time
w = wupdate_eid(w);

%==========================================================================
% EXAMPLES

if 0==1
    % load sac files into Matlab
    bdir = '/home/carltape/REPOSITORIES/GEOTOOLS/python_util/util_data_syn/';
    %idir = strcat(bdir,'20160124103030230_iris/');
    idir = strcat(bdir,'20160124103029557_aec/');
    ftag = '*.sac';
    [w,fnames] = loadsac_all(idir,ftag);
    
    % plot record section
    rssort = 2;
    iabs = 0;
    tshift = 100;   % time before origin time
    tmark = [];
    T1 = [];
    T2 = [];
    pmax = 30;
    iintp = [];
    inorm = 1;
    tlims = [];
    nfac = [];
    azstart = [];
    iunit = [];
    imap = 1;
    plotw_rs(w,rssort,iabs,tshift,tmark,T1,T2,pmax,iintp,inorm,tlims,nfac,azstart,iunit,imap);

    %------------------------------------------------
    %% iex=30 in run_getwaveform.py (MFFZ)
    bdir = '/home/carltape/REPOSITORIES/pyseis/';
    %idir = strcat(bdir,'20160124103037400_iex21/');
    %idir = strcat(bdir,'20160124103037400_iex22/');
    %idir = strcat(bdir,'20161106092910579_iex30/');
    %idir = strcat(bdir,'20160114190410727_iex31/');
    idir = strcat(bdir,'20151031025635572/');        % step response
    idir = strcat(bdir,'20171210122855089/'); oshift = 50;
    %idir = strcat(bdir,'20171108064911000/ENZ/');
    ftag = {'*.z','*.r','*.t'};
    %ftag = {'*.z','*.e','*.n'};
    %oshift = 200;
    %ftag = {'*.sac'};
    [w,fnames] = loadsac_all(idir,ftag);
    %w3 = w2tc_rot(w,0);
    
    % use default plotting settings (plus tshift, pmax, inorm, imap)
    plotw_rs(w,[],[],oshift,[],[],[],30,[],1,[],[],[],[],1);
    
    %% low-pass filter to identify corrupted records
    T1 = 1; T2 = 4;
    plotw_rs(w,[],[],oshift,[],T1,T2,30,[],1,[],[],[],[],1);  % normalization
    %%plotw_rs(w,[],[],oshift,[],T1,T2,30,[],0,[],[],[],[],0);  % no normalization
    
    %% COLA and KDAK (iex=22)
    nfac = 2; tlims = [0 200]; iintp = 0;
    plotw_rs(wkeep(w,{'KDAK'}),[],[],oshift,[],T1,T2,30,iintp,1,tlims,nfac,[],[],0);
    nfac = 2; tlims = [50 400]; iintp = 0;
    plotw_rs(wkeep(w,{'COLA'}),[],[],oshift,[],T1,T2,30,iintp,1,tlims,nfac,[],[],0);
    
    % list sensor types
    for ii=1:length(w)/3
        [s1,s2,s3,s4,s5] = getm(w(ii),'network','station','location','channel','RESPFILE');
        disp(sprintf('%2s %4s %2s %3s %s',s1,s2,s3,s4,s5));
    end
    
end

%==========================================================================
