function write_seis_obspy(ftag,lon,lat,dep,mag,otime)
%WRITE_SEIS_OBSPY write a text file to load into a python script
%
% This is a (hopefully) temporary script to write a text file containing
% input parameters for a python script to fetch waveforms from the IRIS DMC
% and then process them using obspy.
%

n = length(lon);

filename = [ftag '_obspy.txt'];

disp(['write_seis_obspy.m: writing ' filename]);

stfmt = '%3i  %s  %s %12.6f %12.6f %8.1f %6.2f';

eid = otime2eid(otime);

fid = fopen(filename,'w');
for ii = 1:n
    % note: depth is in meters
    fprintf(fid,[stfmt '\n'],...
        ii,eid{ii},datestr(otime(ii),'yyyy-mm-ddTHH:MM:SS.FFF'),...
        lon(ii),lat(ii),dep(ii)*1000,mag(ii));
end
fclose(fid);

%==========================================================================
% EXAMPLE

if 0==1
    % origin times of target events in the Alaska catalog
    otar = [    datenum('2015/11/20 10:53:48')
                datenum('2015/10/22 13:16:15')
                datenum('2014/12/13 15:47:31')
                datenum('2015/10/31 02:56:35')
                datenum('2015/11/06 01:20:12')
                datenum('2014/10/21 00:36:58')
                datenum('2014/10/23 16:30:23')
                datenum('2016/01/14 19:04:10')
                datenum('2015/09/28 11:50:13')
                datenum('2016/01/24 10:30:29')
                datenum('2016/11/06 09:29:10')
                datenum('2016/12/08 10:18:13')
            ];
    n = length(otar);
    % get the source parameters from the AEC catalog
    eid = repmat(cellstr(''),n,1);
    for ii=1:n
    [otime(ii),lon(ii),lat(ii),dep(ii),mag(ii),eid(ii)] = read_eq_AEC([otar(ii) 1],[],[]);
    end

    ftag = '~/MFFZ';
    write_seis_obspy(ftag,lon,lat,dep,mag,otime);
end

%==========================================================================

