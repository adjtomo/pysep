function [w,respfile_all] = wupdate_sensor(w)
%WUPDATE_SENSOR read sensor information, update fields in waveform object
%
% For sac files produced by run_getwaveform.py, we stick the sensor
% information into several open sac fields. Ideally we would just read the
% KINST sac header, but it only holds 8 characters, which is insufficient
% to describe most sensors.
%
% The spaces on either side of the sac header do not seem to be saved. This
% results in missing spaces in the concatenated strings. I was unable to
% resolve this.
%
% EXAMPLES
%    Episensor 200 Hz 5 Volt per g/Quanterra 330 Linear
%    Guralp 5T/Quanterra 330 Linear Phase Below 100 Com
%
% called by loadsac_all.m
% Carl Tape, 2016/11/11
%

n = length(w);

respfile_all = repmat(cellstr(''),n,1);

for ii=1:n
    kt3 = []; kt4 = []; kt5 = []; kt6 = []; kt7 = []; kt8 = [];
    try
        kt3 = get(w(ii),'KT3');
        if length(kt3)==7, kt3x = [kt3 ' ']; else kt3x = kt3; end
    end
    try, kt4 = get(w(ii),'KT4'); kt4x = addspace(kt4,kt3); end
    try, kt5 = get(w(ii),'KT5'); kt5x = addspace(kt5,kt4); end
    try, kt6 = get(w(ii),'KT6'); kt6x = addspace(kt6,kt5); end
    try, kt7 = get(w(ii),'KT7'); kt7x = addspace(kt7,kt6); end
    try, kt8 = get(w(ii),'KT8'); kt8x = addspace(kt8,kt7); end
    %whos kt3 kt4 kt5 kt6 kt7 kt8
    %whos kt3x kt4x kt5x kt6x kt7x kt8x
    %kt3,kt4,kt5,kt6,kt7,kt8
    %kt3x,kt4x,kt5x,kt6x,kt7x,kt8x
    ka = [kt3 kt4 kt5 kt6 kt7 kt8];
    kb = [kt3x kt4x kt5x kt6x kt7x kt8x];
    %disp([ii length(kt3) length(kt4) length(kt5) length(kt6) length(kt7) length(kt8)]);
    respfile_all(ii) = cellstr(kb);
    
    if ~isempty(respfile_all{ii})
        w(ii) = addfield(w(ii),'RESPFILE',respfile_all{ii});
    end
end

%==========================================================================

function x2 = addspace(x2,x1)

if length(x2)==8
    return
elseif length(x2)==7
    % no easy way to tell which side to put the space on
    %x2 = [' ' x2 ' '];
    return
elseif length(x2)==6
    x2 = [' ' x2 ' '];
end
        
%==========================================================================
