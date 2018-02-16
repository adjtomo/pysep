function [w,units_all] = wupdate_units(w)
%WUPDATE_UNITS read sensor information
%
% For sac files produced by run_getwaveform.py, we stick the units
% into KUSER0 and any scale factor into SCALE.
%
% called by loadsac_all.m
% Carl Tape, 2016/11/11
%

n = length(w);

units_all = repmat(cellstr(''),n,1);

for ii=1:n
    kuser0 = []; fscale = [];
    try, kuser0 = cellstr(get(w(ii),'KUSER0')); end
    try, fscale = get(w(ii),'SCALE'); end

    if ~isempty(kuser0)
        if strcmp(kuser0,'M/S')
            units_all{ii} = 'm/s';
            % update units based on scale factor
            if ~isempty(fscale)
               switch fscale
                   case 1,      units_all{ii} = 'm/s';
                   case 1e2,    units_all{ii} = 'cm/s';
                   case 1e4,    units_all{ii} = 'mm/s';
                   case 1e9,    units_all{ii} = 'nm/s';  
                   otherwise
                      error('check KUSER1 (10^N, N = 0, 2, 4, 9, allowed)');
               end
            end
        else
            if ~strcmp(kuser0,'COUNTS')
                warning(sprintf('input units on KUSER0 are %s, not M/S or COUNTS',char(kuser0)));
            end
        end
    end

    if ~isempty(units_all{ii})
        w(ii) = set(w(ii),'units',units_all{ii});
    end
end
