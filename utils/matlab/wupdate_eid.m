function w = wupdate_eid(w)
%WUPDATE_EID update event ID (since sac header KEVNM is not large enough)
%
% This script should be removed once we fix the issue in run_getwaveform.py
%
% called by loadsac_all.m
% Carl Tape, 2016/11/11
%

warning('adding a missing digit to the event id (KEVNM)');
warning('this assumes your event id is baded on origin time');

n = length(w);

[eid1,msec] = getm(w,'KEVNM','NZMSEC');

eid2 = eid1;
for ii=1:n
   eid2(ii) = {sprintf('%s%s',eid1{ii},num2str(mod(msec(ii),10)))};
   w(ii) = set(w(ii),'KEVNM',eid2{ii});
end
