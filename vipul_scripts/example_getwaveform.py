import obspy
from obspy.clients.fdsn import Client

c = Client("IRIS")

net = 'II'
stn = 'KDAK'
net = 'AK'
stn = 'PPLA'
loc = '*'
chan = 'BH*'
tbefore = 100
tafter = 600 # for extracting the waveform
sec_before_after_event = 10   # for looking up the event
otime = obspy.UTCDateTime("2016-01-24T10:30:29.557")
tstart = otime - tbefore
tend = otime + tafter

# https://docs.obspy.org/packages/autogen/obspy.clients.fdsn.client.Client.get_waveforms.html#obspy.clients.fdsn.client.Client.get_waveforms
# https://docs.obspy.org/master/packages/autogen/obspy.clients.fdsn.client.Client.get_waveforms_bulk.html
st = c.get_waveforms(network = net, station = stn, location = loc, channel = chan, starttime = tstart, endtime = tend)
# check the headers
for tr in st.traces:
    print(tr.stats)

# Get stations (inventory)
# http://docs.obspy.org/archive/0.9.2/packages/autogen/obspy.fdsn.client.Client.get_stations.html
# http://docs.obspy.org/archive/0.10.2/packages/obspy.station.html
stations = c.get_stations(network=net, starttime=tstart, endtime=tend, station=stn, channel=chan, level="response")
print(stations)
# print(stations.networks[0].stations[0].channels[1])

# Get the event from the catalog
cat = c.get_events(starttime = otime - sec_before_after_event, endtime = otime + sec_before_after_event)

