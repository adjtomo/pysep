import obspy
import matplotlib.backends

ddir = '/home/vipul/REPOSITORIES/GEOTOOLS/python_util/util_data_syn/'
eid = '20090407201255351'
ftag = '*.z'
sacfiles = ddir + eid + '/' + ftag
streams_per_page = 20

sacstream = obspy.read(sacfiles)
# To plot a record section the ObsPy header trace.stats.distance must be defined in meters (Default)
for tr in sacstream.traces:
    tr.stats.distance = tr.stats.sac['dist']*1000.0

# https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.plot.html#obspy.core.stream.Stream.plot
sacstream.plot(type='section',orientation='horizontal', alpha = 1)
