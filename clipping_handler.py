import os
import obspy



def remove_clipped(path):
    ddir = path+'/'
    rawdir = os.path.join(ddir, 'RAW/')
    # Load all sac files as one stream
    data_dict = obspy.read(os.path.join(rawdir+'*.sac'))\
        ._groupby('{network}.{station}.{location}')
    # Split stream into individual station streams
    data = [data_dict[key] for key in data_dict.keys()]
    # Define clipping factor for 24 bits signal
    q = ((2**(24-1))**2)**0.5

    # Empty list to contain all the clipped traces ids
    clipped_stations = []
    # Loop in station streams
    for station in data:
        # loop in traces
        for trace in station:
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
