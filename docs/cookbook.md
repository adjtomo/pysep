# Cookbook

Here we provide some code snippets using PySEP routines or utilities to perform
tasks useful for manipulating or preparing data for moment tensor and waveform 
inversion software.

## Generating Source Receiver Maps

PySEP generates QuakeML and StationXML files during data gathering, which can 
be used to plot source receiver maps. By default, these files are saved to the
main output directory under the names `event.xml` and `stations.xml`.

```python
from obspy import read_events, read_inventory
from pysep.utils.plot import plot_source_receiver_map

event = read_events("event.xml")[0]
inv = read_inventory("stations.xml")

plot_source_receiver_map(inv=inv, event=event)
```

Stations can also be subset, e.g., when plotting subsetted record sections. The
`subset` argument takes a list of trace IDs which can be retrieved from each
trace.

```python
st_new = st.select(network="II")  
subset = [tr.get_id() for tr in st_new]

plot_source_receiver_map(inv=inv, event=event, subset=subset)
```

## Create STATIONS file for SPECFEM

SPECFEM requires a STATIONS file that defines geographical coordinates of stations that
will be used in a simulation. Using PySEP and ObsPy, you can quickly generate a STATIONS
file for a given region

```python
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from pysep.utils.io import write_stations_file

c = Client("IRIS")
# This is for the region of northern Alaska, broadband seismometers only
# collecting online stations from Jan. 1, 2000 until present time
inv = c.get_stations(network="*", station="*", 
                     channel="BH?,BL?,HH?,HL?",
                     minlatitude=64.5, maxlatitude=72., 
                     minlongitude=-168., maxlongitude=-140.,
                     level="station", 
                     starttime=UTCDateTime("2000-01-01T00:00:00"), 
                     endtime=UTCDateTime()
                     )

# Write 'STATIONS' file, ordered alphabetically by station name
write_stations_file(inv, fid="STATIONS", order_by="station")
```

```
$ head STATIONS
  A19K    AK     70.2043   -161.0713    0.0    0.0
  A21K    AK     71.3221   -156.6175    0.0    0.0
  A22K    AK     71.0033   -154.9742    0.0    0.0
   ANM    AK     64.5646   -165.3732    0.0    0.0
  B18K    AK     69.3641   -161.8016    0.0    0.0
  B20K    AK     70.0079   -157.1599    0.0    0.0
  B22K    AK     70.3400   -153.4196    0.0    0.0
  C16K    AK     68.2746   -165.3436    0.0    0.0
  C18K    AK     68.6483   -161.1943    0.0    0.0
  C19K    AK     69.1049   -159.5874    0.0    0.0

```

## Fetch moment tensors as CMTSOLUTION files

It's often useful to gather a catalog of events and their moment tensors as a starting point 
for waveform inversion. PySEP has some built in moment tensor routines that can help us 
with that.

### From the GCMT catalog

This routine queries the GCMT catalog for moment tensors matching a given origin time and
magnitude. Given an example event, we can return a Catalog object containing the GCMT 
moment tensor. Taking this [M6.4 in northern Alaska](https://earthquake.alaska.edu/event/018aap2cqu) as an example.

```python
from obspy import UTCDateTime
from pysep.utils.mt import get_gcmt_moment_tensor

origintime = UTCDateTime("2018-08-12T14:58:53")
magnitude=6.4
cat = get_gcmt_moment_tensor(origintime=origintime, magnitude=magnitude)
cat.write("CMTSOLUTION_GCMT", format="CMTSOLUTION")
```

```
$ cat CMTSOLUTION_GCMT
 PDE 2018 08 12 14 58 54.30   69.5600 -145.3000   2.2 6.4 6.4 NORTHERN ALASKA
event name:   C201808121458A
time shift:           7.9000
half duration:        4.0000
latitude:            69.7400
longitude:         -144.7800
depth:               12.0000
Mrr:           -7.690000E+24
Mtt:           -1.990000E+25
Mpp:            2.760000E+25
Mrt:            4.890000E+24
Mrp:           -1.510000E+25
Mtp:           -4.730000E+25
```

### From USGS Catalog

Using the event we grabbed in the previous section (from GCMT), we can query 
the USGS catalog for their moment tensor solution.

>__NOTE__: USGS moment tensor searches are more strict and require an existing 
  Event object from which it takes the hypocentral location, origin time and magnitude

```python
# Continuing from codeblock above
from pysep.utils.mt import get_usgs_moment_tensor

cat = get_usgs_moment_tensor(event=cat[0])
cat.write("CMTSOLUTION_USGS", format="CMTSOLUTION")
UserWarning: Moment tensor has no source time function. Half duration will be set to 1.0.
```

```
$ cat CMTSOLUTION_USGS
 PDE 2018 08 12 14 58 53.50   69.5762 -145.2910  15.8 6.4 6.4 NORTHERN ALASKA    
event name:89 km SW of Kaktovik, Alaska                                          
time shift:          -0.2030                                                     
half duration:        1.0000                                                     
latitude:            69.6864                                                     
longitude:         -144.8187                                                     
depth:               11.5000                                                     
Mrr:           -1.015100E+25                                                     
Mtt:           -1.644300E+25                                                     
Mpp:            2.659300E+25                                                     
Mrt:           -8.487000E+24                                                     
Mrp:           -1.468700E+25                                                     
Mtp:           -4.792400E+25  
```

### From GeoNet Catalog (New Zealand)

For those working with New Zealand data, it is possible to grab 
moment tensor data from [John Ristau's catalog](https://github.com/GeoNet/data/tree/main/moment-tensor). 
For these events you will need to provide a corresponding GeoNet event id. 
We'll take [event 2018p130600](https://www.geonet.org.nz/earthquake/2018p130600) as an example.

```python
from obspy.clients.fdsn import Client
from pysep.utils.mt import get_geonet_mt

c = Client("GEONET")
cat = c.get_events(eventid="2018p130600")
focal_mechanism = get_geonet_mt(event_id="2018p130600", units="nm")  # units can also be dynecm
cat[0].focal_mechanisms = [focal_mechanism]
cat.write("CMTSOLUTION_GEONET", format="CMTSOLUTION")
```
```
$ cat CMTSOLUTION_GEONET
 PDE 2018 02 18 07 43 48.13  -39.9490  176.2995  20.6 5.2 5.2 NORTH ISLAND, NEW ZEALAND
event name:           8E23C1
time shift:           0.0000
half duration:        0.6989
latitude:           -39.9490
longitude:          176.2995
depth:               20.5946
Mrr:           -2.479380E+23
Mtt:            1.314880E+23
Mpp:            1.164500E+23
Mrt:            5.032500E+22
Mrp:            6.607700E+22
Mtp:            9.359300E+22
```
-------------------------

## Convert SPECFEM-generated synthetics to SAC files

Some may find it useful to convert synthetic seismograms into SAC files to simplify later processing. 
The following code snippet will write out SAC files on a per-component basis for synthetics generated 
in a geographic coordinate system (lat/lon) from SPECFEM2D, 3D or 3D_GLOBE.

>__NOTE:__ The `read_sem` function automatically appends SAC headers using the provided metadata

```python
# cd path/to/specfem_workdir
from glob import glob
from pysep.utils.io import read_sem

for fid in glob("./OUTPUT_FILES/*sem*):
    st = read_sem(fid=fid, source="./DATA/CMTSOLUTION", stations="./DATA/STATIONS")
    st.write(os.path.basename(fid), format="SAC")
```

If your synthetics were generated in a Cartesian coordinate system (XYZ) you will need to use a 
separate function, as ObsPy does not like coordinate systems that do not adhere to WGS84

```python
# cd path/to/specfem_workdir
from glob import glob
from pysep.utils.io import read_sem_cartesian

for fid in glob("./OUTPUT_FILES/*sem*):
    st = read_sem_cartesian(fid=fid, source="./DATA/CMTSOLUTION", stations="./DATA/STATIONS")
    st.write(f"{st[0].get_id()}.sac", format="SAC")
```

-------------------------

## Append SAC headers to existing Streams

To append SAC headers to your own seismic data, you can directly use the
`PySEP` utility functions `append_sac_headers` and 
`format_sac_header_w_taup_traveltimes`

```python
from obspy import read, read_events, read_inventory
from pysep.utils.cap_sac import append_sac_headers, format_sac_header_w_taup_traveltimes

st = read()
inv = read_inventory()
event = read_events()[0]
st = append_sac_headers(st=st, inv=inv, event=event)
st = format_sac_header_w_taup_traveltimes(st=st, model="ak135")
```

-------------------------

##  Set custom order on output station file 

PySEP creates a text file for stations gathered. This list contains station IDs, coordinates,
epicentral distance and azimuth from a given event. By default, this list is ordered by the 
internal order of the station metadata (usually alphabetical). 

Users can set the order of this station list by specifying the `order_stations_by` parameter.
This is provided as a keyword argument to the initiation of PySEP and can be set with:

```bash
pysep -c pysep_config.yaml --order_stations_list_by distance
```

Acceptable arguments for this are: 'network', 'station', 'latitude', 'longitude', 'distance', 'azimuth'

You can also call this from a Python environment

```python
from pysep import Pysep

sep = Pysep(..., order_stations_list_by="distance")
```


-------------------------

##  Pointing PySEP to custom, local databases 

Data are often stored in custom databases that we cannot predict the 
structure of. To point PySEP at your local databases, the simplest method would
be to find a way to read your data and metadata into ObsPy objects, which 
you can then feed into the PySEP machinery. 

```python
from pysep import Pysep
from obspy import read, read_events, read_inventory
st = read()
inv = read_inventory()
event = read_events()[0]
sep = Pysep(st=st, inv=inv, event=event, config_file="config.yaml")
```

