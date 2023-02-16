# Declust: catalog declustering

>__WARNING__: This docs page is still under construction

`Declust` is a core module which can decluster event catalogs and apply
geographical weighting to events and stations. These tools are mainly used for
curating input data for tomographic inversions.

## Motivation (Seismic Tomography)

Declustering: Say you gather a general catalog of events and a list of stations for a given region using ObsPy. 
If you directly used this catalog without any many of the simulations run have the possibility of being redundant, 
as clusters of seismicity (e.g., foreshock aftershock sequences) tend to sample the same spatial area, providing a 
limited amount of new knowledge at the price of an entire waveform simulation (expensive!).


Weighting: For this same dataset, let's say 90% of your stations are focused in one part of your domain, but you 
have a select few (10%) that are geographically isolated. If all equally weighted, then the 90% of stations will
tend to overshadow the 10%. Source receiver weighting helps alleviate this by penalizing stations or events that
have many nearby counterparts, and promoting more geographically isolated sources.

## Scripting Declust

`Declust` must be scripted or used in a Python interactive environment, it currently does not have a command line interface. 

### Initializing

Declust requires two inputs, an Obspy `Catalog` defining the events you would like to include, and an ObsPy `Inventory` defining 
stations.

```python
from pysep import Declust
from obspy import read_events, read_inventory

# Example data
cat = read_events()
inv = read_inventory()

declust = Declust(cat=cat, inv=inv)
```

Internally, `Declust` will comb through the Catalog and Inventory for relevant metadata (locations, IDs, data availability).  
Users may provide their own Catalog and Inventory gathered via ObsPy or other methods. 

### Data Availability

`Declust` calculates data availability based on station uptime. That is, for each event and its corresponding origin time, 
`Declust` will check if each station in the Inventory was "on" during the event origin time. 

```python
declust.data_avail
```
```bash
{'quakeml:eu.emsc/event/20120404_0000041': ['GR.FUR', 'GR.WET', 'BW.RJOB'], 'quakeml:eu.emsc/event/20120404_0000038': ['GR.FUR', 'GR.WET', 'BW.RJOB'], 'quakeml:eu.emsc/event/20120404_0000039': ['GR.FUR', 'GR.WET', 'BW.RJOB']}
```

Users may provide their own data availability dictionary to `Declust` that can be based on more sophisticated measurements,
such as the presence of a P-wave arrival, or a measurement window compared with synthetics. The keys of such a dictionary
must match the resource ID of each event in the catalog (i.e., `event.resource_id.id`) and the values are lists of 
stations which are 'available', and in the format NN.SSS (N=network, S=station).

### Threshold Catalog

`Declust` has a function `threshold_catalog` which removes events in the Catalog based on depth values. This is very similar 
to the `Catalog.filter()` function provided natively in ObsPy. Users input a number of depth boundaries (units: km) and the
minimum magnitude and minimum data availability allowed for all Events within those boundaries. 

This function updates the internal Catalog attribute of the `Declust` class

```python
declust.threshold_catalog(zedges=[0, 10, 30, 100], min_mags=[4, 5, 6], min_data=10)
```

The above example partitions all events in the Catalog between depths of 0-10km, 10-30km and 30-100km. Within each depth range, 
it defines the minimum magnitude allowed (e.g., M4 for 0-10km), and the minimum data availability (i.e., an event must have 
atleast 10 stations 'available' to be considered).

### Declustering

Declustering breaks up the domain into a number of cells and groups events based on which cell they are in. Within each cell 
a group of events are sorted by a specific criteria (e.g., magnitude, depth), and only a chosen number of events are retained.





