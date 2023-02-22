# PySEP: data download

The PySEP class is a tool for downloading seismic waveform and metadata 
en-masse. It applies some preprocessing, instrument response removal and
rotation to all data for use in moment tensor and tomographic inversions. 
PySEP can be used via the command line, or through scripting.

> __Note__: See [PySEP class API documentation](
  https://adjtomo.github.io/pysep/autoapi/pysep/pysep/index.html#pysep.pysep.Pysep)
  for a list of available input parameters

### What it does
1. Create or gather event metadata (QuakeML) with user-defined event parameters 
2. Gather station metdata (StationXML) in a bounding box surrounding an event, 
3. Curtail station list due to distance or azimuth constraints
4. Gather three-component waveform data for chosen station catalog
5. Quality check waveforms: remove gaps, null out missing components, address 
  clipped amplitudes
6. Preprocess waveforms: detrend, remove response, amplitude scaling, 
  standardizing time series
7. Rotate streams to desired components: ZNE, RTZ, UVW (triaxial orthogonal)
8. Append SAC headers with event and station metadata, and TauP arrivals
9. Write per-component SAC files, plus StationXML, QuakeML and MSEED files
10. Write CAP (Cut-and-Paste) weight files for moment tensor inversions
11. Write config YAML files which can be used to reproduce data gathering/processing


--------------------------------------------------------------------------------

## PySEP via Command Line

Normal users should use PySEP as a command line tool. 

To bring up the command line tool help message:

```bash
pysep -h 
```

To list out available pre-defined configuration files

```bash
pysep -l  # or pysep --list

-p/--preset -e/--event
-p MTUQ2022_workshop -e 2009-04-07T201255_ANCHORAGE.yaml
-p MTUQ2022_workshop -e 2014-08-25T161903_ICELAND.yaml
-p MTUQ2022_workshop -e 2017-09-03T033001_NORTH_KOREA.yaml
-p MTUQ2022_workshop -e 2020-04-04T015318_SOUTHERN_CALIFORNIA.yaml
```

To run one of the pre-defined configuration files

``` bash
pysep -p MTUQ2022_workshop -e 2017-09-03T033001_NORTH_KOREA.yaml 
```

To create a template configuration file that you must fill out with your own
parameters

```bash
pysep -W  # or pysep --write
```

To run this newly created configuration file

```bash
pysep -c pysep_config.yaml
```


### Input Parameters and YAML Parameter File

> __Note__: See [PySEP class API documentation](
  https://adjtomo.github.io/pysep/autoapi/pysep/pysep/index.html#pysep.pysep.Pysep)
  for a list of available input parameters

When using PySEP from the command line, all input variables are controlled by a 
[YAML](https://yaml.org/) parameter file, which is a text file containing 
a set of parameter names and corresponding values, for example:

```yaml
# key: value
origin_time: '2009-04-07T20:12:55.351000Z'
```

These parameters control everything from the hypocentral location of 
your earthquake, to the specific waveform data you want to collect, to the
types of preprocessing steps to be applied.

To generate a template parameter file you can run:

```bash
pysep -W  # or pysep --write
```


### Multiple Event Input

To use the same configuration file with multiple events, you can use an event 
file passed to PySEP through the command line.

When using this option, the event parameters inside the config file will be
ignored, but all the other parameters will be used to gather data and metadata.

Event input files should be text files where each row describes one event with 
the following parameters as columns:

> ORIGIN_TIME LONGITUDE LATITUDE DEPTH[KM] MAGNITUDE

For an example event input file called 'event_input.txt', call structure is:

```bash
pysep -c pysep_config.yaml -E event_input.txt
```

--------------------------------------------------------------------------------
### Legacy Filenaming Schema

The new version of PySEP uses a file naming schema that is incompatible with 
previous versions, which may lead to problems in established workflows. 

To honor the legacy naming schema of PySEP, simply use the ``--legacy_naming`` 
argument in the command line. This will change how the event tag is formatted,
how the output directory is structured, and how the output SAC files are named.

```bash
pysep -c pysep_config.yaml --legacy_naming
```


### Output Filename Control

The event tag used to name the output directory and written SAC files can be set
manually by the user using the ``--event_tag`` argument. If not given, the tag 
will default to a string consisting of event origin time and Flinn-Engdahl 
region (or just origin time if ``--legacy_naming`` is used). Other output files 
such as the config file and ObsPy objects can be set as in the following: 

```bash
pysep -c pysep_config.yaml \
    --overwrite_event_tag event_abc \
    --config_fid event_abc.yaml \
    --stations_fid event_abc_stations.txt \
    --inv_fid event_abc_inv.xml \
    --event_fid event_abc_event.xml \
    --stream_fid event_abc_st.ms
```

Please note: the output SAC file names are hardcoded and cannot be changed 
by the user. If this is a required feature, please open up a GitHub issue, and 
the developers will address this need.


### ObsPy Mass Downloader

[ObsPy's Mass Download](https://docs.obspy.org/packages/autogen/obspy.clients.fdsn.mass_downloader.html)
feature allows for large data downloads over all available data services. This may be useful if 
you don't care where your data comes from and just want to download all data available. This 
ignores the ``client`` option and downloads waveform and station metadata for all available 
data services.

To use the mass download option from the command line, you will first need to add the following
parameter to your YAML config file:

```yaml
use_mass_download: true     
```

Then from the command line, you can run things as normal, PySEP will know to use the mass download option 
to grab waveform and station metadata. There are two additional keyword arguments which you can provide 
from the command line. 

- domain_type (str): Define the search region domain as 
    - ``rectangular``: rectangular bounding box defined by `minlatitude`,
       `minlongitude,` `maxlatitude` and `maxlongitude`
    - ``circular``: circular bounding circle defined by the `event_latitude`,
      `event_longitude` and min and max radii defined by `mindistance` and 
      `maxdistance`
- delete_tmpdir (bool): Removes the temporary directories that store the MSEED and
  StationXML files which were downloaded by the mass downloader.
  Saves space but also if anything fails prior to saving data,
  the downloaded data will not be saved. Defaults to True.

```bash
pysep -c config.yaml --domain_type circular --delete_tmpdir False
```
-------------------------------------------------------------------------------

## Scripting PySEP

Check out the Cookbook section for additional useful functions and routines
that can be implemented.

All of the gathered data/metadata are saved as attributes of the Pysep class 
with typical ObsPy naming schema

```python
from pysep import Pysep
sep = Pysep(config_file='pysep_config.yaml')
sep.run()
print(sep.st[0].stats.sac)
sep.inv
sep.event
```

You can also pass parameters directly to the instantiation of the PySEP 
class. See the PySEP docstring for input parameter types and definitions.

```python
from pysep import Pysep
sep = Pysep(origin_time="2000-01-01T00:00:00", event_latitude=64.8596,
            event_longitude=-147.8498, event_depth_km=15., ....
            )
```

> __Note__: See [PySEP class API documentation](
  https://adjtomo.github.io/pysep/autoapi/pysep/pysep/index.html#pysep.pysep.Pysep)
  for a list of available input parameters
