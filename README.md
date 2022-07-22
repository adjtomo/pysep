Python Seismogram Extraction and Processing
===========================================

`PySEP` uses ObsPy tools to request seismic data and metadata, process, 
standardize and format the data, and write out SAC files with the underlying 
motivation of preparing data for moment tensor inversions, although the 
written files can be used for other purposes.

The main processing steps taken by `PySEP` are:

1. Create or gather event metadata (QuakeML) with user-defined event parameters 
2. Gather station metdata (StationXML) in a bounding box surrounding an event, 
3. Curtail station list due to distance or azimuth constraints
4. Gather three-component waveform data for chosen station catalog
5. Quality check waveforms: remove gaps, null out missing components, address 
  clipped amplitudes
6. Preprocess waveforms: detrend, remove response, amplitude scaling, 
  standardizing time series.
7. Rotate streams to desired components: ZNE, RTZ, UVW (triaxial orthogonal)
8. Append SAC headers with event and station metadata, and TauP arrivals
9. Write per-component SAC files, and StationXML, QuakeML and MSEED files
10. Write CAP (Cut-and-Paste) weight files for moment tensor inversions
11. Write config YAML files which can be used to re-run data gathering/processing
12. Plot a record section and source-receiver map

PySEP also has the capacity to:

* Interface with a Lawrence Livermore National Laboratory database of nuclear
  explosion and earthquake waveforms (see note)
* Allow access to embargoed waveform data and PASSCAL HDF5 (PH5) datasets
* Access pre-defined configuration files for data used in previous studies 
* Input custom TauP models for arrival time estimation

### Installation

We recommend installing PySEP into a Conda environment. Dependencies are 
installed via Conda where possible, with Pip used to install PySEP itself. 
This will set two command line tools `pysep` and `recsec`
```bash
$ conda create -n pysep python=3.10
$ conda activate pysep
$ git clone --branch devel https://github.com/uafgeotools/pysep.git
$ cd pysep
$ conda install --file requirements.txt
$ pip install -e .
```

### Running Tests

PySEP comes with unit testing which should be run before and after making any
changes to see if your changes have altered the expected code behavior.
```bash
$ cd tests
$ pytest
```

### Gallery

An example record section produced by the `recsec` tool inside PySEP

![](docs/images/record_section.png)

An example station map generated from collected metadata during a PySEP run

![](docs/images/station_map.png)

--------------------------------------------------------------------------------

## Command Line Usage

Normal users should use PySEP as a command line tool. 

To bring up the command line tool help message:

```bash
$ pysep -h 
```

To list out available pre-defined configuration files

```bash
$ pysep -l  # or pysep --list
-p/--preset -e/--event
-p MTUQ2022_workshop -e 2009-04-07T201255_ANCHORAGE.yaml
-p MTUQ2022_workshop -e 2014-08-25T161903_ICELAND.yaml
-p MTUQ2022_workshop -e 2017-09-03T033001_NORTH_KOREA.yaml
-p MTUQ2022_workshop -e 2020-04-04T015318_SOUTHERN_CALIFORNIA.yaml
```

To run one of the pre-defined configuration files

``` bash
$ pysep -p MTUQ2022_workshop -e 2017-09-03T033001_NORTH_KOREA.yaml 
```

To create a template configuration file that you must fill out with your own
parameters

```bash
$ pysep -W  # or pysep --write
```

To run this newly created configuration file

```bash
$ pysep -c pysep_config.yaml
```

--------------------------------------------------------------------------------
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
$ pysep -c pysep_config.yaml -E event_input.txt
```

--------------------------------------------------------------------------------
### Legacy Filenaming Schema

The new version of PySEP uses a file naming schema that is incompatible with 
previous versions, which may lead to problems in established workflows. 

To honor the legacy naming schema of PySEP, simply use the ``--legacy_naming`` 
argument in the command line. This will change how the event tag is formatted,
how the output directory is structured, and how the output SAC files are named.

```bash
$ pysep -c pysep_config.yaml --legacy_naming
```

--------------------------------------------------------------------------------
### Output Filename Control

The event tag used to name the output directory and written SAC files can be set
manually by the user using the ``--event_tag`` argument. If not given, the tag 
will default to a string consisting of event origin time and Flinn-Engdahl 
region (or just origin time if ``--legacy_naming`` is used). Other output files 
such as the config file and ObsPy objects can be set as in the following: 

```bash
$ pysep -c pysep_config.yaml \
    --event_tag event_abc \
    --config_fid event_abc.yaml \
    --stations_fid event_abc_stations.txt \
    --inv_fid event_abc_inv.xml \
    --event_fid event_abc_event.xml \
    --stream_fid event_abc_st.ms
```

Please note: the output SAC file names are hardcoded and cannot be changed 
by the user. If this is a required feature, please open up a GitHub issue, and 
the developers will address this need.

--------------------------------------------------------------------------------

### Record Section plotter

PySEP also comes with a pretty sophisticated record section tool, which plots
seismic data acquired by PySEP. When you have successfully collected your data,
it will reside in the /SAC folder of the PySEP output directory. 


To see available record section plotting commands

```bash
$ recsec -h  # RECordSECtion
```

To plot the waveform data in a record section with default parameters

```bash
$ recsec --pysep_path ./SAC
```

To plot a record section with a 7km/s move out, high-pass filtered at 1s

```bash
$ recsec --pysep_path ./SAC --move_out 7 --min_period_s 1
```

To sort your record section by azimuth and not distance (default sorting)

```bash
$ recsec --pysep_path ./SAC --sort_by azimuth
```

Have a look at the -h/--help message and the docstring at the top of `recsec.py`
for more options.

#### Plotting SPECFEM synthetics

RecSec can plot SPECFEM-generated synthetic seismograms in ASCII format. 
These can be plotted standalone, or alongside observed seismograms to look at
data-synthetic misfit. 

To access metadata, RecSec requires the CMTSOLUTION and STATIONS file that were 
used by SPECFEM to generate the synthetics. Based on a standard 
SPECFEM3D_Cartesian working directory, plotting synthetics only would have 
the following call structure:

```bash
$ recsec --syn_path OUTPUT_FILES/ --cmtsolution DATA/CMTSOLUTION --stations DATA/STATIONS
```

To compare observed and synthetic data, you would have name the --pysep_path
and tell RecSec to preprocess both data streams identically

```bash
$ recsec --pysep_path ./SAC --syn_path OUTPUT_FILES/ --cmtsolution DATA/CMTSOLUTION --stations DATA/STATIONS --preprocess both
```

Preprocessing flags such as filtering and move out will be applied to both
observed and synthetic data.

--------------------------------------------------------------------------------

## Scripting PySEP

More advanced users can use PySEP as a scripting tool rather than a command 
line tool. 

All of the gathered data/metadata are saved as attributes of the Pysep class 
with typical ObsPy naming schema

```python
>>> from pysep import Pysep
>>> sep = Pysep(config_file='pysep_config.yaml')
>>> sep.run()
>>> print(sep.st[0].stats.sac)
>>> sep.inv
>>> sep.event
```

Although not the preferred method of interacting with PySEP, you can forgo the 
config file and pass parameters directly to the instantiation of the PySEP 
class, making PySEP a bit more flexible.

```python
>>> from pysep import Pysep
>>> sep = Pysep(origin_time="2000-01-01T00:00:00", event_latitude=64.8596,
                event_longitude=-147.8498, event_depth_km=15., ....
                )
```

Check out the Pysep.run() function for other API options for using PySEP.

### Scripting RecSec

The RECord SECtion tool can also be scripted. It simply requires an ObsPy Stream
object as input. Tunable parameters can be fed in as input variables.

```python
>>> from obspy import read
>>> from pysep.recsec import plotw_rs
>>> st = read()
>>> plotw_rs(st=st, sort_by="distance")
```

### Miscellaneous Functionality

#### Append SAC headers to existing Streams

To append SAC headers to your own seismic data, you can directly use the
`PySEP` utility functions `append_sac_headers` and 
`format_sac_header_w_taup_traveltimes`

```python
>>> from pysep.utils.cap_sac import append_sac_headers, format_sac_header_w_taup_traveltimes
>>> from obspy import read, read_events, read_inventory
>>> st = read()
>>> inv = read_inventory()
>>> event = read_events()[0]
>>> st = append_sac_headers(st=st, inv=inv, event=event)
>>> st = format_sac_header_w_taup_traveltimes(st=st, model="ak135")
```

#### Reading in SPECFEM-generated synthetics

PySEP contains a utility function `read_synthetics` to read in 
SPECFEM-generated synthetics with appropriately crafted SAC headers. 
Given a standard SPECFEM3D working directory, reading in SPECFEM synthetics 
might look something like:

```python
>>> from pysep.utils.io import read_synthetics
>>> st = read_synthetics(fid="OUTPUT_FILES/NZ.BFZ.BXE.semd", 
                         cmtsolution="DATA/CMTSOLUTION", 
                         stations="DATA/STATIONS")
>>> print(st)
```

####  Pointing PySEP to custom, local databases 
Data are often stored in custom databases that we cannot predict the 
structure of. To point PySEP at your local databases, the simplest method would
be to find a way to read your data and metadata into ObsPy objects, which 
you can then feed into the PySEP machinery. 

```python
>>> from pysep import Pysep
>>> from obspy import read, read_events, read_inventory
>>> st = read()
>>> inv = read_inventory()
>>> event = read_events()[0]
>>> sep = Pysep(st=st, inv=inv, event=event, config_file="config.yaml")
```



--------------------------------------------------------------------------------

### LLNL Note

`PySEP` interfaces with the databases of:

* W. Walter et al. (2006)
  An assembled western United States dataset for regional seismic analysis
  ISSO 9660 CD, LLNL release UCRL-MI-222502
  
  Download link: https://ds.iris.edu/mda/18-001

