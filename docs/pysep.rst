PySEP: data download
====================

The PySEP class is a tool for downloading seismic waveform and metadata 
en-masse. It applies some preprocessing, instrument response removal and
rotation to all data for use in moment tensor and tomographic inversions. 
PySEP can be used via the command line, or through scripting.

.. note::

    See `PySEP class API documentation
    <autoapi/pysep/pysep/index.html#pysep.pysep.Pysep>`__ for a list of
    available input parameters

What it does

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

PySEP via Command Line
----------------------

The most convenient way of using PySEP is as a command line tool. To bring up
the command line tool help message:

.. code:: bash

    pysep -h


To list out available pre-defined configuration files which can be used as an
example of how PySEP operates:

.. code:: bash

    pysep -l  # or pysep --list

    -p/--preset -e/--event
    -p mtuq_workshop_2022 -e 2009-04-07T201255_ANCHORAGE.yaml
    -p mtuq_workshop_2022 -e 2014-08-25T161903_ICELAND.yaml
    -p mtuq_workshop_2022 -e 2017-09-03T033001_NORTH_KOREA.yaml
    -p mtuq_workshop_2022 -e 2020-04-04T015318_SOUTHERN_CALIFORNIA.yaml
    ...

To run one of the pre-defined configuration files to gather and process data:

.. code:: bash

    pysep -p mtuq_workshop_2022 -e 2017-09-03T033001_NORTH_KOREA.yaml


When performing your own data collection, you can create an empty configuration
file that you must fill out with your own parameters:

.. code:: bash

    pysep -W  # or pysep --write


To run this configuration file:

.. code:: bash

    pysep -c pysep_config.yaml


Input Parameters and YAML Parameter File
````````````````````````````````````````

.. note::

    See `PySEP class API documentation
    <autoapi/pysep/pysep/index.html#pysep.pysep.Pysep>`__ for a list of
    available input parameters

When using PySEP from the command line, all input variables are controlled by a 
`YAML <https://yaml.org/>`__ parameter file, which is a text file containing
a set of parameter names and corresponding values, for example:

.. code:: yaml

    # key:       value
    origin_time: '2009-04-07T20:12:55.351000Z'


These parameters control everything from the hypocentral location of 
your earthquake, to the specific waveform data you want to collect, to the
types of preprocessing steps to be applied.

To generate a template parameter file you can run:

.. code:: bash

    pysep -W  # or pysep --write

Scripting PySEP
---------------

.. note::

    See `PySEP class API documentation
    <autoapi/pysep/pysep/index.html#pysep.pysep.Pysep>`__ for a list of
    available input parameters when scripting PySEP.

PySEP can similarly be scripted into other tools. When using a YAML parameter
file, this looks like:

.. code:: python

    from pysep import Pysep
    sep = Pysep(config_file='pysep_config.yaml')
    sep.run()

The results of the `run` function are stored as internal attributes. The most
important attributes are

- `st`: ObsPy Stream object with all waveforms gathered and processed
- `inv`: Station metadata and response information
- `event`: Event object defining the event which the waveforms recorded

.. code:: python

    sep.st
    print(sep.st[0].stats.sac)
    sep.inv
    sep.event


You can also pass parameters directly to the instantiation of the PySEP
class. See the PySEP docstring for input parameter types and definitions.

.. code:: python

    from pysep import Pysep
    sep = Pysep(origin_time="2000-01-01T00:00:00", event_latitude=64.8596,
                event_longitude=-147.8498, event_depth_km=15., ....
                )

-------------------------------------------------------------------------------


PySEP Outputs
--------------

PySEP uses a default directory structure when saving files:

* ``output_dir``: By default, PySEP writes all files to the User-defined
  parameter ``output_dir``, which defaults to the current working directory.
* ``event_tag``: Files are written into a sub-directory defined by the event
  origin time, and Flinn-Engdahl region. For example:
  ``2009-04-07T201255_SOUTHERN_ALASKA``
* ``sac_subdir``: All waveform files are saved in a further sub-directory
  (default: `SAC/`), to avoid cluttering up the output directory. Waveform
  files are tagged by the `event_tag` and trace ID.

Users can use the parameters ``write_files`` and ``plot_files`` to control
exactly what files are produced during a PySEP (see `API documentation
<autoapi/pysep/pysep/index.html#pysep.pysep.Pysep>`__ for details).

By default, PySEP will write SAC files, StationXML, QuakeML and config files,
and create a source-receiver map and record section.

Override Directory Names
````````````````````````

In some cases it may be useful for Users to save files directly to their
working directory, without all the automatically generated sub-directories.

* To ignore the automatically generated event tag, you can set the
  ``overwrite_event_tag`` parameter as an empty string. Via the command line:

  .. code:: bash

      pysep -c pysep_config.yaml --overwrite_event_tag ''

  or via scripting:

  .. code:: python

      sep = Pysep(overwrite_event_tag="")

* To set your own event tag, use a string value

  .. code:: bash

        pysep -c pysep_config.yaml --overwrite_event_tag event_abc

* To ignore the SAC subdirectory and save waveform files directly in the
  output directory, use the ``sac_subdir`` parameter, which should be input in
  your YAML parameter file:

  .. code:: yaml

      sac_subdir: ''

  or via scripting

  .. code:: python

      sep = Pysep(sac_subdir="")  # or use a string value to define your own

* `Example`: if a User only wants to save SAC waveforms for the rotated RTZ
  component within their current working directory, ignoring all automatically
  generated sub directories, all other written files and all plots:

  .. code:: python

      sep = Pysep(overwrite_event_tag="", sac_subdir="", write_files="sac_rtz",
                  plot_files="")


Override Filenames
``````````````````

.. note::

    The output SAC file names are hardcoded as trace IDs (with or without the
    event tag). If additional control over file IDs is a required feature,
    please open up a `GitHub issue <https://github.com/adjtomo/pysep/issues>`__

The event tag used to name the output directory and written SAC files can be set
manually by the user using the ``overwrite_event_tag`` argument.

Other output file names can also be changed from their default values, see the
:meth:`write function <pysep.pysep.Pysep.write>` for write file options and
arguments to use for changing file names.

An example of this via the command line:

.. code:: bash

    pysep -c pysep_config.yaml \
        --overwrite_event_tag event_abc \
        --config_fid event_abc.yaml \
        --stations_fid event_abc_stations.txt \
        --inv_fid event_abc_inv.xml \
        --event_fid event_abc_event.xml \
        --stream_fid event_abc_st.ms

Or with scripting

.. code:: python

    sep = Pysep(overwrite_event_tag="event_abc",
                config_fid="event_abc.yaml", ...)


Legacy Filenaming Schema
````````````````````````

The new version of PySEP uses a file naming schema that is incompatible with
previous versions, which may lead to problems in established workflows.

To honor the legacy naming schema of PySEP, simply use the ``legacy_naming``
parameter. This will change how the event tag is formatted, how the output
directory is structured, and how the output SAC files are named.

.. code:: bash

    pysep -c pysep_config.yaml --legacy_naming

Or with scripting

.. code:: python

    sep = Pysep(legacy_naming=True, ...)


Multiple Event Input
--------------------

To use the same configuration file with multiple events, you can use an event
file passed to PySEP through the command line.

.. note::

    Multiple event input is only available for command line usage of PySEP.
    We suggest using a for loop if you would like to script multiple event
    input using PySEP


When using this option, the event parameters inside the config file will be
ignored, but all the other parameters will be used to gather data and metadata.

Event input files should be text files where each row describes one event with
the following parameters as columns:

.. ORIGIN_TIME LONGITUDE LATITUDE DEPTH[KM] MAGNITUDE

For an example event input file called 'event_input.txt', call structure is:

.. code:: bash

    pysep -c pysep_config.yaml -E event_input.txt


ObsPy Mass Downloader
---------------------

`ObsPy's Mass Download
<https://docs.obspy.org/packages/autogen/obspy.clients.fdsn.mass_downloader.html>`__
feature allows for large data downloads over all available data services. This
may be useful if you don't care where your data comes from and just want to
download all data available.

.. note::

    This ignores the `client` parameter and downloads waveform and station
    metadata for all available data services.

To use the mass download option from the command line, you will first need to
add the following parameter to your YAML config file:

.. code:: yaml

    use_mass_download: true


Then from the command line, you can run things as normal, PySEP will know to
use the mass download option to grab waveform and station metadata. There are
two additional keyword arguments which you can provide  from the command line.

- domain_type (str): Define the search region domain as 
    - `rectangular`: rectangular bounding box defined by `minlatitude`,
       `minlongitude,` `maxlatitude` and `maxlongitude`
    - `circular`: circular bounding circle defined by the `event_latitude`,
      `event_longitude` and min and max radii defined by `mindistance_km` and
      `maxdistance_km`
- delete_tmpdir (bool): Removes the temporary directories that store the MSEED and
  StationXML files which were downloaded by the mass downloader.
  Saves space but also if anything fails prior to saving data,
  the downloaded data will not be saved. Defaults to True.

.. code:: bash

    pysep -c config.yaml --domain_type circular --delete_tmpdir False

Or with scripting

.. code:: python

    sep = Pysep(config_file="config.yaml", use_mass_download=True,
                domain_type="circular", delete_tmpdir=False, ...)



Select Trace ID Input
---------------------

It may be useful to not download data for all possible combinations of
`network`, `station`, `location` and `channel`, but rather to gather data
by trace id only.

PySEP accepts a list of station IDs under the parameter `station_ids`, which
allows users to selectively gather data and metadata for a particular event.

.. note::

    Using the `station_ids` option will **ignore** any values provided to
    `network`, `station`, `location` and `channel`. Trace IDs must be in the
    form: 'NN.SSS.LL.CCC'

.. note::

    Select trace ID input is only available through scripting. If this is
    desirable as a command line input, please open up a GitHub issue.

.. code:: python

    from pysep import Pysep

    station_ids = ["II.DAV.00.LHZ", "IU.XMAS.*.LHZ", "IU.SDV.10.LHZ"]
    sep = Pysep(station_ids=station_ids, ...)
