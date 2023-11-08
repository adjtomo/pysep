RecSec: record sections
=======================

RecordSection (abbreviated RecSec) is a record section tool used to plot 
a number of waveforms on a single figure. Waveforms can be processed, and sorted
by various attributes, such as epicentral distance or azimuth. Data can be 
further modified, for example amplitudes can be scaled by expected geometrical
spreading factor to visualize site amplification effects.

.. note::

    See `RecordSection class API documentation
    <autoapi/pysep/recsec/index.html#pysep.recsec.RecordSection>`__ for a list
    of available input parameters

What it does:

* Plot waveform data and synthetic seismograms (by themselves or together)
* Sort source-receiver pairs by absolute or relative distance or (back)azimuth
* Apply time shifts, zero padding, move out and amplitude scaling
* Scale amplitudes by global maximum, trace maximum or empirical geometric spreading factor
* Preprocess data with filters and detrends prior to plotting 


Record Section via the Command Line
-----------------------------------

To see available record section plotting commands

.. code:: bash

    recsec -h  # RECordSECtion


To plot the waveform data in a record section with default parameters

.. code:: bash

    recsec --pysep_path ./SAC

To plot a record section with a 7km/s move out, high-pass filtered at 1s

.. code:: bash

    recsec --pysep_path ./SAC --move_out 7 --min_period_s 1


To sort your record section by azimuth and not distance (default sorting)

.. code:: bash

    recsec --pysep_path ./SAC --sort_by azimuth

Have a look at the -h/--help message and the `RecordSection class API
documentation <autoapi/pysep/recsec/index.html#pysep.recsec.RecordSection>`__
for more options.


Scripting RecSec
----------------
The RECord SECtion tool can also be scripted. You can feed it the same
`pysep_path` as in the command line usage, or directly inject ObsPy Streams.

The `plotw_rs` function can be used to call RecSec with some logic that allows
for multi-page record sections.

.. note::

    If using ObsPy Streams, all traces must have SAC headers with appropriate
    event and station locations. It is recommended that data fed into RecSec
    is gathered directly via PySEP to ensure appropriate SAC header values.

.. code:: python

    from obspy import read, read_inventory, read_events
    from pysep import plotw_rs
    from pysep.utils.cap_sac import append_sac_headers

    st = read()
    inv = read_inventory()
    event = read_events()[0]
    st = append_sac_headers(st, event, inv)  # RecordSection requires SAC header

    plotw_rs(st=st, sort_by="distance", scale_by="normalize")

Have a look at the `RecordSection class API
documentation <autoapi/pysep/recsec/index.html#pysep.recsec.RecordSection>`__
for available options.

Subsetting Stations
```````````````````

When plotting Streams, it may be useful to subset data, e.g., by network.
This can be achieved directly using ObsPy's Stream object

.. code:: python

    st_new = st.select(network="II")
    plotw_rs(st=st)

.. note::

    This can be combined with source receiver map plotting to get a station
    map of the subsetted traces in the record section. See the
    `map plotting documentation
    <cookbook.html#generating-source-receiver-maps>`__ for an example of how to
    do this.


Customizing RecSec figures
--------------------------

Much of the aesthetic look of RecSec figures is hardcoded, however there are 
some keyword arguments that you can provide as flags which may help to achieve
publication-ready figures. Some of these parameters include:

- ytick_fontsize: Fontsize of the Y-axis tick labels
- xtick_fontsize: Fontsize of the X-axis tick labels
- tick_linewidth: Linewidth of axes tick marks
- tick_length: Length of axes tick marks
- label_fontsize: Fontsize of X and Y axis labels
- axis_linewidth: Linewidth of border around figure
- title_fontsize: Fontsize of the title
- xtick_minor: Frequency of minor ticks on the X axis
- xtick_major: Frequency of major ticks on the X axis

To set one of these parameters, just set as a flag, e.g.,

.. code:: bash

    recsec --pysep_path ./SAC --xtick_minor 100 --xtick_major 500

Or when scripting RecSec

.. code:: python

    plotw_rs(pysep_path="./SAC", xtick_minor=100, xtick_major=500)

See :meth:`set_plot_aesthetic <pysep.utils.plot.set_plot_aesthetic>` function
for the entire list of available tuning options for Record Sections.

Plotting SPECFEM synthetics
---------------------------

RecSec can plot SPECFEM-generated synthetic seismograms in ASCII format. Here 
the domain should be defined by geographic coordinates (latitude/longitude). If 
your domain is defined in Cartesian, see the next section.

.. note::

    Record sections  can be plotted standalone, or alongside observed seismograms
    to look at data-synthetic misfit.

To access metadata, RecSec requires the `CMTSOLUTION` and `STATIONS` file that
are used by SPECFEM to generate the synthetics. Based on a standard SPECFEM
working directory, plotting synthetics only would have the following call
structure:

.. code:: bash

    recsec --syn_path OUTPUT_FILES/ --source DATA/CMTSOLUTION --stations DATA/STATIONS

Or when scripting,

.. code:: python

    plotw_rs(syn_path="OUTPUT_FILES", source="DATA/CMTSOLUTION",
             stations="DATA/STATIONS")

You can also directly feed in an ObsPy stream containing your synthetic data
with appropriate SAC headers

.. code:: python

    from glob import glob
    from obspy import Stream
    from pysep.utils.io import read_sem

    st_syn = Stream()
    for fid in glob("OUTPUT_FILES/*.semd"):
        st_syn += read_sem(fid=fid, source="DATA/CMTSOLUTION",
                           stations="DATA/STATIONS")

    plotw_rs(st_syn=st_syn)

To compare observed and synthetic data, you would have name the --pysep_path
and tell RecSec to preprocess both data streams identically

.. code:: bash

    recsec --pysep_path ./SAC --syn_path OUTPUT_FILES/ --source DATA/CMTSOLUTION --stations DATA/STATIONS --preprocess both

Preprocessing flags can be applied to the observed data only (`st`), synthetic
data only (`st_syn`) or both (`both`). See the `preprocess` parameter in the
`RecordSection class API
documentation <autoapi/pysep/recsec/index.html#pysep.recsec.RecordSection>`__

While scripting, Streams for both observed and synthetic data can be injected
together:

.. code:: python

    plotw_rs(st=st, st_syn=st_syn, preprocess="both")


Plotting SPECFEM synthetics in Cartesian
````````````````````````````````````````

Under the hood, RecSec uses some of ObsPy's geodetic
functions to calculate distances and azimuths. Because of this, RecSec will 
fail if coordinates are defined in a Cartesian coordinate system (XYZ), which 
may often be the case when working in SPECFEM2D or in a local domain of 
SPECFEM3D_Cartesian.

To circumvent this, RecSec has a flag `--cartesian`, which will swap out the 
read functions to work with a Cartesian coordinate system. The call is very 
similar to the above:

For SPECFEM3D_Cartesian this would look like

.. code:: bash

    recsec --syn_path OUTPUT_FILES --source DATA/CMTSOLUTION --stations DATA/STATIONS --cartesian


For SPECFEM2D, the source file may not be a CMTSOLUTION. Additionally, the 
default seismogram components may not be defined in ZNE

.. code:: bash

    recsec --syn_path OUTPUT_FILES --source DATA/SOURCE --stations DATA/STATIONS --components Y --cartesian


While scripting, the input parameter `cartesian` can be used:

.. code:: python

    plotw_rs(..., cartesian=True)


Plotting two sets of synthetics (synsyn)
````````````````````````````````````````

It is often useful to compare two sets of synthetic seismograms, where one set
represents 'data', while the other represents synthetics. For example, during
a tomographic inversion, a Target model may be used to generate data. 

RecSec can plot two sets of synthetics in a similar vein as plotting 
data and synthetics together. The User only needs to add the `--synsyn` flag 
and provide paths to both `--pysep_path` and `--syn_path`. 

.. note::

    RecSec makes the assumption that both sets of synthetics share the
    same metadata provided in the `--source` and `--stations` flags.

Let's say you've stored your 'data' in a directory called 'observed/' and your
synthetics in a directory called 'synthetic/'

.. code:: bash

    recsec --pysep_path observed/ --syn_path synthetic/ --source DATA/CMTSOLUTION --stations DATA/STATIONS --synsyn




