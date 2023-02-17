# RecSec: record sections

RecordSection (abbreviated RecSec) is a record section tool used to plot 
a number of waveforms on a single figure. Waveforms can be processed, and sorted
by various attributes, such as epicentral distance or azimuth. Data can be 
further modified, for example amplitudes can be scaled by expected geometrical
spreading factor to visualize site amplification effects.

> __Note__: See [RecordSection class API documentation]( 
  https://adjtomo.github.io/pysep/autoapi/pysep/recsec/index.html#pysep.recsec.RecordSection ) for a list of available input parameters

### What it does:
* Plot waveform data and synthetic seismograms (by themselves or together)
* Sort source-receiver pairs by absolute or relative distance or (back)azimuth
* Apply time shifts, zero padding, move out and amplitude scaling
* Scale amplitudes by global maximum, trace maximum or empirical geometric spreading factor
* Preprocess data with filters and detrends prior to plotting 



## Record Section via the Command Line

PySEP also comes with a pretty sophisticated record section tool, which plots
seismic data acquired by PySEP. When you have successfully collected your data,
it will reside in the /SAC folder of the PySEP output directory. 


To see available record section plotting commands

```bash
recsec -h  # RECordSECtion
```

To plot the waveform data in a record section with default parameters

```bash
recsec --pysep_path ./SAC
```

To plot a record section with a 7km/s move out, high-pass filtered at 1s
s
```bash
recsec --pysep_path ./SAC --move_out 7 --min_period_s 1
```

To sort your record section by azimuth and not distance (default sorting)

```bash
recsec --pysep_path ./SAC --sort_by azimuth
```

Have a look at the -h/--help message and the docstring at the top of `recsec.py`
for more options.


> __Note__: See [RecordSection class API documentation]( 
  https://adjtomo.github.io/pysep/autoapi/pysep/recsec/index.html#pysep.recsec.RecordSection ) for a list of available input parameters

### Customizing RecSec figures

Much of the aesthetic look of RecSec figures is hardcoded, however there are 
some keyword arguments that you can provide as flags which may help. The 
following parameters correspond to Matplotlib plot adjustments. 

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

```bash
recsec --pysep_path ./SAC --xtick_minor 100 --xtick_major 500
```


### Plotting SPECFEM synthetics

RecSec can plot SPECFEM-generated synthetic seismograms in ASCII format. Here 
the domain should be defined by geographic coordinates (latitude/longitude). If 
your domain is defined in Cartesian, see below.

Record sections  can be plotted standalone, or alongside observed seismograms 
to look at data-synthetic misfit. 

To access metadata, RecSec requires the CMTSOLUTION and STATIONS file that were 
used by SPECFEM to generate the synthetics. Based on a standard 
SPECFEM3D_Cartesian working directory, plotting synthetics only would have 
the following call structure:

```bash
recsec --syn_path OUTPUT_FILES/ --cmtsolution DATA/CMTSOLUTION --stations DATA/STATIONS
```

To compare observed and synthetic data, you would have name the --pysep_path
and tell RecSec to preprocess both data streams identically

```bash
recsec --pysep_path ./SAC --syn_path OUTPUT_FILES/ --cmtsolution DATA/CMTSOLUTION --stations DATA/STATIONS --preprocess both
```

Preprocessing flags such as filtering and move out will be applied to both
observed and synthetic data.


### Plotting SPECFEM synthetics in Cartesian

Under the hood, RecSec uses some of ObsPy's geodetic
functions to calculate distances and azimuths. Because of this, RecSec will 
fail if coordinates are defined in a Cartesian coordinate system (XYZ), which 
may often be the case when working in SPECFEM2D or in a local domain of 
SPECFEM3D_Cartesian.

To circumvent this, RecSec has a flag `--cartesian`, which will swap out the 
read functions to work with a Cartesian coordinate system. The call is very 
similar to the above:

For SPECFEM3D_Cartesian this would look like

```bash
recsec --syn_path OUTPUT_FILES --cmtsolution DATA/CMTSOLUTION --stations DATA/STATIONS --cartesian
```

For SPECFEM2D, the source file may not be a CMTSOLUTION. Additionally, the 
default seismogram components may not be defined in ZNE

```bash
recsec --syn_path OUTPUT_FILES --cmtsolution DATA/SOURCE --stations DATA/STATIONS --components Y --cartesian
```


### Plotting two sets of synthetics (synsyn)

It is often useful to compare two sets of synthetic seismograms, where one set
represents 'data', while the other represents synthetics. For example, during
a tomographic inversion, a Target model may be used to generate data. 

RecSec can plot two sets of synthetics in a similar vein as plotting 
data and synthetics together. The User only needs to add the `--synsyn` flag 
and provide paths to both `--pysep_path` and `--syn_path`. 

>__NOTE__: RecSec makes the assumption that both sets of synthetics share the 
> same metadata provided in the `--cmtsolution` and `--stations` flags.

Let's say you've stored your 'data' in a directory called 'observed/' and your
synthetics in a directory called 'synthetic/'

```bash
recsec --pysep_path observed/ --syn_path synthetic/ --cmtsolution DATA/CMTSOLUTION --stations DATA/STATIONS --synsyn
```


--------------------------------------------------

## Scripting RecSec

The RECord SECtion tool can also be scripted. It simply requires an ObsPy Stream
object as input. Tunable parameters can be fed in as input variables.

```python
from obspy import read, read_inventory, read_events
from pysep import RecordSection
from pysep.utils.cap_sac import append_sac_headers
st = read()
inv = read_inventory()
event = read_events()[0]
st = append_sac_headers(st, event, inv)  # RecordSection requires SAC headers
rs = RecordSection(st=st, sort_by="distance")
rs.run()
```
