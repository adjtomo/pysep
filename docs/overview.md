Python Seismogram Extraction and Processing
===========================================

## What is PySEP?

`PySEP` provides a wrapper on common ObsPy data gathering routines to request seismic data and metadata for a given seismic event or events.   It also features routines to process, standardize and format waveform data, ultimately writing out SAC files ready for moment tensor or waveform inversions. 

## What does PySEP do?

A standard PySEP data gathering workflow proceeds as follows:

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

PySEP can also

* Interface with a Lawrence Livermore National Laboratory database of nuclear
  explosion and earthquake waveforms (see note)
* Allow access to embargoed waveform data and PASSCAL HDF5 (PH5) datasets
* Access pre-defined configuration files for data used in previous studies 
* Input custom TauP models for arrival time estimation

## What is RecSec?

RecSec is a record section plotter bundled with the PySEP package. It features a suite of options that allow the 
User to control the details of the record section.

* Plot waveform data and synthetic seismograms (by themselves or together)
* Sort source-receiver pairs by absolute or relative distance or (back)azimuth
* Apply time shifts, zero padding, move out and amplitude scaling
* Scale amplitudes by global maximum, trace maximum or empirical geometric spreading factor
* Preprocess data with filters and detrends prior to plotting 

