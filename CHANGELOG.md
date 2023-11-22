# PySEP Changelog

## Version 0.2.0 

- RecSec plot aesthetics altered to provide more information efficiently  
- Shifted documentation to GitHub Wiki page
- Added moment tensor gathering routines from Pyatoa

## Version 0.3.0 

- Removed llnl_db_client dependency from package and made it optional 
- Re-instated 'mass_downloader' option 
- Added `Declust` class for declustering and weighting
- Completed and re-organized PySEP and RecSec docstrings 
- Bugfix: rotation was not able to be set as null 
- Bugfix: some flag check parameters were not being used during processing
- Removed hard restriction on requiring event depth and magnitude for default event selection
- Add feature 'tmarks' to record section to plot vertical lines at reference times
- Fix ordering of multi-page record sections when using plotw_rs 

## Version 0.3.1 

- Renames package pysep->pysep-adjtomo (only relevant for PyPi, everything in the package remains PySEP)
- Removes llnl_db_client as a dependency of PySEP (PyPi does not allow unpublished dependencies, i.e., via GitHub)
- Adds MANIFEST.in file to retain example config files during source code building and packaging


## Version 0.3.2 

- \#81, \#84
	- Introduce parameters `fill_data_gaps` and `gap_fraction` to address data gaps
	- New hidden parameter `extra_download_pct` handles waveform start and end time 
	- Resampling now occurs before trace start and end time trimming
- Bugfix: bulk query request locations was incorrectly hardocded as a wildcard
- Waveform start and end trim now performed trace by trace with option to fill
  empty boundaries with `fill_data_gaps`
- Introduce new parameter `remove_masked_data` which is a toggle flag to remove
  data which has been merged with no fill value. 

## Version 0.3.3 

- Fixed failing test which was not removing masked data after a trim function
- Improved source receiver map plotting functionality (#91)
- Revamp TauP theoretical arrival time appending to SAC headers (#94)
- RecSec `time_shift_s` now allows shifting by phase arrivals in SAC headers (#94)
- Remove hard no-NoneType restriction on PySEP parameters `water_level` and `output_unit`


## Version 0.4.0 

- Migrates docs from GitHub pages to ReadTheDocs (#102)
- Adds version number to init to allow User to import version
- General Bugfixes (#97, #100, #105)
- **API change** PySEP input Parameters `mindistance` -> `mindistance_km` and 
  `maxdistance` -> `maxdistance_km` (#105)
- File writing occurs throughout PySEP data download and not only at the end,
  and reduced number of default files written. Writes log file (#106)
- Bugfixes not associated with PRs:
  - Fixes multi-page record sections
  - Corrects azimuth definitions for cartesian domains
  - Source-receiver map now handles repeat station names from different networks
- Features not associated with PRs:
  - RecSec kwarg `title` allows overwriting title
  - PySEP warns when config parameters are not used by the program


## Version 0.4.1 
- Adds Tutorial documentation following GEOS626 lab (thanks, Aakash!)
- Adds version release documentation
- Slightly modifies pysep-docs conda environment to accomodate converted nbooks

## Version 0.5.0 
- Improves functions 'read_forcesolution' and 'read_source', which now return
  `obspy.core.event.Event` objects, rather than the makeshift Source objects 
- 'read_forcesolution' can now handle FORCESOLUTION files from both SPECFEM3D
  and SPECFEM3D_GLOBE
- Added function `read_events_plus` that provides additional support to 
  Obspys `read_events` function by allowing for support of FORCESOLUTION and
  SOURCE files from SPECFEM2D/3D/3D_GLOBE
- Remove the `Source` class from `Pysep.utils.io.mt`. This was a remnant of the
  old Pyatoa approach to building an Obspy Event-like object to mimic certain
  behaviors. This has been replaced by read functions which simply return Events
- #117: New quality control function that removes Traces with array length <= 1,
  which would cause preprocessing to fail
- #116: RecSec now logs absmax amplitudes and absmax amplitude ratios IFF both
  `st` and `st_syn` are provided
- #120: Version number is now only sourced from `pyproject.toml`, other 
  locations now reference this file to determine version number
- #124:
  - API Change: RecordSection parameter `cmtsolution` has been **renamed** to 
	`source`. 
  - RecordSection now only expects readable files in  --pysep_path or 
	--syn_path.
	- New `RecordSection.read_data()`function which handles data reading logic 
	and can read both obs data (.SAC from PySEP) and syn data 
	(SPECFEM ASCII files or SAC files)
	- Bugfix: Added an exit catch in RecordSection to stop the workflow if no 
	data is available
- Bugfix: RecSec unable to read different `source` file formats. New parameter
  'srcfmt' allows User to set this manually. If not given, RecSec will attempt
  to guess the file format based on the name of the file.
- Bugfix: `read_specfem3d_cmtsolution_cartesian` was unable to handle 
	Flinn-Engdahl regions that had spaces in them. Also it was unable to 
	handle extra lines in the file. 
- Bugfix #126: Fixed some incorrect parameters in example config files
- PySEP variable namechange: `_legacy_naming` -> 'legacy_naming`
- New RecSec Parameters `wildcard` and `syn_wildcard` to specify how to 
	search for data to plot with RecSec
- Removed unnused parameters 'legacy_naming' and 'log_level' from 
	`Pysep.write_config`

## Version 0.5.1

- Bugfix: RecSec subset streams, which checked that 'st' and 'st_syn' had the same stations, would not run for streams of the same length, leading to edge case where same length streams would plot out of order because they had not been sorted. removed the criteria and now subset streams runs at all times
