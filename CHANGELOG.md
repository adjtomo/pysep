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

- #81, #84
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

## Version 0.4.0 (Master/Devel)
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

