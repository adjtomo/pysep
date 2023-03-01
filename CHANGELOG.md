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

## Version 0.3.3 (Devel)

- Source receiver map plot now plots event first so that figures are centered
  on the epicenter and not on the middle of the stations



