# PySEP Changelog

## Version 0.2.0 (Master)

- RecSec plot aesthetics altered to provide more information efficiently  
- Shifted documentation to GitHub Wiki page
- Added moment tensor gathering routines from Pyatoa

## Version 0.3.0 (Devel)

- Removed llnl_db_client dependency from package and made it optional 
- Re-instated 'mass_downloader' option 
- Added `Declust` class for declustering and weighting
- Completed and re-organized PySEP and RecSec docstrings 
- Bugfix: rotation was not able to be set as null 
- Bugfix: some flag check parameters were not being used during processing
- Removed hard restriction on requiring event depth and magnitude for default event selection

