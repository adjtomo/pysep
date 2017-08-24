#!/bin/bash
#
# Run this script to compare the number of extracted files with those in the saved copy 
#
# Warning: If the example crashed before it writes out a file list,
#          then no comparison of files can be made.
#

# all events
events=("20080418093659110" "20100516063454464" "20150928115012000" "20160124103037400" "20090407201255351" "20140831030657111" "20151106012012712" "20160729211826000" "HOYA")
# indices in event_input.py
iex="7 8 6 1 0 2 3 5 4"
nevent=${#events[@]}
nmax=$(($nevent - 1))
echo "$nevent events in check_getwaveform.bash (index 0 to $nmax)"

# all events except HOYA
nmin=0
nmax=7
# HOYA (requires LLNL client)
#nmin=8
#nmax=8
# all events including HOYA
#nmin=0
#nmax=8

ncheck=$(($nmax - $nmin + 1))
echo "checking $ncheck events in check_getwaveform.bash (index $nmin to $nmax)"

save_dir='~/REPOSITORIES/GEOTOOLS/python_util/util_data_syn/getwaveform_saved/check_filenames/'

# Run example
for ii in `seq $nmin $nmax`
do
    python run_getwaveform.py event_input $iex[$ii]
done

# Check number of files generated
for ii in `seq $nmax $nmax`
do
    diff $save_dir/${events[$ii]}_all_filenames ${events[$ii]}/${events[$ii]}_all_filenames
    error=$?
    if [ $error -eq 0 ]
    then
	echo "SUCCESS ${events[$ii]}: Extracted number of files match with the saved version"
    else
	echo "WARNING: file1 and file2 differ"
    fi
done
