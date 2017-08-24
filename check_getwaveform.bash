#!/bin/bash
#
# Run this script to compare the number of extracted files with those in the saved copy 
#
# Warning: If the example crashed before it writes out a file list,
#          then no comparison of files can be made.
#

# all events
events=("20090407201255351" "20160124103037400" "20140831030657111" "20151106012012712" "20160729211826000" "20151106012012712" "20080418093659110" "20100516063454464" "20160518032548320" "HOYA")
# indices in event_input.py
iex=("0" "1" "2" "3" "4" "5" "6" "7" "8" "11")

nevent=${#events[@]}
imax=$(($nevent - 1))
echo "$nevent events in check_getwaveform.bash (index 0 to $imax)"

# all events except HOYA
imin=0
imax=8
# HOYA (requires LLNL client)
#imin=9
#imax=9
# all events including HOYA
#imin=0
#imax=9

ncheck=$(($imax - $imin + 1))
echo "checking $ncheck events in check_getwaveform.bash (index $imin to $imax)"

save_dir=$GEOTOOLS'/python_util/util_data_syn/getwaveform_saved/check_filenames/'

# Run example
for ii in `seq $imin $imax`
do
    python run_getwaveform.py event_input ${iex[$ii]}
done

# Check number of files generated
for ii in `seq $imin $imax`
do
    diff $save_dir/${events[$ii]}_all_filenames ${events[$ii]}/${events[$ii]}_all_filenames
    error=$?
    if [ $error -eq 0 ]
    then
	echo "$ii SUCCESS ${events[$ii]}: Extracted number of files match with the saved version"
    else
	echo "$ii WARNING: file1 and file2 differ"
    fi
done
