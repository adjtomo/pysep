#!/bin/sh
#
# combine_cap_weights.sh
#
# usage:
#   sh combine_cap_weights.sh weight_file_aec.dat weight_file_iris.dat
#
# NOTE order of input files!
#
# This script combines CAP weight files from IRIS and AEC (extracted with matlab)
# Main steps
#   1. get stations from IRIS weight file (python tools)
#   2. check if the stations are in the AEC weights (matlab tools)
#   3. If YES then replace station name, use old weights.
#   4. If NOT then use new station name, use new weights.
#
# 20160924 - cralvizuri <cralvizuri@alaska.edu>
#-----------------------------------------------------------

# TEST FILES (make a copy of them!)
#weight_file_aec="/home/vipul/CAP/inv/scak/MOOS/input_files/L1/M111/20090407201255351_weight111.dat"
#weight_file_iris="/home/alvizuri/shared/llnl/waveforms_alaska/20090407201255351/weight.dat"

#weight_file_aec="/home/alvizuri/manuscripts/2015/fmt_uturuncu/uaf_scholarworks/weights_utuhalf_P01_V10_R01_S10/files/20100516063454464_weight_utuhalf_P01_V10_R01_S10.dat"
#weight_file_iris="/home/alvizuri/shared/llnl/waveforms_uturuncu/20100516063454464/weight.dat"

# note order of input files
weight_file_aec=$1
weight_file_iris=$2
out_weight_file=`echo ${weight_file_iris} | sed 's/.dat/_revised.dat/g'`
out_weight_file_save="${out_weight_file}_SAVE"
logfile="${weight_file_iris}_conversion.log"

if [ -e "${out_weight_file}" ] ; then
    echo "*** WARNING. File already exists. Will not overwrite. EXIT ***"
    echo ${out_weight_file}
    exit
fi
printf "" > $out_weight_file

# MAIN
# Check stations in the weight files (python tools, IRIS)
# These are stations generated with the python tools
printf "Stations available from IRIS but not available in AEC weight files\n">\
    ${logfile}
nsta_na_at_aec=""
echo "processing" $weight_file_iris
while read sta_net_iris weights_iris; do
    stname_iris=`echo $sta_net_iris | sed 's/\./ /g' | awk '{print $3}'`
    stnet=`echo $sta_net_iris | sed 's/\./ /g' | awk '{print $2}'`
    stname_pol=`grep ${stname_iris} ${weight_file_aec} | grep $stnet |  awk '{print $1}'`
    if [ -n "$stname_pol" ] ; then
        pol=`echo $stname_pol | sed 's/\// /' | awk '{print $2}'`
        if [ -n "$pol" ]; then
            pol=`printf "/%-3s" ${pol}`
        fi
        #weights_aec=`grep -e \"\<${stname_iris}\>\" ${weight_file_aec} \
            weights_aec=`grep ${stname_iris} ${weight_file_aec} | grep $stnet \
            | awk '{printf "%5s %3s %3s %3s %3s %3s %6s %6s %6s %6s %6s",\
            $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}'`
        printf "%32s ${weights_aec}\n" ${sta_net_iris}${pol} >> ${out_weight_file}
    else
        echo "Warning. IRIS station ${stname_iris} not available in old weight file ${weight_file_aec}" \
            >> ${logfile}
        printf "%32s ${weights_iris}\n" ${sta_net_iris} >> ${out_weight_file}
        nsta_na_at_aec=$((nsta_na_at_aec + 1))
    fi
done < ${weight_file_iris}
cp ${out_weight_file} ${out_weight_file_save}

if [ -n "${nsta_na_at_aec}" ]; then
    echo "Warning. $nsta_na_at_aec IRIS stations not found in the AEC weight files."
    echo "See ${logfile}"
fi

# SECONDARY CHECK
# Check stations in the weight files (matlab tools, AEC databases)
# Check if station is available in the weights from IRIS data.
echo "processing" $weight_file_aec
printf "\nStations available from AEC but not available in IRIS weight files\n"\
    >> ${logfile}
nsta_na_at_iris=""
while read sta_net_aec weights_aec; do
    stname_aec=`echo $sta_net_aec | sed 's/_/ /g' | awk '{print $1}'`
    stname_pol=`echo $sta_net_aec | sed 's/\// /g' | awk '{print $2}'`
    if [ -n "${stname_pol}" ]; then
        stname_pol="/${stname_pol}"
    fi

    # if weights_iris is not empty then match between AEC and IRIS stations
    weights_iris=`grep -e \"\<${stname_aec}\>\" ${weight_file_iris}`
    if [ -n "${weights_iris}" ]; then
        stname_new=`grep -e \"\<${stname_aec}\>\" ${weight_file_iris} | awk '{print $1}'`
        stname_new_pol="${stname_new}${stname_pol}"
    else
        nsta_na_at_iris=$((nsta_na_at_iris + 1))
        echo "Warning.  AEC station ${stname_aec} not available in new weight file ${weight_file_iris}" \
            >> ${logfile}
    fi
done < ${weight_file_aec}

if [ -n "${nsta_na_at_iris}" ]; then
    echo "Warning. $nsta_na_at_iris AEC stations not found in the IRIS weight files."
    echo "See ${logfile}"
fi

echo "Done. New weight file: ${out_weight_file}"
echo
