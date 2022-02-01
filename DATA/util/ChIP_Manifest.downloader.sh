#!/usr/bin/env bash
# run in parent directory as util/ChIP_Manifest.downloader.sh

skip_header=1

while IFS=$'\t' read -r -a x; 
do 
    if [ "$skip_header" -eq 1 ]
    then
        skip_header=0
        continue
    fi
        
    if [ -z "$x" ]
    then
        continue
    fi
    url="${x[6]}"
    filename="${x[5]}"
    cmd="wget -O $filename $url"
    echo $cmd
    eval $cmd
done < <(util/exportXLsheet.RScript ChIP_Source_Data_Manifest.xlsx)
