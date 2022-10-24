#!/usr/bin/env bash
# run in parent directory as util/ChIP_Manifest.downloader.sh
# Run an Rscript to parse the manifest Excel sheet,
# download the parsed data using wget

skip_header=1

function debug {
    continue
    echo -e  $@ >&2
}

linenum=0
while IFS=$'\t' read -r -a x; 
do 
    linenum=$((linenum+1))

    if [ "$skip_header" -eq 1 ]
    then
        skip_header=0
        continue
    fi
        
    if [ -z "$x" ]
    then
        debug "Skipping blank line $linenum"
        continue
    fi

    if [ ${#x[@]} -lt 6 ]
    then
        debug "Error for line $linenum: not enough fields"
        continue
    fi

    url="${x[6]}"

    if [ -z "$url" ]
    then
        debug "Error for line $linenum: URL missing for '$filename'"
        debug "0: ${x[0]}"
        debug "1: ${x[1]}"
        debug "2: ${x[2]}"
        debug "3: ${x[3]}"
        debug "4: ${x[4]}"
        debug "5: ${x[5]}"
        debug "6: ${x[6]}"
    fi

    echo -e "URL:\t$url" 

    filename="${x[5]}"
    echo -e "FNAME:\t$filename" 

    if [ -z "$filename" ]
    then
        debug "Error for line $linenum: Missing filename $filename"
        continue
    fi
    if ! [ -e "$filename" ]
    then
        cmd="wget -O $filename $url"
        echo "Running: $cmd"
        eval $cmd
    fi
    echo
done < <(util/exportXLsheet.RScript ChIP_Source_Data_Manifest.xlsx)
