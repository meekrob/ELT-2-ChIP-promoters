#!/usr/bin/env bash


orig_bw=$1
base=${orig_bw%.*}
bedGraph=${base}.bedGraph
if ! [ -e $bedGraph ]
then
    cmd="bigWigToBedGraph $orig_bw $bedGraph"
    echo "Creating bedGraph file: $cmd"
    eval $cmd
fi 

interp=${base}.interp.bedGraph

scriptdir=$(dirname ${BASH_SOURCE[0]})
cmd="python $scriptdir/remove_duplicated.py <(python $scriptdir/connect_encode_dots.py $bedGraph) > $interp"
echo $cmd
time eval $cmd

# convert back to bigwig to save space (delete bg file on success)
interp_bw=${interp/.bedGraph/.bw}
cmd="bedGraphToBigWig $interp <(fetchChromSizes ce11) $interp_bw && rm $interp"
echo $cmd
time eval $cmd

