#!/usr/bin/env bash

#                           --samplesLabel 'IDR+,int. enriched' 'IDR+,NOT enriched' 'IDR-,int. enriched' 'IDR-,NOT enriched'\
#                           -S ../01_input/ELT2_LE_combined_subtracted.bw\
#                           --missingDataAsZero 
#                           -bs 100\
#                           -R genes.classA.bed genes.classB.bed genes.classC.bed genes.classD.bed\

# For the class A,B,C,D breakdown

UPSTREAM=$1     # i.e. 1000
DOWNSTREAM=$2   # i.e. 200
SIGNAL=$3       # i.e. ELT2_LE_combined_subtracted.interp.bigWig


if false
then
    computeMatrix scale-regions --regionBodyLength 1200 \
                                --startLabel 'up-1Kb' \
                                --endLabel down+200 \
                                --beforeRegionStartLength 1000\
                                --afterRegionStartLength 500\
                                -R promoters.hilo.classA.bed promoters.hilo.classC.bed promoters.hilo.classB.bed promoters.hilo.classD.bed\
                                -S ELT2_LE_combined_subtracted.interp.bigWig\
                                -p 4 -o promoters.olap100.hilo.mx

    plotHeatmap  --matrixFile promoters.olap100.hilo.mx\
                 -out promoters.olap100.hilo.pdf\
                 --sortRegions no\
                 --colorMap RdYlBu_r\
                 --startLabel '' --endLabel ''\
                 --regionsLabel 'peak+int. enrich.' 'peak+ NOT int. enrich.' 'NO peak + int. enrich.' 'NO peak + NOT int. enrich.'\
                 --samplesLabel 'ELT-2 signal (reps. combined subtracted)'
fi

# For hi versus lo
set -x
if true
then
    computeMatrix scale-regions --regionBodyLength 1200 \
                                --startLabel 'up-1Kb' \
                                --endLabel down+200 \
                                --beforeRegionStartLength 1000\
                                --afterRegionStartLength 500\
                                -R promoters.hilo.up.bed promoters.hilo.down.bed\
                                -S ELT2_LE_combined_subtracted.interp.bigWig\
                                -p 4 -o promoters.hilo.updown.mx

    plotHeatmap  --matrixFile promoters.hilo.updown.mx\
                 -out promoters.updown.pdf\
                 --sortRegions no\
                 --colorMap RdYlBu_r\
                 --startLabel '' --endLabel ''\
                 --regionsLabel 'log2FC > 0' 'log2FC < 0'\
                 --samplesLabel 'ELT-2 signal (reps. combined subtracted)'
fi

