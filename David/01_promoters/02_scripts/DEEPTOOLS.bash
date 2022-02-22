#!/usr/bin/env bash

#                           --samplesLabel 'IDR+,int. enriched' 'IDR+,NOT enriched' 'IDR-,int. enriched' 'IDR-,NOT enriched'\
#                           -S ../01_input/ELT2_LE_combined_subtracted.bw\
#                           --missingDataAsZero 
#                           -bs 100\
#                           -R genes.classA.bed genes.classB.bed genes.classC.bed genes.classD.bed\



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

