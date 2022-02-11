#!/usr/bin/env bash

#                           --samplesLabel 'IDR+,int. enriched' 'IDR+,NOT enriched' 'IDR-,int. enriched' 'IDR-,NOT enriched'\
#                           -S ../01_input/ELT2_LE_combined_subtracted.bw\
#                           --missingDataAsZero 
#                           -bs 100\
#                           -R genes.classA.bed genes.classB.bed genes.classC.bed genes.classD.bed\

if true
then
    computeMatrix scale-regions --regionBodyLength 1200 \
                                --startLabel 'up-1Kb' \
                                --endLabel down+200 \
                                --beforeRegionStartLength 1000\
                                --afterRegionStartLength 500\
                                -R promoters.hilo.classA.bed promoters.hilo.classC.bed promoters.hilo.classB.bed promoters.hilo.classD.bed\
                                -S ELT2_LE_combined_subtracted.interp.bigWig\
                                -p 4 -o promoters.olap100.hilo.mx
fi

plotHeatmap  --matrixFile promoters.olap100.hilo.mx\
             -out promoters.olap100.hilo.pdf\
             --sortRegions no\
             --colorMap RdYlBu_r\
             --startLabel '' --endLabel ''\
             --regionsLabel 'peak+int. enrich.' 'peak+ NOT int. enrich.' 'NO peak + int. enrich.' 'NO peak + NOT int. enrich.'\
             --samplesLabel 'ELT-2 signal (reps. combined subtracted)'
