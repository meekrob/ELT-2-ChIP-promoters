#!/usr/bin/env bash
set -e # die on error

#                           --samplesLabel 'IDR+,int. enriched' 'IDR+,NOT enriched' 'IDR-,int. enriched' 'IDR-,NOT enriched'\
#                           -S ../01_input/ELT2_LE_combined_subtracted.bw\
#                           --missingDataAsZero 
#                           -bs 100\
#                           -R genes.classA.bed genes.classB.bed genes.classC.bed genes.classD.bed\

# For the class A,B,C,D breakdown

UPSTREAM=$1     # i.e. 1000
DOWNSTREAM=$2   # i.e. 200
SIGNAL=$3       # i.e. ELT2_LE_combined_subtracted.interp.bigWig

STAGE=$1

set -x
if true
then
    computeMatrix scale-regions --regionBodyLength 1200 \
                                --startLabel 'up-1Kb' \
                                --endLabel down+200 \
                                --beforeRegionStartLength 1000\
                                --afterRegionStartLength 200\
                                -S ../01_input/ELT2_${STAGE}_combined_subtracted.interp.bw\
                                -R $STAGE.class{A,B,C,D}.v2.bed\
                                -p 4 -o $STAGE.classes.v2.mx
                                #-R $STAGE.classA.v2.bed\
                                #-p 4 -o $STAGE.classA.v2.mx

                 #--binSize 1\
                 #$STAGE.classA.v2.mx\
                 #-out $STAGE.classA.v2.pdf\
                 #--regionsLabel 'peak+int. enrich.' \
                 #--colorMap RdYlBu_r\
fi
#cbind.classes.mx\
if true
then
    plotHeatmap  --matrixFile \
                 $STAGE.classes.v2.mx\
                 --samplesLabel 'ELT-2 signal (reps. combined subtracted)'\
                 -out $STAGE.classes.v2.pdf\
                 --sortRegions no\
                 --colorMap Greys\
                 --startLabel '' --endLabel ''\
                 --regionsLabel 'peak+int. enrich.' 'peak+ NOT int. enrich.' 'NO peak + int. enrich.' 'NO peak + NOT int. enrich.'\
                 --yMin 0\
                 --yMax 450\
                 --zMin 0\
                 --zMax 450
fi
