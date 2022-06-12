#!/bin/bash

# USAGE: selectBEDs.sh <folder with BED files> <selection file e.g. conf95_filtered_selectedMotifs.txt>

mkdir ${1}_filtered

module load python/x86_64/3.6.5
OMP_NUM_THREADS=1 python3 /ngsprojects/blsmaize/results/BLSspeller_pipeline/step5_runBlsOnHPC/pyUtils/filterBEDFiles.py $1 $2 ${1}_filtered
