#!/bin/bash

# USAGE: getUniqueTrueAndShuffledMotifs.sh <path to "confXX_selectedMotifs.txt"> <folder with shuffled motifs: confXX_selectedMotifs.txt>

module load python/x86_64/3.6.5
python3 /ngsprojects/blsmaize/results/BLSspeller_pipeline/step5_runBlsOnHPC/pyUtils/substractShuffledOverlap.py $1 $2
