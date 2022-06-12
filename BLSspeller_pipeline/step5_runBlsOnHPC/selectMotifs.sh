#!/bin/bash

# USAGE: selectMotifs.sh <folder that contains "allMotifs.txt" file> <confidence score cutoff>

cd $1
CONF_CUTOFF=$2

module load python/x86_64/3.6.5
find tableOutput/* -type f -name 'part-*' -exec sh -c 'cat $0 >> motifTable.txt' {} \;
OMP_NUM_THREADS=1 python3 /ngsprojects/blsmaize/results/BLSspeller_pipeline/step5_runBlsOnHPC/pyUtils/processBlsOutput.py motifTable.txt . . $CONF_CUTOFF
