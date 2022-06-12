#!/bin/bash

# USAGE: ./processBlsOutput.sh <folder with intermediate files (e.g. named "part-00000", "part-00001", ...)>

cd $1

if [ ! -f ../allMotifs.txt ]; then
    find * -type f -name 'part-*' -exec sh -c 'cat "$0" >> ../allMotifs.txt' {} \;
fi
module load python/x86_64/3.6.5
OMP_NUM_THREADS=1 python3 /ngsprojects/blsmaize/results/BLSspeller_pipeline/step5_runBlsOnHPC/pyUtils/processBlsOutput.py ../allMotifs.txt .. .. None
