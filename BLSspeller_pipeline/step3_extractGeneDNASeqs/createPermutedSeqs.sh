#!/bin/bash

# USAGE: ./createPermutedSeqs.sh <output folder with "masked and "unmasked" dir>
# -> in the output folder, an extra folder "permutedSeqs" will be created

cd $1
mkdir permutedSeqs
cd -
module load python/x86_64/3.6.5
python3 PyUtils/createPermutedSeqs.py $1
