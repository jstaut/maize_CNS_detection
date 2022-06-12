#!/bin/bash

module load python/x86_64/3.6.5
OMP_NUM_THREADS=1 python3 clusterSelection.py /ngsprojects/blsmaize/results/BLSspeller_pipeline/step2_selectClusters/input_Pac/Orthogroups_cleaned.tsv /ngsprojects/blsmaize/results/BLSspeller_pipeline/step2_selectClusters/results_Pac True 4
