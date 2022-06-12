#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_vmem=5G
module load nextflow
nextflow -DXmx=5G -C blsInputPreparation_Pac.config run blsInputPreparation.nf
