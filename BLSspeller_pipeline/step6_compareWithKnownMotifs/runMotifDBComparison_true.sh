#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_vmem=5G

module load nextflow

nextflow -DXmx=5G -C motifDBComparison2_true.config run motifDBComparison2.nf
