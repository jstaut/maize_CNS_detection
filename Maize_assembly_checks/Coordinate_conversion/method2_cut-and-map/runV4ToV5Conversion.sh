#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_vmem=10G

module load nextflow

nextflow -DXmx=10G -C v4ToV5Conversion.config run v4ToV5Conversion.nf -resume
