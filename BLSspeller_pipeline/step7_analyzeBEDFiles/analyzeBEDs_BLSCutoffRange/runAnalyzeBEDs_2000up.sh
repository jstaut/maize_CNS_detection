#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_vmem=5G

module load nextflow

nextflow -DXmx=5G -C analyzeBEDs_2000up.config run analyzeBEDs.nf
