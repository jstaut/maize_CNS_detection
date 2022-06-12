#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_vmem=10G

module load nextflow

nextflow -DXmx=10G -C getOrthoGroups_Pac.config run getOrthoGroups.nf
