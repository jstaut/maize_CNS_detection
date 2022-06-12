#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_vmem=5G
module load nextflow
nextflow -DXmx=5G -C geneDNAExtraction.config run geneDNAExtraction.nf -resume silly_bell
cd temp_results/5000up1000down
mkdir masked
mkdir unmasked
find * -maxdepth 0 -type f -name 'masked_*' -exec sh -c 'mv "$0" masked' {} \;
find * -maxdepth 0 -type f -name 'unmasked_*' -exec sh -c 'mv "$0" unmasked' {} \;
