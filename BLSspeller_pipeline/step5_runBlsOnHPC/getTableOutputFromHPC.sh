#!/bin/bash

# USAGE: ./getOutputFromHPC.sh <location of files on HPC> <target location>

HPC_LOC=$1
TARGET_LOC=$2

# Make target location with 'tableOutput' folder
mkdir -p $TARGET_LOC/tableOutput
# Write original HPC location to target folder
echo $HPC_LOC > $TARGET_LOC/tableLocationOnHPC.txt

# Retrieve files from HPC
scp -i /home/jasst/.ssh/id_rsa vsc44506@login.hpc.ugent.be:$HPC_LOC/* $TARGET_LOC/tableOutput
