#!/bin/bash

# USAGE: ./getBEDOutputFromHPC.sh <location of files on HPC> <target location> <confidence score>

HPC_LOC=$1
TARGET_LOC=$2
CONF=$3

# Make target location with 'tableOutput' folder
mkdir -p $TARGET_LOC/BEDOutput_conf$CONF
# Write original HPC location to target folder
echo $HPC_LOC > $TARGET_LOC/BEDLocationOnHPC.txt

# Retrieve files from HPC
scp -i /home/jasst/.ssh/id_rsa vsc44506@login.hpc.ugent.be:$HPC_LOC/* $TARGET_LOC/BEDOutput_conf$CONF
