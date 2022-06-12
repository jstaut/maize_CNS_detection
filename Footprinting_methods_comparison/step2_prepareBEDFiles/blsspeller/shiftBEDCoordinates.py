#!/usr/bin/env python3

import os
os.environ["OMP_NUM_THREADS"] = "1"
import sys
import pandas as pd

args = sys.argv
if len(args) != 3:
    print("wrong number of args")
    sys.exit()
else:
    BEDFile = args[1]
    outputName = args[2]

# Shift coordinates by 1
bed = pd.read_csv(BEDFile, sep='\t', header=None)
bed[1] = bed[1]+1
bed[2] = bed[2]+1
bed.to_csv(outputName, sep='\t', header=False, index=False)
