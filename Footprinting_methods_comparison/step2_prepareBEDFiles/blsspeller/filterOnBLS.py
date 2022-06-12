#!/usr/bin/env python3

from os import listdir
from os.path import isfile, join
import sys

args = sys.argv
if len(args) != 4:
    print("wrong number of args")
    sys.exit()
else:
    BEDFolder = args[1]
    outputFile = args[2]
    threshold = int(args[3])

motifFiles = [join(BEDFolder, f) for f in listdir(BEDFolder) if (isfile(join(BEDFolder, f)) and f[-4:]==".bed")]

with open(outputFile, 'w') as outFile:
	for motifFile in motifFiles:
		targetGenes = set()
		with open(motifFile) as bedFile:
			for line in bedFile:
				fields = line.split('\t')
				if 100*float(fields[6]) >= threshold:
					targetGenes.add(fields[4])
					outFile.write(line)
