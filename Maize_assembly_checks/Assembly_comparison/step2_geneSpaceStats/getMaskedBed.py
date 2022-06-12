#!/usr/bin/env python3

# USAGE: getMaskedBed.py <fasta file>

import sys
import os

if len(sys.argv) == 2:

    if sys.argv[1].endswith(".fa") or sys.argv[1].endswith(".txt"):
        input_fasta = open(sys.argv[1],'r')
    else:
        raise Exception("Unsupported File Type")
else:
    print("Wrong number of arguments")
    raise SystemExit


n, state = 0, 0 # state: 0 is not masked, 1 is masked
chrom, start, end = None, None, None

with input_fasta as f:
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            # Print end range
            if state == 1:
                print('\t'.join([chrom ,str(start), str(n)]))
                start, end, state  = 0, 0, 0
            n = 0 # Reset character
            chrom = line.split(" ")[0].replace(">","")
        else:
            for char in line:
                if state == 0 and char == "N":
                    state = 1
                    start = n
                elif state == 1 and char != "N":
                    state = 0
                    end = n
                    print('\t'.join([chrom ,str(start), str(end)]))
                else:
                    pass

                n += 1 # First base is 0 in bed format.

# Print mask close if on the last chromosome.
if state == 1:
            print('\t'.join([chrom ,str(start), str(n)]))
            start, end, state  = 0, 0, 0