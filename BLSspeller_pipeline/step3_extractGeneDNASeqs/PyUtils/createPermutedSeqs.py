
from random import sample, seed
from os import listdir
from os.path import isfile, join
import os
import sys

args = sys.argv
if len(args) != 2:
    print("wrong number of args")
    sys.exit()
else:
    loc_output = args[1]

os.chdir(loc_output)

# This script assumes the fasta's to be in 2 folders: "masked" and "unmasked"
# It also assumes file naming of "(un)masked_promoter_{species}.fasta
maskedFiles = [f for f in listdir("masked") if (isfile(join("masked", f)) and f[-6:]=='.fasta')]
species = list()
for maskedFile in maskedFiles:
    species.append(maskedFile.split("_")[-1].split(".")[0])

# This script will take the unmasked sequences and shuffle the bases, then it will mask the original positions as exons
for spec in species:
    seed(0)
    with open("masked/masked_promoter_"+spec+".fasta") as maskedFile:
        with open("unmasked/unmasked_promoter_"+spec+".fasta") as unmaskedFile:
            with open("permutedSeqs/permutedSeqs_promoter_"+spec+".fasta", 'w') as outFile:
                while True:
                    maskedLine = maskedFile.readline().rstrip()
                    unmaskedLine = unmaskedFile.readline().rstrip()
                    if not maskedLine:
                        break
                    if maskedLine[0] == '>':
                        outFile.write(maskedLine+"\n")
                    else:
                        shuffledSeqList = sample(unmaskedLine, len(unmaskedLine))
                        for i in range(len(maskedLine)):
                            if maskedLine[i] == "N":
                                shuffledSeqList[i] = "N"
                        outFile.write("".join(shuffledSeqList)+"\n")
