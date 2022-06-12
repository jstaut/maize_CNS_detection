
import sys
import shutil
from Bio.Seq import Seq

args = sys.argv
if len(args) != 4:
    print("wrong number of args")
    sys.exit()
else:
    BEDFolderLoc = args[1]
    selectionFile = args[2]
    outFolderLoc = args[3]

# Try to copy the BED file with the selected motif name
# If not found, try to copy the file of the reverse complement motif
with open(selectionFile) as file:
    for line in file:
        motif = line.rstrip()
        try:
            shutil.copy(BEDFolderLoc+"/"+motif+".bed", outFolderLoc)
        except FileNotFoundError:
            motif_rc = str(Seq(motif).reverse_complement())
            shutil.copy(BEDFolderLoc+"/"+motif_rc+".bed", outFolderLoc)