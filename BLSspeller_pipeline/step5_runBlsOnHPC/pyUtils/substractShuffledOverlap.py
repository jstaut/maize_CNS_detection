
import sys

args = sys.argv
if len(args) != 3:
    print("wrong number of args")
    sys.exit()
else:
    loc_true = args[1]
    loc_shuffled = args[2]

# Read in the selected motifs as a set (both for the true and shuffled data set)
trueMotifs = set()
shuffledMotifs = set()
with open(loc_true) as file:
    for motif in file:
        trueMotifs.add(motif.rstrip())
with open(loc_shuffled) as file:
    for motif in file:
        shuffledMotifs.add(motif.rstrip())

# Separate common and unique motifs between the true and shuffled motif sets
commonMotifs = trueMotifs & shuffledMotifs
uniquelyTrueMotifs = trueMotifs - shuffledMotifs
uniquelyShuffledMotifs = shuffledMotifs - trueMotifs

# Write the separated motifs to separate files
with open(loc_true[:-4]+"_uniquelyTrue.txt", 'w') as outFile:
    for motif in uniquelyTrueMotifs:
        outFile.write(motif+"\n")
with open(loc_shuffled[:-4]+"_uniquelyShuffled.txt", 'w') as outFile:
    for motif in uniquelyShuffledMotifs:
        outFile.write(motif+"\n")
with open(loc_true[:-4]+"_commonMotifs.txt", 'w') as outFile:
    for motif in commonMotifs:
        outFile.write(motif+"\n")
