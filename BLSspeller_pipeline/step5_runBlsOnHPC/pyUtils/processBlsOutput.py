
import pandas as pd
import sys

args = sys.argv
if len(args) != 5:
    print("wrong number of args")
    sys.exit()
else:
    loc_motifs = args[1]
    loc_out_QC = args[2]
    loc_out_motifs = args[3]
    selectedConfCutoff = None if args[4] == 'None' else int(args[4])

loc_out_QC_filtered = loc_out_QC+"/confidence_filtered_QC.txt"
loc_out_QC = loc_out_QC+"/confidence_QC.txt"

# Read in data and separate motifs and confidence scores
table = pd.read_csv(loc_motifs, sep='\t', header=None)
motifs = table[0]
# Select confidence score columns as columns that only have values between 0 and 1
confidence = table.iloc[:,1:].loc[:,((table.iloc[:,1:]<1) & (table.iloc[:,1:]>0)).sum(axis=0) > 0]
# Isolate conserved family count columns
nbOfConfCols = confidence.shape[1]
conservedFamilyCount = table.loc[:,1:nbOfConfCols]

# Filter the motif table based on a minimal number of gene families: 25
filteredConf = confidence.loc[conservedFamilyCount.iloc[:,-1]>=25,:]
filteredMotifs = motifs[conservedFamilyCount.iloc[:,-1]>=25]

# Write the number of motifs for different confidence cutoffs into a QC file
with open(loc_out_QC, 'w') as outFile:
    outFile.write("Confidence score\tNumber of motifs\n")
for confCutoff in range(70, 100):
    with open(loc_out_QC, 'a') as outFile:
        outFile.write(str(int(confCutoff))+"\t"+str(sum((confidence>confCutoff/100).sum(axis=1)>0))+"\n")

# Write the number of motifs for the filtered confidence cutoffs into a QC file
with open(loc_out_QC_filtered, 'w') as outFile:
    outFile.write("Confidence score\tNumber of motifs\n")
for confCutoff in range(70, 100):
    with open(loc_out_QC_filtered, 'a') as outFile:
        outFile.write(str(int(confCutoff))+"\t"+str(sum((filteredConf>confCutoff/100).sum(axis=1)>0))+"\n")

# Given a confidence cutoff, write the selected motifs into a file
if selectedConfCutoff is not None:
    loc_out_file = loc_out_motifs+"/conf"+str(selectedConfCutoff)+"_selectedMotifs.txt"
    selectedMotifs = motifs[(confidence>selectedConfCutoff/100).sum(axis=1)>0]
    with open(loc_out_file, 'w') as outFile:
        for motif in selectedMotifs.values:
            outFile.write(motif+"\n")

# Given a confidence cutoff, write the filtered, selected motifs into a file
if selectedConfCutoff is not None:
    loc_out_file = loc_out_motifs+"/conf"+str(selectedConfCutoff)+"_filtered_selectedMotifs.txt"
    selectedMotifs = filteredMotifs[(filteredConf>selectedConfCutoff/100).sum(axis=1)>0]
    with open(loc_out_file, 'w') as outFile:
        for motif in selectedMotifs.values:
            outFile.write(motif+"\n")
