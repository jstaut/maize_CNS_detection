
import pandas as pd
import sys

args = sys.argv
if len(args) != 5:
    print("wrong number of args")
    sys.exit()
else:
    clusterFile = args[1]
    resultFolder = args[2]
    onlyMaize = bool(args[3])
    minimalNbOfSpecies = int(args[4])

def countClusterSizes(clusterSizes):
    sizeCounts = dict()
    for i in range(1, max(clusterSizes)+1):
        sizeCounts[i] = (clusterSizes == i).sum()
    return sizeCounts

# Read in Orthogroup file
orthoGroups = pd.read_csv(clusterFile, sep='\t', header=0)
numberOfSpecies = orthoGroups.shape[1]-1 # first column is OG ID
if minimalNbOfSpecies > numberOfSpecies:
    sys.exit('The minimal number of species is larger than the total amount of species')
# Calculate the number of clusters with a given amount of species
clusterSizes = numberOfSpecies - orthoGroups.iloc[:,1:].isnull().sum(axis=1)
sizeCountDict = countClusterSizes(clusterSizes)

# Write cluster size info to QC file
with open(resultFolder+"/clusters_QC.txt", 'w') as outFile:
    outFile.write("All multi-species groups retained:\n"+("\n".join([str(items[0])+"\t"+str(items[1]) for items in sizeCountDict.items()])+"\nTotal:\t"+str(sum(sizeCountDict.values()))+"\nTotal with minimal "+str(minimalNbOfSpecies)+" species:\t"+str(sum(clusterSizes >= minimalNbOfSpecies))+"\n"))

# Write clusters to file
orthoGroups.to_csv(resultFolder+"/Orthogroups_all.tsv", sep='\t', index=False)

# Only retain clusters with at least a given amount of species
orthoGroups_reduced = orthoGroups.iloc[(clusterSizes >= minimalNbOfSpecies).values,:]
# Write reduced clusters to file
orthoGroups_reduced.to_csv(resultFolder+"/Orthogroups_minimal"+str(minimalNbOfSpecies)+"Species.tsv", sep='\t', index=False)

# Only keep groups with maize included
if onlyMaize:
    orthoGroups_onlyMaize = orthoGroups.iloc[(~orthoGroups["zma"].isnull()).values,:]
    clusterSizes = numberOfSpecies - orthoGroups_onlyMaize.iloc[:,1:].isnull().sum(axis=1)
    sizeCountDict = countClusterSizes(clusterSizes)
    
    # Write cluster size info to QC file
    with open(resultFolder+"/clusters_QC.txt", 'a') as outFile:
        outFile.write("\nOnly multi-species groups that include maize retained:\n"+("\n".join([str(items[0])+"\t"+str(items[1]) for items in sizeCountDict.items()])+"\nTotal:\t"+str(sum(sizeCountDict.values()))+"\nTotal with minimal "+str(minimalNbOfSpecies)+" species:\t"+str(sum(clusterSizes >= minimalNbOfSpecies))+"\n"))
    
    orthoGroups_onlyMaize.to_csv(resultFolder+"/Orthogroups_onlyMaize.tsv", sep='\t', index=False)

    # Only retain clusters with at least a given amount of species
    orthoGroups_onlyMaize_reduced = orthoGroups_onlyMaize.iloc[(clusterSizes >= minimalNbOfSpecies).values,:]
    # Write reduced clusters to file
    orthoGroups_onlyMaize_reduced.to_csv(resultFolder+"/Orthogroups_onlyMaize_minimal"+str(minimalNbOfSpecies)+"Species.tsv", sep='\t', index=False)
