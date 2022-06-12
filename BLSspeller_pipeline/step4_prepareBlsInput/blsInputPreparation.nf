
geneDNASequences_ch = Channel.fromPath("$params.input/*.fasta", type: 'file')

Channel.fromPath("$params.input/*.tsv", type: 'file')
       .into { orthogroups1_ch; orthogroups2_ch }

tree_ch = Channel.fromPath("$params.input/*.txt", type: 'file' )

process addClusterInfo {

    module 'python/x86_64/3.6.5'

    input:
    tuple file(geneDNASequence), file(orthogroups) from geneDNASequences_ch.combine(orthogroups1_ch)

    output:
    file "*.fasta" into seqsWithClusterInfo_ch

    """
    #!/usr/bin/env python3

    import pandas as pd
    from Bio import SeqIO

    # Get clusters
    clusters = pd.read_csv("$orthogroups", sep='\\t')

    # Get organism name, based on file name structure
    organism = ("_".join("$geneDNASequence".split("/")[-1].split("_")[2:])).split(".")[0]

    # Write promoters with OG info to file
    with open(organism+"_geneDNAWithCluster.fasta", 'a') as outFile:
        if organism in clusters.columns[1:]:
            promDict = dict() # Get promoter sequences
            for prom in SeqIO.parse(open("$geneDNASequence"),'fasta'):
                promDict[prom.id.split(':')[0]] = [prom.id, str(prom.seq)] # Gene ID as key, full header and sequence as value
            empty = clusters[organism].isnull()
            for i in range(len(clusters.index)):
                if not empty[i]:
                    clusterGenes = clusters[organism][i].split(', ')
                    for gene in clusterGenes:
                        info = promDict[gene]
                        header = info[0].split(":") # Split the header
                        header[1] = clusters["Orthogroup"][i] # Add the OG ID as second field in the header
                        outFile.write(">"+(":".join(header))+"\\n"+info[1]+"\\n")
    """
}

process makeOrthoFiles {

    module 'python/x86_64/3.6.5'

    input:
    file orthogroups from orthogroups2_ch
    file tree from tree_ch
    file geneDNAWithCluster from seqsWithClusterInfo_ch.collect()

    output:
    file "*.txt" into blsInputFiles_ch

    """
    #!/usr/bin/env python3

    import os
    os.environ["OMP_NUM_THREADS"] = "1"
    import pandas as pd
    from Bio import SeqIO
    from os import listdir, path
    import importlib.machinery as imp_mach
    import importlib.util as imp_util
    from random import seed, shuffle

    # Import own module
    loader = imp_mach.SourceFileLoader('treeModule', '$params.treeUtils/treeModule.py')
    spec = imp_util.spec_from_loader('treeModule', loader)
    treeUtils = imp_util.module_from_spec(spec)
    loader.exec_module(treeUtils)

    # A function that takes a list of lists, containing genes per orthogroup and shuffles the IDs while keeping the structure
    def shuffleGenes(genesPerCluster):
        seed(1)
        flattened = [gene for cluster in genesPerCluster for gene in cluster]
        shuffle(flattened)
        for i in range(len(genesPerCluster)):
            for j in range(len(genesPerCluster[i])):
                genesPerCluster[i][j] = flattened.pop(0) # replace the original elements by the shuffled list
        return genesPerCluster
    
    # Load in tree + clean it up (rescale full tree to 1 & remove child labels)
    with open("$tree") as treeFile:
        tree = treeFile.readline()
    tree = treeUtils.cleanTree(tree)

    # Get clusters
    clusters = pd.read_csv("$orthogroups", sep='\\t')
    clusterSpecies = list(clusters.columns[1:])
    
    # Get number of species in each cluster
    nbOfSpecies = len(clusterSpecies)
    nbOfClusterSpecies_all = (nbOfSpecies - clusters.iloc[:,1:].isnull().sum(axis=1)).values
    nbOfClusterSpecies = dict() # Keys are the OG ID, values are the number of species in the OG
    for i in range(len(clusters.index)):
        nbOfClusterSpecies[clusters["Orthogroup"][i]] = str(nbOfClusterSpecies_all[i])
    
    # Write the cluster ID, tree and number of species into orthogroup files
    for orthogroup, nbOfSpec in nbOfClusterSpecies.items():
        with open(orthogroup+".txt", 'a') as orthoFile:
            orthoFile.write(orthogroup+"\\n"+tree+"\\n"+nbOfSpec+"\\n")
        if $params.createFalseOrthoGroups:
            with open(orthogroup+"_fake.txt", 'a') as orthoFile:
                orthoFile.write(orthogroup+"\\n"+tree+"\\n"+nbOfSpec+"\\n")
    
    # Get promotor file info (sequences)
    promFiles = dict() # Key is a species name that can be found in the cluster file, value is the file name of the promoter sequences
    for fileName in "$geneDNAWithCluster".split(" "):
        if fileName[-6:].lower() == ".fasta" and path.getsize(fileName)>0: # read non-empty fasta files
            species = fileName.split("/")[-1].replace("_geneDNAWithCluster.fasta", "").split(".")[0]
            if species in clusterSpecies:
                promFiles[species] = fileName
            else:
                raise Exception("Mismatch between species in promoter files (= all FASTA files!) vs orthogroup file!\\n"+species+" not present in orthogroup file")
    
    for species, promFile in promFiles.items():
        promDict = dict() # Keys are gene IDs, values are promoter sequences
        
        if $params.createFalseOrthoGroups:
            genesPerCluster = list()
        
        # Store all promoter sequences of a given species in a dict
        for prom in SeqIO.parse(open(promFile),'fasta'):
            promDict[prom.id.split(':')[0]] = str(prom.seq) # Keys are gene IDs, values are promoter sequences
        empty = clusters[species].isnull()
        for i in range(len(clusters.index)):
            if not empty[i]:
                orthogroup = clusters["Orthogroup"][i]
                clusterGenes = clusters[species][i].split(', ')
                
                if $params.createFalseOrthoGroups:
                    genesPerCluster.append(clusterGenes) # create list of lists, containing gene IDs per cluster
                
                promSeqs = [] # A list with the promoter sequences, in the same order as 'clusterGenes'
                for gene in clusterGenes:
                    promSeqs.append(promDict[gene])
                with open(orthogroup+".txt", 'a') as orthoFile:
                    orthoFile.write(" ".join(clusterGenes)+"\\t"+species+"\\n"+" ".join(promSeqs)+"\\n")
        
        if $params.createFalseOrthoGroups:
            fakeClusters = shuffleGenes(genesPerCluster)
            for i in range(len(clusters.index)):
                if not empty[i]:
                    orthogroup = clusters["Orthogroup"][i]
                    clusterGenes = fakeClusters.pop(0) # The first list of genes
                    promSeqs = []
                    for gene in clusterGenes:
                        promSeqs.append(promDict[gene])
                    with open(orthogroup+"_fake.txt", 'a') as orthoFile:
                        orthoFile.write(" ".join(clusterGenes)+"\\t"+species+"\\n"+" ".join(promSeqs)+"\\n")
    """
}

process sortAndTarFiles {

    publishDir "$params.outdir", mode: 'move'

    input:
    file blsInputFiles from blsInputFiles_ch.collect()

    output:
    file "*.tar.gz" into final_ch

    """
    mkdir input_true
    if [ $params.createFalseOrthoGroups == True ]
    then
        mkdir input_fake
    fi
    for BLS_INPUT_FILE in $blsInputFiles
    do
        if [ \${BLS_INPUT_FILE:(-9)} == "_fake.txt" ]
        then
            cp \$BLS_INPUT_FILE input_fake
        else
            cp \$BLS_INPUT_FILE input_true
        fi
    done

    tar -czf input_true.tar.gz input_true
    if [ $params.createFalseOrthoGroups == True ]
    then
        tar -czf input_fake.tar.gz input_fake
    fi
    
    """
}
