
motifTableParts_ch = Channel
                        .fromPath("$params.input/part-*", type: 'file')
                        .collect()

Channel
    .fromPath("$params.precomputed/precomputedMotifs.pickle", type: 'file')
    .into{ precomputedMotifs1_ch; precomputedMotifs2_ch }

process getNbOfMotifsInDB {

    output:
    env DB_MOTIF_NUMBER into nbOfMotifsInDB_ch

    """
    DB_MOTIF_NUMBER=\$(grep '//' $params.motifDB | wc -l)
    """
}

process createFullTable {

    input:
    file motifTableParts from motifTableParts_ch

    output:
    file "motifTable.txt" into motifTable_ch

    """
    for MOTIF_TABLE_PART in $motifTableParts
    do
        cat \$MOTIF_TABLE_PART >> motifTable.txt
    done
    """
}

process selectMotifs {

    module 'python/x86_64/3.6.5'

    input:
    tuple file(motifTable), file(precomputedMotifs) from motifTable_ch.combine(precomputedMotifs1_ch)

    output:
    tuple file("motifsToBeCalculated.pickle"), file("motifSets.pickle") into pythonMotifs_ch

    """
    #!/usr/bin/env python3
    
    import os
    os.environ["OMP_NUM_THREADS"] = "1"
    import pandas as pd
    import pickle
    from random import shuffle, seed
    from difflib import SequenceMatcher
    from shutil import copy

    def shuffleMotif(motif):
        similarity = 1
        leastSimilarShuffle = motif
        for i in range(10):
            seed(i)
            charList = list(motif)
            shuffle(charList)
            shuffledMotif = ''.join(charList)
            newSimilarity = SequenceMatcher(None, motif, shuffledMotif).ratio()
            if newSimilarity < similarity:
                similarity = newSimilarity
                leastSimilarShuffle = shuffledMotif
        return leastSimilarShuffle
    
    # Read in input
    table = pd.read_csv("$motifTable", sep='\\t', header=None)
    blsThresholds = "$params.blsThresholds".split(",")
    confRange = [int(i) for i in "$params.confRange".split(",")]
    
    # Separate motifs and confidence scores
    motifs = table[0]
    # Select confidence score columns as columns that only have values between 0 and 1
    confidence = table.iloc[:,1:].loc[:,((table.iloc[:,1:]<1) & (table.iloc[:,1:]>0)).sum(axis=0) > 0]
    nbOfConfCols = confidence.shape[1]
    # Select family count columns score columns
    conservedFamilyCount = table.loc[:,1:nbOfConfCols]

    # Only keep the highest BLS threshold columns (as given by the BLS thresholds parameter)
    confidence = confidence.iloc[:,-len(blsThresholds):]
    conservedFamilyCount = conservedFamilyCount.iloc[:,-len(blsThresholds):]
    confidence.columns = conservedFamilyCount.columns
    nbOfConfCols = confidence.shape[1]
    assert nbOfConfCols == len(blsThresholds)

    # Write the number of motifs for different confidence cutoffs into a QC file
    with open("motifsConf_QC.txt", 'w') as outFile:
        outFile.write("Confidence score\\tNumber of motifs\\n")
        for confCutoff in range(confRange[0], confRange[1]):
            outFile.write(str(int(confCutoff))+"\\t"+str(sum((confidence>confCutoff/100).sum(axis=1)>0))+"\\n")
    
    # Same as above, but filtered based on minimal family count filter
    with open("motifsConfFiltered_QC.txt", 'w') as outFile:
        outFile.write("Confidence score\\tNumber of motifs\\n")
        for confCutoff in range(confRange[0], confRange[1]):
            outFile.write(str(int(confCutoff))+"\\t"+str(sum(((confidence>confCutoff/100) & (conservedFamilyCount>=int($params.famCountThreshold))).sum(axis=1)>0))+"\\n")
    
    # Select sets of motifs for different threshold values
    motifSets = dict() # Dictionary with (bls,conf) thresholds as key and a set of motifs as value
    motifSets_null = dict() # Dictionary with (bls,conf) thresholds as key and a set of (shuffled) null motifs as value
    if blsThresholds[-1] == 0:
        # Only based on confidence, not taking BLS into account
        for confCutoff in range(confRange[0], confRange[1]):
            selectedMotifs = motifs[((confidence>confCutoff/100) & (conservedFamilyCount>=int($params.famCountThreshold))).sum(axis=1)>0]
            motifSet = set(selectedMotifs)
            motifSet_null = set()
            for motif in motifSet:
                motifSet_null.add(shuffleMotif(motif))
            motifSets[(0,confCutoff)] = motifSet
            motifSets_null[(0,confCutoff)] = motifSet_null
    else:
        # For multiple specific BLS and confidence thresholds
        with open("motifsConfAndBLS_QC.txt", 'w') as outFile:
            outFile.write("\\t"+("\\t".join([str(int(i)) for i in range(confRange[0], confRange[1])]))+"\\n")
            for bls_i in range(1, len(blsThresholds)+1):
                QCLine = [blsThresholds[-bls_i]]
                for confCutoff in range(confRange[0], confRange[1]):
                    selectedMotifs = motifs[(confidence.iloc[:,-bls_i]>(confCutoff/100)) & (conservedFamilyCount.iloc[:,-bls_i]>=int($params.famCountThreshold))]
                    motifSet = set(selectedMotifs)
                    motifSet_null = set()
                    for motif in motifSet:
                        motifSet_null.add(shuffleMotif(motif))
                    motifSets[(blsThresholds[-bls_i],confCutoff)] = motifSet
                    motifSets_null[(blsThresholds[-bls_i],confCutoff)] = motifSet_null
                    QCLine.append(str(len(motifSets[(blsThresholds[-bls_i],confCutoff)])))
                outFile.write(("\\t".join(QCLine))+"\\n")

    # Make a set with all selected motifs
    totalSelectedMotifs = set()
    for motifSet in motifSets.values():
        totalSelectedMotifs.update(motifSet)
    for motifSet in motifSets_null.values():
        totalSelectedMotifs.update(motifSet)
    
    # Load a python dictionary with motifs that already have their Ncor precomputed
    with open("$precomputedMotifs", 'rb') as file:
        precomputedDict = pickle.load(file)
    # Get set of precomputed motifs
    precomputedSet = set(precomputedDict.keys())
    
    # Get the motifs that still need Ncor calculation and write them into a file
    motifsToBeCalculated = totalSelectedMotifs - precomputedSet
    with open("motifsToBeCalculated.pickle", 'wb') as outFile:
        pickle.dump(motifsToBeCalculated, outFile)

    # Save dictionaries of motif sets for both the real and null motifs
    allMotifSets = [motifSets, motifSets_null]
    with open("motifSets.pickle", 'wb') as outFile:
        pickle.dump([motifSets, motifSets_null], outFile)
    
    # Copy QC files to output directory
    copy("motifsConf_QC.txt", "$params.outdir")
    copy("motifsConfFiltered_QC.txt", "$params.outdir")
    copy("motifsConfAndBLS_QC.txt", "$params.outdir")
    """
}

pythonMotifs_ch
            .multiMap { it ->
            toBeCalculated: it[0]
            sets: it[1]
            }
            .set { motifs_ch }

process motifsToTransfac {

    module 'python/x86_64/3.5.1'

    input:
    file motifsToBeCalculated from motifs_ch.toBeCalculated

    output:
    file "*_transfac_chunk*" into transfacMotifs_ch

    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import numpy as np
    import pickle
    from Bio.Seq import Seq
    from Bio import Seq as Seqq
    from itertools import product
    from Bio import motifs as motifsBioython

    def extend_ambiguous_dna(seq):
        d = Seqq.IUPAC.IUPACData.ambiguous_dna_values
        return list(map("".join, product(*map(d.get, seq))))
    
    # Get the motifs to convert
    with open("$motifsToBeCalculated", 'rb') as file:
        motifSet = pickle.load(file)
    
    motifs = list(motifSet)
    chunkSize = $params.chunkSize
    chunkNb = 0
    while chunkNb < len(motifSet)/chunkSize or chunkNb == 0:
        fileName_motifs = 'motifs_transfac_chunk'+''.join(['0']*(3-len(str(chunkNb))))+str(chunkNb)
        motifSubset = motifs[(chunkNb*chunkSize):(chunkNb*chunkSize+chunkSize)]
        with open(fileName_motifs, 'w') as f:
            for motif in motifSubset:
                seqs = []
                for j in extend_ambiguous_dna(motif):
                    seqs.append(Seq(j))
                m = motifsBioython.create(seqs)
                m.name = motif
                f.write("AC "+motif+"\\nXX\\n")
                f.write(m.format('transfac'))
        chunkNb += 1
    """
}

process comparePFMs {

    module "rsa-tools"

    input:
    file motifChunk from transfacMotifs_ch.flatten()

    output:
    file "${motifChunk}_cors.tab" into newNcor_ch

    """
    
    compare-matrices -file1 $motifChunk -file2 $params.motifDB -format transfac -o ${motifChunk}_cors -lth Ncor 0

    """
}

process updateNcorsAndCalculateMetrics {

    module 'python/x86_64/3.6.5'

    publishDir "$params.outdir", mode: 'copy'

    input:
    file inputFiles from newNcor_ch.collect().combine(precomputedMotifs2_ch).combine(motifs_ch.sets).combine(nbOfMotifsInDB_ch)

    output:
    file "*.txt" into final_ch

    """
    #!/usr/bin/env python3

    import os
    os.environ["OMP_NUM_THREADS"] = "1"
    import pandas as pd
    import numpy as np
    import pickle

    def getNullMotifFDR(ncorList, ncorList_null, ncorCutoff):
        positives = sum([i>ncorCutoff for i in ncorList])
        if len(ncorList) == 0 or len(ncorList_null) == 0 or positives == 0:
            return None
        else:
            falsePositives = sum([i>ncorCutoff for i in ncorList_null])
            positives_norm = positives/len(ncorList) # Normalize based on total amount of motifs
            falsePositives_norm = falsePositives/len(ncorList_null) # Normalize based on total amount of motifs
            return falsePositives_norm/positives_norm # Calculate the FDR

    def getNbOfMatches(ncorList, ncorList_null, ncorCutoff, estimated=False):
        positives = sum([i>ncorCutoff for i in ncorList])
        if estimated:
            nullMotifFDR = getNullMotifFDR(ncorList, ncorList_null, ncorCutoff)
            # Calculate the estimated number of true matches by substracting the FDR
            return None if nullMotifFDR is None else (1-nullMotifFDR)*positives 
        else:
            return positives

    def getNbOfMatchedDBMotifs(motifSet, ncorDict, ncorCutoff, familyInfo = None):
        # Retrieve database hits (motifs in DB that are matched by a BLS motif)
        databaseHits = set()
        for motif in motifSet:
            for i in range(len(ncorDict[motif][0])):
                if ncorDict[motif][0][i] > ncorCutoff:
                    DBMotifName = ncorDict[motif][1][i]
                    if familyInfo is None:
                        databaseHits.add(DBMotifName)
                    else:
                        DBMotifFamily = familyInfo[DBMotifName]
                        databaseHits.add(DBMotifFamily)
        return len(databaseHits)

    def getNonEstimatedMetrics(motifSet, ncorDict, ncorList, ncorList_null, dbSize, ncorCutoff, familyInfo):
        
        nbOfMotifs = len(motifSet)
        
        M1_nbOfMatches = getNbOfMatches(ncorList, ncorList_null, ncorCutoff)
        M2_precision = None if nbOfMotifs == 0 else round(100*M1_nbOfMatches/nbOfMotifs,2)
        M3_nullMotifFDR = getNullMotifFDR(ncorList, ncorList_null, ncorCutoff)
        M3_nullMotifFDR = None if M3_nullMotifFDR is None else round(100*M3_nullMotifFDR,2)
        M4_nbOfMatchedDBMotifs = getNbOfMatchedDBMotifs(motifSet, ncorDict, ncorCutoff)
        M5_recall = round(100*M4_nbOfMatchedDBMotifs/dbSize,2)
        M6_F1 = None if nbOfMotifs == 0 else (2*(M2_precision/100)*(M5_recall/100))/((M2_precision/100)+(M5_recall/100))
        M6_F1 = None if M6_F1 is None else round(100*M6_F1,2)
        M7_nbOfMatchedDBFamilies = getNbOfMatchedDBMotifs(motifSet, ncorDict, ncorCutoff, familyInfo=familyInfo)
        M8_familyRecall = round(100*M7_nbOfMatchedDBFamilies/(len(set(familyInfo.values()))),2)
        M9_familyF1 = None if nbOfMotifs == 0 else (2*(M2_precision/100)*(M8_familyRecall/100))/((M2_precision/100)+(M8_familyRecall/100))
        M9_familyF1 = None if M9_familyF1 is None else round(100*M9_familyF1,2)
        
        naming = ["nbOfMatches", "precision", "nullMotifFDR", "nbOfMatchedDBMotifs", "recall", "F1", "nbOfMatchedDBFamilies", "familyRecall", "familyF1"]
        naming = [metricName+"At"+str(round(100*ncorCutoff)) for metricName in naming]
        
        return naming, [M1_nbOfMatches, M2_precision, M3_nullMotifFDR, M4_nbOfMatchedDBMotifs, M5_recall, M6_F1, M7_nbOfMatchedDBFamilies, M8_familyRecall, M9_familyF1]
        

    def calculateMetrics(motifSet, motifSet_null, ncorDict, dbSize, familyInfo):
        
        # Retrieve Ncors
        ncorList = list()
        ncorList_null = list()
        for motif in motifSet:
            ncorList.extend(ncorDict[motif][0])
        for motif in motifSet_null:
            ncorList_null.extend(ncorDict[motif][0])
        
        M1_nbOfMotifs = len(motifSet)
        
        M2_1_estNbOfMatches = getNbOfMatches(ncorList, ncorList_null, 0.75, estimated=True)
        M2_1_estNbOfMatches = None if M2_1_estNbOfMatches is None else round(M2_1_estNbOfMatches)
        M2_2_estPrecision = None if (M1_nbOfMotifs == 0 or M2_1_estNbOfMatches is None) else round(100*M2_1_estNbOfMatches/M1_nbOfMotifs,2)
        realNbMatched = getNbOfMatchedDBMotifs(motifSet, ncorDict, 0.75)
        nullNbMatched = getNbOfMatchedDBMotifs(motifSet_null, ncorDict, 0.75)
        M2_3_estNbOfMatchedDBMotifs = realNbMatched-nullNbMatched
        M2_4_estRecall = round(100*M2_3_estNbOfMatchedDBMotifs/dbSize,2)
        M2_5_estF1 = None if M2_2_estPrecision is None else (2*(M2_2_estPrecision/100)*(M2_4_estRecall/100))/((M2_2_estPrecision/100)+(M2_4_estRecall/100))
        M2_5_estF1 = None if M2_5_estF1 is None else round(100*M2_5_estF1,2)
        # Calculate the TF family based metrics
        realNbMatched = getNbOfMatchedDBMotifs(motifSet, ncorDict, 0.75, familyInfo=familyInfo)
        nullNbMatched = getNbOfMatchedDBMotifs(motifSet_null, ncorDict, 0.75, familyInfo=familyInfo)
        M2_6_estNbOfMatchedDBMotifs = realNbMatched-nullNbMatched
        M2_7_estFamilyRecall = round(100*M2_6_estNbOfMatchedDBMotifs/(len(set(familyInfo.values()))),2)
        M2_8_estFamilyF1 = None if M2_2_estPrecision is None else (2*(M2_2_estPrecision/100)*(M2_7_estFamilyRecall/100))/((M2_2_estPrecision/100)+(M2_7_estFamilyRecall/100))
        M2_8_estFamilyF1 = None if M2_8_estFamilyF1 is None else round(100*M2_8_estFamilyF1,2)
        
        namingAt80, nonEstMetricsAt80 = getNonEstimatedMetrics(motifSet, ncorDict, ncorList, ncorList_null, dbSize, 0.80, familyInfo=familyInfo)
        namingAt85, nonEstMetricsAt85 = getNonEstimatedMetrics(motifSet, ncorDict, ncorList, ncorList_null, dbSize, 0.85, familyInfo=familyInfo)
        
        naming = ["nbOfMotifs", "estimatedNbOfMatches", "estimatedPrecision", "estimatedNbOfDBMatches", "estimatedRecall", "estimatedF1", "estimatedFamilyNbOfDBMatches", "estimatedFamilyRecall", "estimatedFamilyF1"]
        naming.extend(namingAt80)
        naming.extend(namingAt85)
        
        allMetrics = [M1_nbOfMotifs,M2_1_estNbOfMatches,M2_2_estPrecision,M2_3_estNbOfMatchedDBMotifs,M2_4_estRecall,M2_5_estF1,M2_6_estNbOfMatchedDBMotifs,M2_7_estFamilyRecall,M2_8_estFamilyF1]
        allMetrics.extend(nonEstMetricsAt80)
        allMetrics.extend(nonEstMetricsAt85)
        
        allStringMetrics = [str(i) for i in allMetrics]
        
        return naming, allStringMetrics
    
    # Separate input files
    allInputFiles = "$inputFiles".split(" ")
    DBNbFile = allInputFiles.pop(-1)
    with open(DBNbFile) as file:
        nbOfmotifsInDB = int(file.readline().rstrip())
    motifSetsFile = allInputFiles.pop(-1)
    precomputedMotifsFile = allInputFiles.pop(-1)
    ncorFiles = allInputFiles # what remains are Ncor files

    # Get the TF family info
    with open("$params.familyDict", 'rb') as file:
        familyDict = pickle.load(file)

    # Get the precomputed motif-Ncor dictionary
    with open(precomputedMotifsFile, 'rb') as file:
        motifsWithNcor = pickle.load(file)
    precomputedMotifs = set(motifsWithNcor.keys())
    
    # Update the precomputed dictionary with the new motifs
    for ncorFile in ncorFiles:
        ncorTable = pd.read_csv(ncorFile, sep='\\t')
        for i in range(ncorTable.shape[0]):
            if not ncorTable.iloc[i,0] in precomputedMotifs: # make sure to only add info on a motif if it has not been computed before
                if not ncorTable.iloc[i,0] in motifsWithNcor:
                    motifsWithNcor[ncorTable.iloc[i,0]] = ([ncorTable.iloc[i,2]], [ncorTable.iloc[i,1]]) # store all matches in a list
                else:
                    motifsWithNcor[ncorTable.iloc[i,0]][0].append(ncorTable.iloc[i,2]) # add extra matches to list
                    motifsWithNcor[ncorTable.iloc[i,0]][1].append(ncorTable.iloc[i,1])
    
    # Get the motif sets dictionaries
    with open(motifSetsFile, 'rb') as file:
        allMotifSets = pickle.load(file)
    motifSets, motifSets_null = allMotifSets[0], allMotifSets[1]

    # Very rarely (4 in 100,000 motifs), a motif does not get an Ncor assigned to it by matrix-comparison
    # Remove these motifs from the motif sets, as they cannot be used for Ncor distribution calculations
    ncorMotifs = set(motifsWithNcor.keys())
    for thresholds, motifset in motifSets.items():
        motifSets[thresholds] = motifset & ncorMotifs
    for thresholds, motifset in motifSets_null.items():
        motifSets_null[thresholds] = motifset & ncorMotifs

    # Calculate all needed metrics
    allMetricNames = None
    allCalculatedMetrics =  dict()
    for thresholds in motifSets.keys():
        naming, metrics = calculateMetrics(motifSets[thresholds], motifSets_null[thresholds], motifsWithNcor, nbOfmotifsInDB, familyDict)
        allCalculatedMetrics[thresholds] = metrics
        if allMetricNames is None:
            allMetricNames = naming

    blsThresholds = "$params.blsThresholds".split(",")
    confRange = [int(i) for i in "$params.confRange".split(",")]
    confThresholds = [i for i in range(confRange[0],confRange[1])]

    for i in range(len(allMetricNames)):
        with open("comparisonMetrics_"+allMetricNames[i]+".txt", 'w') as outFile:
            # Header
            outFile.write("\\t"+"\\t".join(blsThresholds)+"\\n")
            for confThreshold in confThresholds:
                confLine = list()
                for blsThreshold in blsThresholds:
                    thresholds = (blsThreshold, confThreshold)
                    confLine.append(allCalculatedMetrics[thresholds][i])
                outFile.write(str(confThreshold)+"\\t"+"\\t".join(confLine)+"\\n")

    # Overwrite precomputed motifs
    if ncorTable.shape[0] > 0:
        with open("$params.precomputed/precomputedMotifs.pickle", 'wb') as outFile:
            pickle.dump(motifsWithNcor, outFile)
    """
}
