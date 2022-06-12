
oldGenome_ch = Channel
                .fromPath("$params.input/*_old.fasta", type: 'file' )

newGenome_ch = Channel
                .fromPath("$params.input/*_new.fasta", type: 'file' )

bed_ch = Channel
                .fromPath("$params.input/*.bed", type: 'file' )

process cutOutSequences {

    module 'bedtools/x86_64/2.2.28'

    input:
    file oldGenome from oldGenome_ch
    file bed from bed_ch

    output:
    file "sequences.fasta" into cutSequences_ch

    """
    bedtools getfasta -fi $oldGenome -bed $bed -fo sequences.fasta
    """
}

process mapToNewGenome {

    module 'last/x86_64/809'

    input:
    file newGenome from newGenome_ch
    file sequences from cutSequences_ch

    output:
    file "alignedSequences.maf" into maf_ch

    """
    lastdb -c -m1111110 newGenomeDB $newGenome
    lastal -P$params.lastThreads -q3 -e30 newGenomeDB $sequences > alignedSequences.maf
    """
}

process processMaf {

    module 'python/x86_64/3.6.5'

    input:
    file maf from maf_ch

    output:
    file "mafTable.txt" into mafTable_ch

    """
    #!/usr/bin/env python3

    loc_in = "$maf"

    names = []
    mafFormat1 = ["name1", "start1", "alnSize1", "strand1", "seqSize1", "alignment1"]
    mafFormat2 = ["name2", "start2", "alnSize2", "strand2", "seqSize2", "alignment2"]
    extra = ["perc_id", "cutStart", "cutEnd"]

    with open(loc_in) as inputFile:
        with open("mafTable.txt", 'a') as outFile:
            record = []
            firstSequence = True
            for line in inputFile:
                if line[0] == "a":
                    metrics = line.split()[1:]
                    if len(names) == 0:
                        for metric in metrics:
                            names.append(metric.split("=")[0])
                        names.extend(mafFormat1)
                        names.extend(mafFormat2)
                        names.extend(extra)
                        outFile.write("\\t".join(names) + "\\n")
                    for metric in metrics:
                        record.append(metric.split("=")[1])
                elif line[0] == "s" and firstSequence:
                    align = line.split()[1:]
                    for i in range(len(align)):
                        record.append(align[i])
                    firstSequence = False
                elif line[0] == "s":
                    align = line.split()[1:]
                    chrom, cutRange = align[0].split(":")
                    align[0] = chrom
                    cutBegin, cutEnd = cutRange.split("-")
                    for i in range(len(align)):
                        record.append(align[i])
                    mismatch = 0
                    al1 = record[names.index("alignment1")]
                    al2 = record[names.index("alignment2")]
                    for i in range(len(al1)):
                        if al1[i] != al2[i]:
                            mismatch += 1
                    record.append(str(100*((len(al1)-mismatch)/float(len(al1)))))
                    record.append(cutBegin)
                    record.append(cutEnd)
                    outFile.write("\\t".join(record) + "\\n")
                    record = []
                    firstSequence = True
    """
}

process filterAndCastIntoBed {

    module 'R/x86_64/3.5.1'

    input:
    file processedMaf from mafTable_ch

    output:
    file "filtered.bed" into originalBed_ch

    """
    #!/usr/bin/env Rscript

    mafTable <- read.table("$processedMaf", header = T, sep = "\\t")
    mafTable_filter <- mafTable[mafTable\$perc_id >= 99,]
    BED <- data.frame(name=mafTable_filter\$name1, start=mafTable_filter\$start1, end=mafTable_filter\$start1+mafTable_filter\$alnSize1, strand=mafTable_filter\$strand1)
    write.table(BED, file = "filtered.bed", sep = "\\t", quote = F, row.names = F, col.names = F)
    """
}

process mergeBedFiles {

    module 'bedtools/x86_64/2.2.28'

    publishDir "$params.outdir", mode: 'copy'

    input:
    file bed from originalBed_ch

    output:
    file "merged.bed" into mergedBed_ch

    """
    bedtools sort -i $bed | mergeBed > merged.bed
    """
}
