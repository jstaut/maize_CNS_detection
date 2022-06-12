
genomeHTML_ch = Channel
                .fromPath("$params.input/*_gen.fa", type: 'file' )

library_ch = Channel
                .fromPath("$params.input/*Lib.fa", type: 'file' )

process getGenomeFromHTML {

    input:
    file genome from genomeHTML_ch

    output:
    tuple env(VERSION), file(genome) into genome_ch

    """
    VERSION=`echo $genome | grep -oP '(?<=v)[45](?=_gen.fa)' | awk '{print "V" \$1}'`
    """
}

genome_library_ch = genome_ch
                            .combine(library_ch)

process runRepeatMasker {

    module 'python/x86_64/3.6.5'
    module 'RepeatMasker/x86_64/4-1-1'
    module 'rmblastn/x86_64/2.10.0'

    publishDir "$params.outdir", mode: 'copy'

    input:
    tuple val(version), file(genome), file(library) from genome_library_ch

    output:
    tuple val(version), file("*.masked") into masked_ch

    """
    RepeatMasker -pa $params.maskerThreads -dir . -q -div 40 -lib $library -e rmblast -cutoff 225 -gff $genome
    """
}

// Duplicate and reverse the order of the two genomes, using the first value as the reference label

(genomes1_ch, genomes2_ch) = masked_ch.into(2)

bothGenomes1_ch = genomes1_ch
                            .flatten()
                            .toList()
                            .map{ it.subList(0,2) + [ it.get(3) ] }

bothGenomes2_ch = genomes2_ch
                            .flatten()
                            .toList()
                            .map{ it.subList(2,4) + [ it.get(1) ] }

maskedGenomes_ch = bothGenomes1_ch.mix(bothGenomes2_ch)

(maskedGenomes1_ch, maskedGenomes2_ch) = maskedGenomes_ch.into(2)

process runLAST {

    module 'last/x86_64/809'

    input:
    tuple val(version), file(genome1), file(genome2) from maskedGenomes1_ch

    output:
    tuple val(version), file("aligned.maf") into originalMaf_ch

    """
    lastdb -c -m1111110 maizeVersionDB $genome1
    lastal -P$params.lastThreads -q3 -e30 maizeVersionDB $genome2 | last-split > aligned.maf
    """
}

process processMaf {

    module 'python/x86_64/3.6.5'

    input:
    tuple val(version), file(maf) from originalMaf_ch

    output:
    tuple val(version), file("processedMaf.txt") into processedMaf_ch

    """
    #!/usr/bin/env python3

    loc_in = "$maf"

    names = []
    mafFormat1 = ["name1", "start1", "alnSize1", "strand1", "seqSize1", "alignment1"]
    mafFormat2 = ["name2", "start2", "alnSize2", "strand2", "seqSize2", "alignment2"]

    with open(loc_in) as inputFile:
        with open("processedMaf.txt", 'a') as outFile:
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
                        names.append("perc_id")
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
                    for i in range(len(align)):
                        record.append(align[i])
                    mismatch = 0
                    al1 = record[names.index("alignment1")]
                    al2 = record[names.index("alignment2")]
                    for i in range(len(al1)):
                        if al1[i] != al2[i]:
                            mismatch += 1
                    record.append(str(100*((len(al1)-mismatch)/float(len(al1)))))
                    outFile.write("\\t".join(record) + "\\n")
                    record = []
                    firstSequence = True
    """
}

process filterAndCastIntoBed {

    module 'R/x86_64/3.5.1'

    input:
    tuple val(version), file(processedMaf) from processedMaf_ch

    output:
    tuple val(version), file("filtered1.bed"), file("filtered2.bed") into originalBed_ch

    """
    #!/usr/bin/env Rscript

    mafTable <- read.table("$processedMaf", header = T, sep = "\\t")
    mafTable_filter <- mafTable[mafTable\$perc_id >= 99,]
    BED1 <- data.frame(name = mafTable_filter\$name1, 
                    start = mafTable_filter\$start1, 
                    end = as.integer(mafTable_filter\$start1+mafTable_filter\$alnSize1),
                    strand = mafTable_filter\$strand1)
    BED2 <- data.frame(name = mafTable_filter\$name2,
                    start = as.integer(ifelse(mafTable_filter\$strand2=="+", mafTable_filter\$start2, mafTable_filter\$seqSize2-mafTable_filter\$start2-mafTable_filter\$alnSize2)), 
                    end = as.integer(ifelse(mafTable_filter\$strand2=="+", mafTable_filter\$start2+mafTable_filter\$alnSize2, mafTable_filter\$seqSize2-mafTable_filter\$start2)), 
                    strand = mafTable_filter\$strand2)
    write.table(BED1, file = "filtered1.bed", sep = "\\t", quote = F, row.names = F, col.names = F)
    write.table(BED2, file = "filtered2.bed", sep = "\\t", quote = F, row.names = F, col.names = F)
    """
}

process mergeBedFiles {

    module 'bedtools/x86_64/2.2.28'

    publishDir "$params.outdir", mode: 'copy'

    input:
    tuple val(version), file(bed1), file(bed2) from originalBed_ch

    output:
    tuple val(version), file("*_merged1.bed"), file("*_merged2.bed") into mergedBed_ch

    """
    bedtools sort -i $bed1 | mergeBed -d 10 > ${version}_merged1.bed
    bedtools sort -i $bed2 | mergeBed -d 10 > ${version}_merged2.bed
    """
}

process calcFinalStats {

    module 'python/x86_64/3.6.5'

    publishDir "$params.outdir", mode: 'copy'

    input:
    tuple val(version), file(bed1), file(bed2), file(maskedGenome1), file(maskedGenome2) from mergedBed_ch.join(maskedGenomes2_ch)

    output:
    file "*_final_stats.txt" into final_ch

    """
    #!/usr/bin/env python3

    total_mapped1 = 0
    with open("$bed1") as inputFile:
        for line in inputFile:
            startEnd = line.rstrip().split("\\t")[1:3]
            total_mapped1 += (int(startEnd[1]) - int(startEnd[0]))
    
    total_mapped2 = 0
    with open("$bed2") as inputFile:
        for line in inputFile:
            startEnd = line.rstrip().split("\\t")[1:3]
            total_mapped2 += (int(startEnd[1]) - int(startEnd[0]))

    total_genome1 = 0
    masked1 = 0
    with open("$maskedGenome1") as inputFile:
        for line in inputFile:
            if line[0] != '>':
                total_genome1 += len(line.rstrip())
                masked1 += line.rstrip().count("N")
    unmasked1 = total_genome1 - masked1

    total_genome2 = 0
    masked2 = 0
    with open("$maskedGenome2") as inputFile:
        for line in inputFile:
            if line[0] != '>':
                total_genome2 += len(line.rstrip())
                masked2 += line.rstrip().count("N")
    unmasked2 = total_genome2 - masked2

    with open("${version}_final_stats.txt", 'w') as outFile:
        outFile.write("LAST mapped the other genome onto maize ${version}\\n$version is " + str(total_genome1) + " bp long,\\nthe other is " + str(total_genome2) + " bp long.\\n$version has " + str(unmasked1) + " bp unmasked,\\nthe other has " + str(unmasked2) + " bp unmasked.\\nThus, $version has " + str(100*(unmasked1/total_genome1)) + "% of its genome unmasked,\\nthe other genome is " + str(100*(unmasked2/total_genome2)) + "% unmasked.\\n$version has " + str(100*(total_mapped1/unmasked1)) + "% of its genome covered by the other genome mapped on it.\\n" + str(100*(total_mapped2/unmasked2)) + "% of the other genome is mapped on the $version genome.\\n")
    """
}
