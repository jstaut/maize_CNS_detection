
htmls_ch = Channel
                .fromPath("$params.input/*.txt", type: 'file' )
                .splitText()
                .map{it -> it.trim()}

fasta_ch = Channel
                .fromPath("$params.input/*.fasta", type: 'file' )

process getHtmlProteome {

    input:
    val html from htmls_ch

    output:
    tuple env(SPECIES), file("*.fasta"), file("*_QC.txt") into proteome_ch

    """
    SPECIES=`echo $html | grep -oP '(?<=Zm-)\\w*(?=-)' | head -n 1 | awk '{print "Zm-" \$1 "-NAM"}'`
    wget -O proteome.fasta.gz $html
    gunzip -c proteome.fasta.gz > \$SPECIES.fasta
    echo 'Number of sequences in original file:' > \${SPECIES}_QC.txt
    grep '>' \$SPECIES.fasta | wc -l >> \${SPECIES}_QC.txt
    """
}

process getFileProteome {

    input:
    file fasta from fasta_ch

    output:
    tuple env(SPECIES), file(fasta), file("*_QC.txt") into alreadyPrimTrans_ch

    """
    SPECIES=`echo $fasta | grep -oP '(?<=proteome.selected_transcript.)[\\w-_]*(?=.fasta)'`
    echo 'Number of sequences in original file:' > \${SPECIES}_QC.txt
    grep '>' $fasta | wc -l >> \${SPECIES}_QC.txt
    """
}

process getPrimaryTranscripts {

    module 'python/x86_64/3.6.5'

    input:
    tuple val(species), file(proteome), file(QC) from proteome_ch

    output:
    tuple val(species), file("*_primTrans*"), file(QC) into newPrimTrans_ch

    """
    #!/usr/bin/env python3

    from Bio import SeqIO

    inputFile = "$proteome"
    outputFile = "$proteome"
    outputFile = outputFile.split(".")
    outputFile[0] += "_primTrans"
    outputFile = ".".join(outputFile)

    protDict = dict()
    count = 0

    for prot in SeqIO.parse(open(inputFile),'fasta'):
        generalID = prot.id.split("_")[0]
        try:
            isoIDNb = prot.id.split("_")[1][-3:]
        except IndexError:
            print(prot.id)
            print(loc_in)
        if generalID not in protDict.keys():
            protDict[generalID] = [str(prot.seq), isoIDNb]
        else:
            if len(protDict[generalID][0]) < len(prot.seq):
                protDict[generalID] = [str(prot.seq), isoIDNb]
            elif len(protDict[generalID][0]) == len(prot.seq) and protDict[generalID][0] != str(prot.seq):
                count += 1
                if protDict[generalID][1] > isoIDNb:
                    protDict[generalID] = [str(prot.seq), isoIDNb]

    with open(outputFile, 'a') as outFile:
        for key, value in protDict.items():
            outFile.write(">"+key+"\\n"+value[0]+"\\n")
    
    with open("$QC", 'a') as QCFile:
        QCFile.write("Extracting primary transcripts...\\n")
        QCFile.write("Detected "+str(count)+" isoforms with same length but different sequence\\n -> Selected the ones with smallest isoform ID number\\n")
        QCFile.write("Number of extracted primary transcripts:\\n"+str(len(protDict))+"\\n")
    """
}

primTrans_ch = alreadyPrimTrans_ch.mix(newPrimTrans_ch)

process cleanProteome {

    module 'python/x86_64/3.6.5'

    input:
    tuple val(species), file(proteome), file(QC) from primTrans_ch

    output:
    tuple val(species), file("*_cleaned*"), file(QC) into cleaned_ch

    """
    #!/usr/bin/env python3

    from Bio import SeqIO

    inputFile = "$proteome"
    outputFile = "$proteome"
    outputFile = outputFile.split(".")
    outputFile[0] += "_cleaned"
    outputFile = ".".join(outputFile)

    with open(outputFile, 'a') as outFile:
        for prot in SeqIO.parse(open(inputFile),'fasta'):
            if "|" in prot.description:
                outFile.write(">"+(prot.description.split("|")[1].lstrip())+"\\n")
            else:
                outFile.write(">"+prot.id+"\\n")
            outFile.write(str(prot.seq).replace("*","")+"\\n")
    """
}

process runOrthoFinder {

    module 'diamond'
    module 'mcl'
    module 'fastme'
    module 'OrthoFinder'

    input:
    file proteomes from cleaned_ch.collect()

    output:
    file 'results/*' into orthoOutput_ch

    """
    mkdir neededProteomes
    mkdir results
    for input in $proteomes
    do 
        if \$(echo \$input | grep -q '^input'); then
            SPECIES=`cat \$input`
        elif \$(echo \$input | grep -q '.fasta\$'); then
            cp \$input ./neededProteomes/\$SPECIES.fasta
        elif \$(echo \$input | grep -q '_QC.txt\$'); then
            cp \$input ./results
        fi
    done
    
    orthofinder -t $params.orthoThreads -f ./neededProteomes
    find . -type f -name "Orthogroups.GeneCount.tsv" -exec sh -c 'cp "\$0" results' {} \\;
    find . -type f -name "Orthogroups.tsv" -exec sh -c 'cp "\$0" results' {} \\;
    find . -type f -name "Orthogroups.txt" -exec sh -c 'cp "\$0" results' {} \\;
    find . -type f -name "SpeciesTree_rooted_node_labels.txt" -exec sh -c 'cp "\$0" results' {} \\;
    find . -type f -name "SpeciesTree_rooted.txt" -exec sh -c 'cp "\$0" results' {} \\;
    tar -zcf ./neededProteomes/OrthoFinder.tar.gz ./neededProteomes/OrthoFinder
    mv ./neededProteomes/OrthoFinder.tar.gz results
    """
}

process cleanOrthoGroups {

    module 'python/x86_64/3.6.5'

    publishDir "$params.outdir", mode: 'copy'

    input:
    file orthoOutput from orthoOutput_ch

    output:
    file 'results/**' into final_ch

    """
    #!/usr/bin/env python3

    import os
    from Bio import SeqIO
    import pandas as pd
    from shutil import copy
    
    # Create final folder structure to be published
    os.makedirs("results/QC_files")

    # Separate QC files from the rest of the files
    QC_files = dict()
    other_files = []
    for someFile in "$orthoOutput".split(" "):
        if someFile.endswith("_QC.txt"):
            QC_files[someFile[:-7]] = someFile
        else:
            other_files.append(someFile)
    
    # Only retain multi-species orthogroups -> "cleaned orthogroups"
    orthoGroups = pd.read_csv("Orthogroups.tsv", sep='\\t', header=0)
    orthoGroups = orthoGroups.iloc[(orthoGroups.iloc[:,1:].isnull().sum(axis=1)<(orthoGroups.shape[1]-2)).values,:]
    orthoGroups.to_csv("Orthogroups_cleaned.tsv", sep='\\t', index=False)
    other_files.append("Orthogroups_cleaned.tsv")

    # Put files in result directory
    for aFile in other_files:
        copy(aFile, "results")

    # Calculate QC stats and add to QC file
    for species in orthoGroups.iloc[:,1:].columns:
        genesInOrthoGroups = sum([gene.count(",")+1 for gene in orthoGroups[species].dropna()])
        with open(QC_files[species], 'a') as QC_file:
            QC_file.write("Number of primary transcripts in (multi-species) orthogroups:\\n"+str(genesInOrthoGroups))
        # Add QC files to result directory
        copy(QC_files[species], "results/QC_files")
    """
}
