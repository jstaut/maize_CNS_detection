
log.info """\
         PROMOTER EXTRACTION NEXTFLOW PIPELINE
         =====================================
         input dir: ${params.input}
         """
         .stripIndent()

plaza_id_ch = Channel
                .fromPath("$params.input/PLAZA_IDs.txt", type: 'file' )
                .splitText()
                .map{it -> it.trim()}

maizegdb_id_ch = Channel
                .fromPath("$params.input/maizeGDB_IDs.txt", type: 'file' )
                .splitText()
                .map{it -> it.trim()}

process getPlazaFiles {

    input:
    val ID from plaza_id_ch

    output:
    tuple val(ID), file("*_genome.fasta"), file("*.gff3") into plaza_files_ch

    """
    wget -O ${ID}_genome.fasta.gz ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_05/Genomes/${ID}.fasta.gz
    # wget -O ${ID}_proteome.fasta.gz ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_05/Fasta/proteome.selected_transcript.${ID}.fasta.gz
    wget -O ${ID}.gff3.gz ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_05/GFF/${ID}/annotation.selected_transcript.all_features.${ID}.gff3.gz
    gunzip ${ID}_genome.fasta.gz
    # gunzip ${ID}_proteome.fasta.gz
    gunzip ${ID}.gff3.gz
    """
}

process getMaizeGDBFiles {

    input:
    val ID from maizegdb_id_ch

    output:
    tuple env(new_ID), file("*_genome.fasta"), file("*_proteome.fasta"), file("*.gff3") into maizeGDB_files_1_ch

    """
    prim_ID=\$(echo $ID | cut -f1 -d_)
    sec_ID=\$(echo $ID | cut -f2 -d_)
    new_ID=Zm-\${prim_ID}-NAM
    wget -O \${new_ID}_genome.fasta.gz https://download.maizegdb.org/Zm-\${prim_ID}-REFERENCE-NAM-1.0/Zm-\${prim_ID}-REFERENCE-NAM-1.0.fa.gz
    wget -O \${new_ID}_proteome.fasta.gz https://download.maizegdb.org/Zm-\${prim_ID}-REFERENCE-NAM-1.0/Zm-\${prim_ID}-REFERENCE-NAM-1.0_\${sec_ID}.1.protein.fa.gz
    wget -O \${new_ID}.gff3.gz https://download.maizegdb.org/Zm-\${prim_ID}-REFERENCE-NAM-1.0/Zm-\${prim_ID}-REFERENCE-NAM-1.0_\${sec_ID}.1.gff3.gz
    gunzip \${new_ID}_genome.fasta.gz
    gunzip \${new_ID}_proteome.fasta.gz
    gunzip \${new_ID}.gff3.gz
    """
}

process filterGFF {

    module 'python/x86_64/3.6.5'

    input:
    tuple val(ID), file(genome), file(proteome), file(gff) from maizeGDB_files_1_ch

    output:
    tuple val(ID), file(genome), file("*_filtered.gff3") into maizeGDB_files_2_ch

    """
    #!/usr/bin/env python3

    from Bio import SeqIO

    loc_gff = "$gff"
    loc_prot = "$proteome"
    loc_out = "${ID}_filtered.gff3"

    protDict = dict()
    count = 0

    for prot in SeqIO.parse(open(loc_prot),'fasta'):
        generalID = prot.id.split("_")[0]
        isoIDNb = prot.id.split("_")[1][-3:]
        if generalID not in protDict.keys():
            protDict[generalID] = [str(prot.seq), isoIDNb]
        else:
            if len(protDict[generalID][0]) < len(prot.seq):
                protDict[generalID] = [str(prot.seq), isoIDNb]
            elif len(protDict[generalID][0]) == len(prot.seq) and protDict[generalID][0] != str(prot.seq):
                count += 1
                if protDict[generalID][1] > isoIDNb:
                    protDict[generalID] = [str(prot.seq), isoIDNb]

    selected_prot = set()
    for key, value in protDict.items():
        selected_prot.add("_".join([key,"P"+str(value[1])]))

    selected_transcript = set()
    with open(loc_gff) as file:
        for line in file:
            if not line.startswith('#'):
                rec = line.rstrip().split('\\t')
                if rec[2] == "CDS":
                    prot_iso = rec[8].split(';')[0].split('=')[1]
                    if prot_iso in selected_prot:
                        transcr_iso = rec[8].split(';')[1].split('=')[1]
                        selected_transcript.add(transcr_iso)

    with open(loc_gff) as file:
        with open(loc_out, 'a') as out:
            for line in file:
                if line.startswith('#'):
                    out.write(line)
                else:
                    rec = line.rstrip().split('\\t')
                    if rec[2] == "CDS":
                        prot_iso = rec[8].split(';')[0].split('=')[1]
                        if prot_iso in selected_prot:
                            out.write(line)
                    elif rec[2] == "mRNA" or rec[2] == "exon" or rec[2] == "five_prime_UTR" or rec[2] == "three_prime_UTR":
                        transcr_iso = rec[8].split(';')[0].split('=')[1]
                        if transcr_iso in selected_transcript:
                            out.write(line)
                    else:
                        out.write(line)
    """
}

required_files_ch = plaza_files_ch.mix(maizeGDB_files_2_ch)

(required_files1_ch, required_files2_ch) = required_files_ch.into(2)

process GetNbrhoodGFF {
    
    module 'python/x86_64/3.6.5'
    module 'bedtools'
    module 'samtools'

    input:
    tuple val(ID), file(genome), file(gff) from required_files1_ch

    output:
    tuple ID, file("${ID}_${params.up_window}up_${params.down_window}down.gff3") into nbrhood_gff_ch

    """

    python3 $params.pyUtilDir/main.py $gff $genome -o${ID} -up $params.up_window -down $params.down_window -iso

    """
}

nbrhood_gff_ch.join(required_files2_ch).set{ files_ch }

process GetNbrhoodFasta {
    
    publishDir "$params.outdir", mode: 'move'

    module 'python/x86_64/3.6.5'
    module 'bedtools'
    module 'samtools'

    input:
    tuple val(ID), file(nbrhood_gff), file(genome), file(gff) from files_ch

    output:
    tuple file("unmasked_promoter_${ID}.fasta"), file("masked_promoter_${ID}.fasta") into fasta_prom_ch

    """
    awk 'BEGIN {FS="\\t"} \$3=="exon"' $gff | sortBed > exonsWithUTR_${ID}.gff3

    awk 'BEGIN {FS="\\t"} \$3=="five_prime_UTR" || \$3=="three_prime_UTR"' $gff | sortBed > UTR_${ID}.gff3

    subtractBed -a exonsWithUTR_${ID}.gff3 -b UTR_${ID}.gff3 -sorted > exons_${ID}.gff3

    python3 $params.pyUtilDir/getExonCoords.py $gff exons_${ID}.gff3

    bedtools maskfasta -fi $genome -bed exons_${ID}.gff3 -fo masked_${ID}.fasta

    samtools faidx masked_${ID}.fasta
    samtools faidx $genome

    bedtools getfasta -fi masked_${ID}.fasta -bed $nbrhood_gff -name+ > masked_promoter_${ID}.fasta
    bedtools getfasta -fi $genome -bed $nbrhood_gff -name+ > unmasked_promoter_${ID}.fasta

    """
}

workflow.onComplete {
	log.info ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
