
# USAGE: ./getFeatureBEDs.sh <full species gff> <fasta index file> <search space BED> <searc space naming label>

FULL_GFF=$1
FAIDX=$2
SEARCH_SPACE=$3
SEARCH_SPACE_LABEL=$4

module load bedops
module load bedtools
grep -P '\tCDS\t' $FULL_GFF | gff2bed | sort-bed - | mergeBed | sort-bed - > CDS.bed
grep -P '\tfive_prime_UTR\t' $FULL_GFF | gff2bed | sort-bed - | mergeBed | sort-bed - > 5_UTR.bed
grep -P '\tthree_prime_UTR\t' $FULL_GFF | gff2bed | sort-bed - | mergeBed | sort-bed - > 3_UTR.bed
grep -P '\tgene\t' $FULL_GFF | gff2bed | sort-bed - | mergeBed | sort-bed - > gene.bed
cat $FAIDX | cut -f1,2 | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,0,$2}' > full.bed

# Create intron, noncoding and nongene files out of the previous basic files
bedtools subtract -a gene.bed -b CDS.bed | sort-bed - | mergeBed | sort-bed - > tmp1.bed
bedtools subtract -a tmp1.bed -b 5_UTR.bed | sort-bed - | mergeBed | sort-bed - > tmp2.bed
bedtools subtract -a tmp2.bed -b 3_UTR.bed | sort-bed - | mergeBed | sort-bed - > introns.bed
rm tmp1.bed
rm tmp2.bed
bedtools subtract -a full.bed -b CDS.bed | sort-bed - | mergeBed | sort-bed - > noncoding.bed
bedtools subtract -a full.bed -b gene.bed | sort-bed - | mergeBed | sort-bed - > nongene.bed

# Sort & merge search space bed + subtract coding regions
bedtools sort -i $SEARCH_SPACE | mergeBed | sort-bed - > ${SEARCH_SPACE_LABEL}.bed
bedtools intersect -a ${SEARCH_SPACE_LABEL}.bed -b noncoding.bed | sort-bed - | mergeBed | sort-bed - > ${SEARCH_SPACE_LABEL}_noncoding.bed

# Intersect with search space
bedtools intersect -a ${SEARCH_SPACE_LABEL}_noncoding.bed -b CDS.bed | sort-bed - | mergeBed | sort-bed - > ${SEARCH_SPACE_LABEL}_CDS.bed
bedtools intersect -a ${SEARCH_SPACE_LABEL}_noncoding.bed -b 5_UTR.bed | sort-bed - | mergeBed | sort-bed - > ${SEARCH_SPACE_LABEL}_5_UTR.bed
bedtools intersect -a ${SEARCH_SPACE_LABEL}_noncoding.bed -b 3_UTR.bed | sort-bed - | mergeBed | sort-bed - > ${SEARCH_SPACE_LABEL}_3_UTR.bed
bedtools intersect -a ${SEARCH_SPACE_LABEL}_noncoding.bed -b gene.bed | sort-bed - | mergeBed | sort-bed - > ${SEARCH_SPACE_LABEL}_gene.bed
bedtools intersect -a ${SEARCH_SPACE_LABEL}_noncoding.bed -b full.bed | sort-bed - | mergeBed | sort-bed - > ${SEARCH_SPACE_LABEL}_full.bed
bedtools intersect -a ${SEARCH_SPACE_LABEL}_noncoding.bed -b introns.bed | sort-bed - | mergeBed | sort-bed - > ${SEARCH_SPACE_LABEL}_introns.bed
bedtools intersect -a ${SEARCH_SPACE_LABEL}_noncoding.bed -b nongene.bed | sort-bed - | mergeBed | sort-bed - > ${SEARCH_SPACE_LABEL}_nongene.bed

# Calculate genome coverage and write to QC file
# Genome
echo "Genome size:" > sizes_QC.txt
GENOME_SIZE=$(cat full.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
echo $GENOME_SIZE >> sizes_QC.txt
# Noncoding
echo "Noncoding size:" >> sizes_QC.txt
NONCODING_SIZE=$(cat noncoding.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
NONCODING_PERC=$(echo "100*$NONCODING_SIZE/$GENOME_SIZE" | bc -l)
echo "$NONCODING_SIZE ($NONCODING_PERC %)" >> sizes_QC.txt
# Gene
echo "Gene size:" >> sizes_QC.txt
GENE_SIZE=$(cat gene.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
GENE_PERC=$(echo "100*$GENE_SIZE/$GENOME_SIZE" | bc -l)
echo "$GENE_SIZE ($GENE_PERC %)" >> sizes_QC.txt
# Nongene
echo "Nongene size:" >> sizes_QC.txt
NONGENE_SIZE=$(cat nongene.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
NONGENE_PERC=$(echo "100*$NONGENE_SIZE/$GENOME_SIZE" | bc -l)
echo "$NONGENE_SIZE ($NONGENE_PERC %)" >> sizes_QC.txt
# 5' UTR
echo "5' UTR size:" >> sizes_QC.txt
UTR_5_SIZE=$(cat 5_UTR.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
UTR_5_PERC=$(echo "100*$UTR_5_SIZE/$GENOME_SIZE" | bc -l)
echo "$UTR_5_SIZE ($UTR_5_PERC %)" >> sizes_QC.txt
# 3' UTR
echo "3' UTR size:" >> sizes_QC.txt
UTR_3_SIZE=$(cat 3_UTR.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
UTR_3_PERC=$(echo "100*$UTR_3_SIZE/$GENOME_SIZE" | bc -l)
echo "$UTR_3_SIZE ($UTR_3_PERC %)" >> sizes_QC.txt
# Introns
echo "Introns size:" >> sizes_QC.txt
INTRON_SIZE=$(cat introns.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
INTRON_PERC=$(echo "100*$INTRON_SIZE/$GENOME_SIZE" | bc -l)
echo "$INTRON_SIZE ($INTRON_PERC %)" >> sizes_QC.txt
# CDS
echo "CDS size:" >> sizes_QC.txt
CDS_SIZE=$(cat CDS.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
CDS_PERC=$(echo "100*$CDS_SIZE/$GENOME_SIZE" | bc -l)
echo "$CDS_SIZE ($CDS_PERC %)" >> sizes_QC.txt
# Full search space
echo "Search space size (noncoding):" >> sizes_QC.txt
SEARCH_SPACE_SIZE=$(cat ${SEARCH_SPACE_LABEL}_noncoding.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
SEARCH_SPACE_PERC=$(echo "100*$SEARCH_SPACE_SIZE/$GENOME_SIZE" | bc -l)
echo "$SEARCH_SPACE_SIZE ($SEARCH_SPACE_PERC %)" >> sizes_QC.txt

# Calculate genome coverage relative to search space and write to QC file
# Full search space
echo "Search space size (noncoding):" > ${SEARCH_SPACE_LABEL}_sizes_QC.txt
echo $SEARCH_SPACE_SIZE >> ${SEARCH_SPACE_LABEL}_sizes_QC.txt
# Gene
echo "Gene size:" >> ${SEARCH_SPACE_LABEL}_sizes_QC.txt
GENE_SIZE=$(cat ${SEARCH_SPACE_LABEL}_gene.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
GENE_PERC=$(echo "100*$GENE_SIZE/$SEARCH_SPACE_SIZE" | bc -l)
echo "$GENE_SIZE ($GENE_PERC %)" >> ${SEARCH_SPACE_LABEL}_sizes_QC.txt
# Nongene
echo "Nongene size:" >> ${SEARCH_SPACE_LABEL}_sizes_QC.txt
NONGENE_SIZE=$(cat ${SEARCH_SPACE_LABEL}_nongene.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
NONGENE_PERC=$(echo "100*$NONGENE_SIZE/$SEARCH_SPACE_SIZE" | bc -l)
echo "$NONGENE_SIZE ($NONGENE_PERC %)" >> ${SEARCH_SPACE_LABEL}_sizes_QC.txt
# 5' UTR
echo "5' UTR size:" >> ${SEARCH_SPACE_LABEL}_sizes_QC.txt
UTR_5_SIZE=$(cat ${SEARCH_SPACE_LABEL}_5_UTR.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
UTR_5_PERC=$(echo "100*$UTR_5_SIZE/$SEARCH_SPACE_SIZE" | bc -l)
echo "$UTR_5_SIZE ($UTR_5_PERC %)" >> ${SEARCH_SPACE_LABEL}_sizes_QC.txt
# 3' UTR
echo "3' UTR size:" >> ${SEARCH_SPACE_LABEL}_sizes_QC.txt
UTR_3_SIZE=$(cat ${SEARCH_SPACE_LABEL}_3_UTR.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
UTR_3_PERC=$(echo "100*$UTR_3_SIZE/$SEARCH_SPACE_SIZE" | bc -l)
echo "$UTR_3_SIZE ($UTR_3_PERC %)" >> ${SEARCH_SPACE_LABEL}_sizes_QC.txt
# Introns
echo "Introns size:" >> ${SEARCH_SPACE_LABEL}_sizes_QC.txt
INTRON_SIZE=$(cat ${SEARCH_SPACE_LABEL}_introns.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
INTRON_PERC=$(echo "100*$INTRON_SIZE/$SEARCH_SPACE_SIZE" | bc -l)
echo "$INTRON_SIZE ($INTRON_PERC %)" >> ${SEARCH_SPACE_LABEL}_sizes_QC.txt
# CDS
echo "CDS size:" >> ${SEARCH_SPACE_LABEL}_sizes_QC.txt
CDS_SIZE=$(cat ${SEARCH_SPACE_LABEL}_CDS.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
CDS_PERC=$(echo "100*$CDS_SIZE/$SEARCH_SPACE_SIZE" | bc -l)
echo "$CDS_SIZE ($CDS_PERC %)" >> ${SEARCH_SPACE_LABEL}_sizes_QC.txt
