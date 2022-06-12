
# NEEDED FILES IN CURRENT FOLDER:
# v4_masked.fa (= v4_gen.fa.masked in assemblyCheck.nf output)
# v5_masked.fa (= v5_gen.fa.masked in assemblyCheck.nf output)
# v4.gff3
# v5.gff3
# v4-v5_overlap.bed (= V4_merged1.bed in assemblyCheck.nf output)
# v5-v4_overlap.bed (= V5_merged1.bed in assemblyCheck.nf output)
# getMaskedBed.py

# Get bed files of the full genome
module load samtools/x86_64/1.13
samtools faidx v4_masked.fa
samtools faidx v5_masked.fa
cat v4_masked.fa.fai | cut -f1,2 | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,0,$2}' > v4_full.bed
cat v5_masked.fa.fai | cut -f1,2 | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,0,$2}' > v5_full.bed

# Get bed files of the masked genome
module load python/x86_64/3.6.5
python3 getMaskedBed.py v4_masked.fa > v4_masked.bed
python3 getMaskedBed.py v5_masked.fa > v5_masked.bed

# Get bed files of the unmasked genome
module load bedtools/x86_64/2.30.0
bedtools sort -i v4_full.bed | mergeBed | sortBed > v4_full_s.bed
bedtools sort -i v4_masked.bed | mergeBed | sortBed > v4_masked_s.bed
bedtools sort -i v5_full.bed | mergeBed | sortBed > v5_full_s.bed
bedtools sort -i v5_masked.bed | mergeBed | sortBed > v5_masked_s.bed
bedtools subtract -nonamecheck -a v4_full_s.bed -b v4_masked_s.bed | sortBed | mergeBed | sortBed > v4_unmasked.bed
bedtools subtract -nonamecheck -a v5_full_s.bed -b v5_masked_s.bed | sortBed | mergeBed | sortBed > v5_unmasked.bed

# Get bed files of the unique fraction of the unmasked genome
bedtools sort -i v4-v5_overlap.bed | mergeBed | sortBed > v4-v5_overlap_s.bed
bedtools sort -i v5-v4_overlap.bed | mergeBed | sortBed > v5-v4_overlap_s.bed
bedtools subtract -nonamecheck -a v4_unmasked.bed -b v4-v5_overlap_s.bed | sortBed | mergeBed | sortBed > v4_unique.bed
bedtools subtract -nonamecheck -a v5_unmasked.bed -b v5-v4_overlap_s.bed | sortBed | mergeBed | sortBed > v5_unique.bed

# Get bed file of the annotated gene space
awk 'BEGIN {FS="\t"}; $3=="gene"' v4.gff3 > v4_genes.gff3
awk 'BEGIN {FS="\t"}; $3=="gene"' v5.gff3 > v5_genes.gff3
module load bedops/x86_64/2.4.37
gff2bed < v4_genes.gff3 | sortBed | mergeBed | sortBed > v4_genes.bed
gff2bed < v5_genes.gff3 | sortBed | mergeBed | sortBed > v5_genes.bed

# Get bed files of gene space in the interesting parts
bedtools intersect -nonamecheck -a v4_genes.bed -b v4_full_s.bed | sortBed | mergeBed > v4_genes_in_full.bed
bedtools intersect -nonamecheck -a v4_genes.bed -b v4_unmasked.bed | sortBed | mergeBed > v4_genes_in_unmasked.bed
bedtools intersect -nonamecheck -a v4_genes.bed -b v4-v5_overlap_s.bed | sortBed | mergeBed > v4_genes_in_overlap.bed
bedtools intersect -nonamecheck -a v4_genes.bed -b v4_unique.bed | sortBed | mergeBed > v4_genes_in_unique.bed
bedtools intersect -nonamecheck -a v5_genes.bed -b v5_full_s.bed | sortBed | mergeBed > v5_genes_in_full.bed
bedtools intersect -nonamecheck -a v5_genes.bed -b v5_unmasked.bed | sortBed | mergeBed > v5_genes_in_unmasked.bed
bedtools intersect -nonamecheck -a v5_genes.bed -b v5-v4_overlap_s.bed | sortBed | mergeBed > v5_genes_in_overlap.bed
bedtools intersect -nonamecheck -a v5_genes.bed -b v5_unique.bed | sortBed | mergeBed > v5_genes_in_unique.bed

# Get all the coverage numbers for the gene space
V4_FULL_GENE=$(cat v4_genes_in_full.bed | sortBed | mergeBed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
V4_UNMASKED_GENE=$(cat v4_genes_in_unmasked.bed | sortBed | mergeBed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
V4_OVERLAP_GENE=$(cat v4_genes_in_overlap.bed | sortBed | mergeBed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
V4_UNIQUE_GENE=$(cat v4_genes_in_unique.bed | sortBed | mergeBed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
V5_FULL_GENE=$(cat v5_genes_in_full.bed | sortBed | mergeBed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
V5_UNMASKED_GENE=$(cat v5_genes_in_unmasked.bed | sortBed | mergeBed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
V5_OVERLAP_GENE=$(cat v5_genes_in_overlap.bed | sortBed | mergeBed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
V5_UNIQUE_GENE=$(cat v5_genes_in_unique.bed | sortBed | mergeBed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
# Get all the coverage numbers for the entire subspace
V4_FULL=$(cat v4_full_s.bed | sortBed | mergeBed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
V4_UNMASKED=$(cat v4_unmasked.bed | sortBed | mergeBed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
V4_OVERLAP=$(cat v4-v5_overlap_s.bed | sortBed | mergeBed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
V4_UNIQUE=$(cat v4_unique.bed | sortBed | mergeBed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
V5_FULL=$(cat v5_full_s.bed | sortBed | mergeBed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
V5_UNMASKED=$(cat v5_unmasked.bed | sortBed | mergeBed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
V5_OVERLAP=$(cat v5-v4_overlap_s.bed | sortBed | mergeBed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
V5_UNIQUE=$(cat v5_unique.bed | sortBed | mergeBed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
# Calculate the gene coverage percentage
V4_FULL_PERC=$(echo "100*$V4_FULL_GENE/$V4_FULL" | bc -l)
V4_UNMASKED_PERC=$(echo "100*$V4_UNMASKED_GENE/$V4_UNMASKED" | bc -l)
V4_OVERLAP_PERC=$(echo "100*$V4_OVERLAP_GENE/$V4_OVERLAP" | bc -l)
V4_UNIQUE_PERC=$(echo "100*$V4_UNIQUE_GENE/$V4_UNIQUE" | bc -l)
V5_FULL_PERC=$(echo "100*$V5_FULL_GENE/$V5_FULL" | bc -l)
V5_UNMASKED_PERC=$(echo "100*$V5_UNMASKED_GENE/$V5_UNMASKED" | bc -l)
V5_OVERLAP_PERC=$(echo "100*$V5_OVERLAP_GENE/$V5_OVERLAP" | bc -l)
V5_UNIQUE_PERC=$(echo "100*$V5_UNIQUE_GENE/$V5_UNIQUE" | bc -l)

# Write everything in a results file
echo "For V4:" > results.txt
echo "In full genome: $V4_FULL bp out of $V4_FULL_GENE are gene --> this is $V4_FULL_PERC%" >> results.txt
echo "In unmasked genome: $V4_UNMASKED bp out of $V4_UNMASKED_GENE are gene --> this is $V4_UNMASKED_PERC%" >> results.txt
echo "In overlapping genome: $V4_OVERLAP bp out of $V4_OVERLAP_GENE are gene --> this is $V4_OVERLAP_PERC%" >> results.txt
echo "In unique genome: $V4_UNIQUE bp out of $V4_UNIQUE_GENE are gene --> this is $V4_UNIQUE_PERC%" >> results.txt
echo "For V5:" >> results.txt
echo "In full genome: $V5_FULL bp out of $V5_FULL_GENE are gene --> this is $V5_FULL_PERC%" >> results.txt
echo "In unmasked genome: $V5_UNMASKED bp out of $V5_UNMASKED_GENE are gene --> this is $V5_UNMASKED_PERC%" >> results.txt
echo "In overlapping genome: $V5_OVERLAP bp out of $V5_OVERLAP_GENE are gene --> this is $V5_OVERLAP_PERC%" >> results.txt
echo "In unique genome: $V5_UNIQUE bp out of $V5_UNIQUE_GENE are gene --> this is $V5_UNIQUE_PERC%" >> results.txt
