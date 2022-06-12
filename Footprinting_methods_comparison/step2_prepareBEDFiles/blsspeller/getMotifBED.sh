
# USAGE: getMotifBED.sh <path to folder with all bed files> <name of output file>

INPUT_FOLDER=$1
OUTPUT_NAME=$2
NONCODING_SEARCH_SPACE=/ngsprojects/blsmaize/results/Footprinting_methods_comparison/step1_getGenomeSubspaces/scaffoldsRemoved/2kb_noncoding.bed

module load bedops
module load bedtools
module load python/x86_64/3.6.5

# Filter for BLS threshold 0.98
python3 ./filterOnBLS.py $INPUT_FOLDER ${OUTPUT_NAME}_bls98_unshifted.bed 98

# Shift coordinates
python3 ./shiftBEDCoordinates.py ${OUTPUT_NAME}_bls98_unshifted.bed ${OUTPUT_NAME}_bls98.bed

# Clean BED by selecting first 3 columns & sort merging it
cut -f -3 ${OUTPUT_NAME}_bls98.bed > ${OUTPUT_NAME}_tmp.bed
sort-bed ${OUTPUT_NAME}_tmp.bed | mergeBed | sort-bed - > ${OUTPUT_NAME}.bed
rm ${OUTPUT_NAME}_tmp.bed
rm ${OUTPUT_NAME}_bls98.bed
rm ${OUTPUT_NAME}_bls98_unshifted.bed

# Intersect with the search space & remove the scaffolds
bedtools intersect -a $NONCODING_SEARCH_SPACE -b ${OUTPUT_NAME}.bed | sort-bed - | mergeBed | sort-bed - > fullSearch_${OUTPUT_NAME}_tmp.bed
grep '^[^s]' fullSearch_${OUTPUT_NAME}_tmp.bed > fullSearch_${OUTPUT_NAME}.bed
rm ${OUTPUT_NAME}.bed
rm fullSearch_${OUTPUT_NAME}_tmp.bed

# Add a label as extra column (for enrichment pipeline)
sed -i "s/$/\t${OUTPUT_NAME}/" fullSearch_${OUTPUT_NAME}.bed
