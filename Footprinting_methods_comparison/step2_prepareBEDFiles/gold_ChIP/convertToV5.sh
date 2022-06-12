
BED_FOLDER=$1
CHAIN_FILE=$2
SEARCH_SPACE=$3

# Lift over to V5
module load liftOver
module load bedtools/x86_64/2.2.28
module load bedops/x86_64/2.4.37
mkdir unlifted
mkdir searchSpaceIntersected
for FILE in $BED_FOLDER/*.bed; do
    BASE=$(basename $FILE _0.01.bed)
    # Sort merge first
    bedtools sort -i $FILE | mergeBed | sort-bed - > sortedFile.tmp
    liftOver sortedFile.tmp $CHAIN_FILE ChIP_${BASE}_V5.bed ChIP_${BASE}_unlifted.bed
    # Clean up
    rm sortedFile.tmp
    mv ChIP_${BASE}_unlifted.bed unlifted
    # Intersect with search space
    bedtools sort -i ChIP_${BASE}_V5.bed | mergeBed | sort-bed - > sortedFile.tmp
    bedtools intersect -a sortedFile.tmp -b $SEARCH_SPACE | sort-bed - | mergeBed | sort-bed - > searchSpaceIntersected/ChIP_${BASE}.bed
done
