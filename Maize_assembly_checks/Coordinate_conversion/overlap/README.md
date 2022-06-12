Resulting BED files of method 1 and 2 were sorted and merged with bedtools. Then the overlapping regions were calculated with bedtools intersect. 

The following command was used to count the total length of all BED file regions:

cat overlap.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
