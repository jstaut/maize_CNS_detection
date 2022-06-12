
motifBEDs_ch = Channel.fromPath("$params.input_motifs/*.bed", type: 'file' )
subSpaceBEDs_ch = Channel.fromPath("$params.input_subSpace/*.bed", type: 'file' )

process calculateCoverage {

	module 'bedtools/x86_64/2.2.28'
	module 'bedops/x86_64/2.4.37'

	input:
	tuple file(motifBED), file(subSpaceBED) from motifBEDs_ch.combine(subSpaceBEDs_ch)

	output:
	file "*.txt" into coverage_ch
	file "*.bed" into subMotifs_ch

	"""
	MOTIF_BED=${motifBED.baseName}
	SPACE_BED=${subSpaceBED.baseName}
	bedtools intersect -a $motifBED -b $subSpaceBED > \${MOTIF_BED}.\${SPACE_BED}.bed
	SIZE_MOTIFS=\$(cat \${MOTIF_BED}.\${SPACE_BED}.bed | awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}')
	SIZE_INTERVAL=\$(cat $subSpaceBED | awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}')
	PERCENTAGE=\$(echo "100*\$SIZE_MOTIFS/\$SIZE_INTERVAL" | bc -l)
	CNS_NUMBER=\$(cat \${MOTIF_BED}.\${SPACE_BED}.bed | wc -l)
	echo "\$SPACE_BED	\$SIZE_MOTIFS	\$SIZE_INTERVAL	\$PERCENTAGE	\$CNS_NUMBER" > \$RANDOM.txt
	echo "\$SPACE_BED	\$SIZE_MOTIFS	\$SIZE_INTERVAL	\$PERCENTAGE	\$CNS_NUMBER" >> $params.outdir/\${MOTIF_BED}_coverage.txt
	"""
}

process calculateSizeInfo {

	module 'python/x86_64/3.5.1'

	input:
	file subMotifs from subMotifs_ch

	"""
	#!/usr/bin/env python3
	
	import os
	os.environ["OMP_NUM_THREADS"] = "1"
	import pandas as pd
	from statistics import mean
	import numpy as np

	loc = "$subMotifs"
	baseName = loc.split("/")[-1].split(".")[0]
	secondName = loc.split("/")[-1].split(".")[1]

	bedFile = pd.read_csv(loc, sep='\\t', low_memory=False, header=None)
	cnsSizes = bedFile.iloc[:,2]-bedFile.iloc[:,1]

	with open("$params.outdir/"+baseName+"_"+secondName+"_CNSSizes.txt", 'w') as outFile:
		for size in cnsSizes.sort_values():
			outFile.write(str(size)+"\\n")

	mean = mean(cnsSizes)
	minSize = min(cnsSizes)
	perc25 = np.percentile(cnsSizes, 25)
	perc50 = np.percentile(cnsSizes, 50)
	perc75 = np.percentile(cnsSizes, 75)
	perc90 = np.percentile(cnsSizes, 90)
	maxSize = max(cnsSizes)
	stats = "\\t".join([secondName, str(mean), str(minSize), str(perc25), str(perc50), str(perc75), str(perc90), str(maxSize)])

	with open("$params.outdir/"+baseName+"_lengthStats.txt", 'a') as outFile:
		outFile.write(stats+"\\n")
	"""
}
