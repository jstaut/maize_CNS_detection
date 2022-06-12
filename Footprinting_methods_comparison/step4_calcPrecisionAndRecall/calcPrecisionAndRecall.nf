
cnsBEDs_ch = Channel.fromPath("$params.CNS_BEDs/*.bed", type: 'file' )
referenceBEDs_ch = Channel.fromPath("$params.reference_BEDs/*.bed", type: 'file' )

process calculateCoverage {

	module 'bedtools/x86_64/2.2.28'
	module 'bedops/x86_64/2.4.37'

	input:
	tuple file(cnsBED), file(referenceBED) from cnsBEDs_ch.combine(referenceBEDs_ch)

	output:
	file "*.txt" into coverage_ch

	"""
	CNS_BED=${cnsBED.baseName}
	REF_BED=${referenceBED.baseName}

	# Intersect with the search space
	bedtools intersect -a $cnsBED -b $params.searchSpace | sort-bed - | mergeBed | sort-bed - > CNSs.bed
	bedtools intersect -a $referenceBED -b $params.searchSpace | sort-bed - | mergeBed | sort-bed - > reference.bed

	# Calculate overlap between the CNSs and reference (ACR/UMR/ChIP)
	bedtools intersect -a CNSs.bed -b reference.bed | sort-bed - | mergeBed | sort-bed - > intersection.bed

	# Calculate the bp coverage
	COV_SEARCH=\$(cat $params.searchSpace | awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}')
	COV_CNS=\$(cat CNSs.bed | awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}')
	COV_REF=\$(cat reference.bed | awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}')
	COV_OVERLAP=\$(cat intersection.bed | awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}')

	PRECISION_PERC=\$(echo "100*\$COV_OVERLAP/\$COV_CNS" | bc -l)
	RECALL_PERC=\$(echo "100*\$COV_OVERLAP/\$COV_REF" | bc -l)
	PRECISION_FRACT=\$(echo "\$COV_OVERLAP/\$COV_CNS" | bc -l)
	RECALL_FRACT=\$(echo "\$COV_OVERLAP/\$COV_REF" | bc -l)
	F1_PERC=\$(echo "100*(2*\$PRECISION_FRACT*\$RECALL_FRACT)/(\$PRECISION_FRACT+\$RECALL_FRACT)" | bc -l)
	ENRICHMENT=\$(echo "(\$COV_OVERLAP/\$COV_REF)/(\$COV_CNS/\$COV_SEARCH)" | bc -l)

	# To keep storage low
	rm CNSs.bed
	rm reference.bed
	rm intersection.bed

	# Write results to file
	echo "\$REF_BED	\$CNS_BED	\$PRECISION_PERC	\$RECALL_PERC	\$F1_PERC	\$ENRICHMENT	\$COV_SEARCH	\$COV_CNS	\$COV_REF	\$COV_OVERLAP" > \${REF_BED}_\${CNS_BED}.txt
	"""
}

process combineStats {

	module 'python/x86_64/3.5.1'

	publishDir "$params.outdir", mode: 'move'

	input:
	file statFiles from coverage_ch.collect()

	output:
	file "allStats.txt" into final_ch

	"""
	echo "RefData	CNSData	Precision	Recall	F1	Enrichment-fold	SearchSpace	CNSCoverage	ReferenceCoverage	OverlapCoverage" > allStats.txt
	for statFile in $statFiles
	do
		cat \$statFile >> allStats.txt
	done

	"""
}
