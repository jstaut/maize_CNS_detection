
originalMotifs_ch = Channel
						.fromPath("$params.input/*.bed", type: 'file' )

process filterBEDs {

	module 'python/x86_64/3.5.1'

	input:
	file bedFiles from originalMotifs_ch.collect()

	output:
	file "filtered.bed" into filteredBED_ch

	"""
	#!/usr/bin/env python3
	
	import shutil

	outLoc = "filtered.bed"
	targetGenesLoc = "targetGenes.txt"

	motifFiles = "$bedFiles".split(" ")
	threshold = $params.BLS
	
	with open(targetGenesLoc, 'w') as targetGenesOut:
		with open(outLoc, 'w') as outFile:
			for motifFile in motifFiles:
				targetGenes = set()
				with open(motifFile) as bedFile:
					for line in bedFile:
						fields = line.split('\\t')
						if 100*float(fields[6]) >= threshold:
							targetGenes.add(fields[4])
							outFile.write(line)
				targetGenesOut.write(str(len(targetGenes))+"\\n")
	
	shutil.copy(targetGenesLoc, "$params.outdir")
	"""
}

process calculateCoverage {

	module 'bedtools/x86_64/2.2.28'
	module 'bedops/x86_64/2.4.37'

	input:
	file filteredBED from filteredBED_ch

	output:
	tuple env(COVERAGE_ORIGINAL), env(COVERAGE_FULL), env(COVERAGE_INTRONS), env(COVERAGE_5_UTR), env(COVERAGE_3_UTR), env(COVERAGE_NONGENE) into gather_ch

	"""
	sort-bed $filteredBED | mergeBed | sort-bed - > filtered_sorted.bed
	bedtools intersect -a filtered_sorted.bed -b $params.fullBED | sort-bed - | mergeBed | sort-bed - > fullSearchSpace.bed
	bedtools intersect -a filtered_sorted.bed -b $params.intronBED | sort-bed - | mergeBed | sort-bed - > introns.bed
	bedtools intersect -a filtered_sorted.bed -b $params.UTR5BED | sort-bed - | mergeBed | sort-bed - > 5_UTR.bed
	bedtools intersect -a filtered_sorted.bed -b $params.UTR3BED | sort-bed - | mergeBed | sort-bed - > 3_UTR.bed
	bedtools intersect -a filtered_sorted.bed -b $params.nonGeneBED | sort-bed - | mergeBed | sort-bed - > nonGene.bed

	COVERAGE_ORIGINAL=\$(cat filtered_sorted.bed | awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}')
	COVERAGE_FULL=\$(cat fullSearchSpace.bed | awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}')
	COVERAGE_INTRONS=\$(cat introns.bed | awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}')
	COVERAGE_5_UTR=\$(cat 5_UTR.bed | awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}')
	COVERAGE_3_UTR=\$(cat 3_UTR.bed | awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}')
	COVERAGE_NONGENE=\$(cat nonGene.bed | awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{print SUM}')

	"""
}

process gatherResults {

	publishDir "$params.outdir", mode: 'copy'

	input:
	tuple val(coverage_original), val(coverage_full), val(coverage_introns), val(coverage_5_utr), val(coverage_3_utr), val(coverage_nongene) from gather_ch

	output:
	file "coverageResults.txt" into results_ch

	"""
	
	echo "Original	$coverage_original	\$(echo '100 * $coverage_original / $params.fullCov' | bc -l)%" > coverageResults.txt
	echo "FullSearch	$coverage_full	\$(echo '100 * $coverage_full / $params.fullCov' | bc -l)%" >> coverageResults.txt
	echo "Introns	$coverage_introns	\$(echo '100 * $coverage_introns / $params.intronCov' | bc -l)%" >> coverageResults.txt
	echo "5' UTR	$coverage_5_utr	\$(echo '100 * $coverage_5_utr / $params.utr5Cov' | bc -l)%" >> coverageResults.txt
	echo "3' UTR	$coverage_3_utr	\$(echo '100 * $coverage_3_utr / $params.utr3Cov' | bc -l)%" >> coverageResults.txt
	echo "Nongene	$coverage_nongene	\$(echo '100 * $coverage_nongene / $params.nongeneCov' | bc -l)%" >> coverageResults.txt

	"""
}

process getTargetGeneStats {

	module 'python/x86_64/3.5.1'

	publishDir "$params.outdir", mode: 'copy'

	input:
	file $results from results_ch

	output:
	file "targetGeneStats.txt" into final_ch

	"""
	#!/usr/bin/env python3
	
	import os
	os.environ["OMP_NUM_THREADS"] = "1"
	from statistics import median, mean
	import numpy

	with open("targetGeneStats.txt", 'w') as outFile:
		nbOfTargetGenes = list()
		with open("$params.outdir/targetGenes.txt") as f:
			for line in f:
				nbOfTargetGenes.append(int(line.rstrip()))
		outFile.write("Median\\t"+str(median(nbOfTargetGenes))+"\\n")
		outFile.write("Mean\\t"+str(mean(nbOfTargetGenes))+"\\n")
		outFile.write("Perc75\\t"+str(numpy.percentile(nbOfTargetGenes, 75))+"\\n")
		outFile.write("Perc90\\t"+str(numpy.percentile(nbOfTargetGenes, 90))+"\\n")
	"""
}
