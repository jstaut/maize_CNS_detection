
originalMotifs_ch = Channel
						.fromPath("$params.input", type: 'dir' )

blsRange_ch = Channel.from( params.minBLS..params.maxBLS )

process filterBEDs {

	module 'python/x86_64/3.5.1'

	publishDir "$params.outdir", pattern: "*.txt", mode: 'copy'

	input:
	tuple val(threshold), file(bedFolder) from blsRange_ch.combine(originalMotifs_ch)

	output:
	tuple val(threshold), file("filtered.bed") into filteredBED_ch

	"""
	#!/usr/bin/env python3
	
	from os import listdir
	from os.path import isfile, join
	import shutil

	folder = "$bedFolder"
	outLoc = "filtered.bed"
	targetGenesLoc = "targetGenes_${threshold}.txt"

	motifFiles = [f for f in listdir(folder) if (isfile(join(folder, f)) and f[-4:]=='.bed')]

	threshold = $threshold

	with open(targetGenesLoc, 'w') as targetGenesOut:
		with open(outLoc, 'w') as outFile:
			for motifFile in motifFiles:
				targetGenes = set()
				with open(join(folder, motifFile)) as bedFile:
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

	input:
	tuple val(threshold), file(filteredBED) from filteredBED_ch

	output:
	tuple val(threshold), env(COVERAGE_ORIGINAL), env(COVERAGE_FULL), env(COVERAGE_INTRONS), env(COVERAGE_NONINTRONS) into gather_ch

	"""
	bedtools intersect -a $filteredBED -b $params.fullGFF > fullSearchSpace.bed
	bedtools intersect -a $filteredBED -b $params.intronGFF > introns.bed
	bedtools intersect -a $filteredBED -b $params.nonIntronGFF > nonIntrons.bed

	COVERAGE_ORIGINAL=\$(cat $filteredBED | sortBed | mergeBed | awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2+1 }END{print SUM}')
	COVERAGE_FULL=\$(cat fullSearchSpace.bed | sortBed | mergeBed | awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2+1 }END{print SUM}')
	COVERAGE_INTRONS=\$(cat introns.bed | sortBed | mergeBed | awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2+1 }END{print SUM}')
	COVERAGE_NONINTRONS=\$(cat nonIntrons.bed | sortBed | mergeBed | awk -F'\\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2+1 }END{print SUM}')

	"""
}

process gatherResults {

	input:
	tuple val(threshold), val(coverage_original), val(coverage_full), val(coverage_introns), val(coverage_nonintrons) from gather_ch

	output:
	tuple val(threshold), val(coverage_original), val(coverage_full), val(coverage_introns), val(coverage_nonintrons) into final_ch

	"""
	
	echo "$threshold	$coverage_original" >> $params.outdir/coverageResults_original.txt
	echo "$threshold	$coverage_full" >> $params.outdir/coverageResults_fullSearch.txt
	echo "$threshold	$coverage_introns" >> $params.outdir/coverageResults_introns.txt
	echo "$threshold	$coverage_nonintrons" >> $params.outdir/coverageResults_nonIntrons.txt

	"""
}

process getTargetGeneStats {

	module 'python/x86_64/3.5.1'

	input:
	file thingsThatNeedToBeCalculated from final_ch.collect()

	"""
	#!/usr/bin/env python3
	
	import os
	os.environ["OMP_NUM_THREADS"] = "1"
	from os import listdir
	from os.path import isfile, join
	from statistics import median, mean
	import numpy

	folder = "$params.outdir"
	medianOut = "$params.outdir/medians.txt"
	meanOut = "$params.outdir/means.txt"
	perc75Out = "$params.outdir/perc75.txt"
	perc90Out = "$params.outdir/perc90.txt"

	files = [f for f in listdir(folder) if (isfile(join(folder, f)) and f[:11]=='targetGenes')]

	with open(medianOut, 'w') as medianFile, open(meanOut, 'w') as meanFile, open(perc75Out, 'w') as perc75File, open(perc90Out, 'w') as perc90File:
		for file in files:
			nbOfTargetGenes = list()
			with open(join(folder, file)) as f:
				for line in f:
					nbOfTargetGenes.append(int(line.rstrip()))
			medianFile.write(file.split(".")[0].split("_")[1]+"\\t"+str(median(nbOfTargetGenes))+"\\n")
			meanFile.write(file.split(".")[0].split("_")[1]+"\\t"+str(mean(nbOfTargetGenes))+"\\n")
			perc75File.write(file.split(".")[0].split("_")[1]+"\\t"+str(numpy.percentile(nbOfTargetGenes, 75))+"\\n")
			perc90File.write(file.split(".")[0].split("_")[1]+"\\t"+str(numpy.percentile(nbOfTargetGenes, 90))+"\\n")
	"""
}
