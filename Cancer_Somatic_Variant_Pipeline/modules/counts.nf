#!/usr/bin/env nextflow

params.outputLocation = ""


process rawCounts {

	label 'python3'

	input:
	tuple val(meta), path(reads)

	tag "${sampleName}"

	output:
	path("${sampleName}.RawCount.txt")

	tag "${sampleName}"

	script:

	sampleName = meta.id
	compressed = reads.collect { it.getExtension() }.contains("gz")

	if(compressed)
		"""
		count=\$(zcat $reads | echo \$((`wc -l`/4)))
		echo -e ${sampleName}'\\t'\$count > ${sampleName}.RawCount.txt
		"""
	else
		"""
		count=\$(cat $reads | echo \$((`wc -l`/4)))
		echo -e ${sampleName}'\\t'\$count > ${sampleName}.RawCount.txt
		"""
}



process processedCounts {

	label 'python3'

	input:
	tuple val(meta), path(reads)

	tag "${sampleName}"

	output:
	path("${sampleName}.ProcessedCount.txt")

	tag "${sampleName}"

	script:
	sampleName = meta.id
	compressed = reads.collect { it.getExtension() }.contains("gz")

	if(compressed)
		"""
		count=\$(zcat $reads | echo \$((`wc -l`/4)))
		echo -e ${sampleName}'\\t'\$count > ${sampleName}.ProcessedCount.txt
		"""
	else
		"""
		count=\$(cat $reads | echo \$((`wc -l`/4)))
		echo -e ${sampleName}'\\t'\$count > ${sampleName}.ProcessedCount.txt
		"""
}


process readStatistics {

	label 'python3'
	publishDir "${params.outputLocation}/Analysis/Read_Statistics", mode: 'copy'

	input:
	path rawReads
	path processedReads

	output:
	path "read_statistics.tsv"

	script:

	"""
	cat $rawReads > raw_read_stats.tsv
	cat $processedReads > processed_read_stats.tsv

	sed -i '1i sample name\traw read count' raw_read_stats.tsv
	sed -i '1i sample name\tprocessed read count' processed_read_stats.tsv

	python3 $projectDir/bin/join_by_first_column.py raw_read_stats.tsv processed_read_stats.tsv > read_statistics.tsv
	"""

}
