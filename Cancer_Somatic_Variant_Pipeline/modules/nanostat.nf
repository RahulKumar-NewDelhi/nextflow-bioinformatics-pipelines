#!/usr/bin/env nextflow

params.outputLocation = ""
params.cpus = ""
params.outputLocation = ""

process nanostatReadStats {

	label 'nanostat'
	
	input:
	tuple val(meta), path(fastq)

	output:	
	path "${sampleName}.statreports"

	tag "${sampleName}"

	script:
	
	sampleName = meta.id

	"""
	NanoStat --fastq $fastq --outdir ${sampleName}.statreports --threads $params.cpus --name ${sampleName}.nanostat.tsv

	"""
}


process nanoplotReadStats {

	label 'nanoplot'
	publishDir "${params.outputLocation}/Analysis/NanoPlot", mode: 'copy'

	input:
	tuple val(meta), path(fastq)

	output:
	tuple val(meta), path("*")

	tag "${sampleName}"

	script:

	sampleName = meta.id

	"""
	NanoPlot -t $params.cpus --fastq $fastq --maxlength 40000 --plots dot --legacy hex -o ${sampleName}_NanoPlot

	"""
}