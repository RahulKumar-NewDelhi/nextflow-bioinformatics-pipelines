#!/usr/bin/env nextflow

params.cpus = 1
params.reference = ""
params.outputLocation = ""


process bwaIndex {

	label 'bwa'
	publishDir "${refPath}", mode: 'copy'

	input:
	path reference

	output:
	path "${reference}*"

	script:

	refPath = file(params.reference).getParent()

	"""
	bwa index $reference

	"""

}


process bwaAlignment {

    label 'bwa'
    cpus params.cpus

    input:
    tuple val(meta), path(reads)
    path bwaIndex

    tag "${sampleName}"

    tag "${sampleName}"

    output:
    path "${sampleName}.sam"

    tag "${sampleName}"

    script:

    sampleName = meta.id
    refPrefix = file(params.reference).getName()

    """
    bwa mem -M -t $task.cpus -R "@RG\tID:$sampleName\tPL:ILLUMINA\tSM:$sampleName" $refPrefix $reads > "${sampleName}.sam"
    """
}
