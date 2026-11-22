#!/usr/bin/env nextflow

params.echo = true
params.outputLocation = ""
params.bed = ""
params.reference = ""

process bedtoolsCoverage {

    label 'bedtools'
    publishDir "${params.outputLocation}/Analysis/Statistics/Horizontal_Coverage", mode: 'copy'

    input:
    path bed
	tuple path(sortedBam), path(sortedBamIndex)

    output:
    path("${sampleName}.coverage.tsv")

    tag "${sampleName}"

    script:
    sampleName =  sortedBam.getSimpleName()
    """
    bedtools coverage -d -a $bed -b $sortedBam > ${sampleName}.coverage.tsv
    """

}

