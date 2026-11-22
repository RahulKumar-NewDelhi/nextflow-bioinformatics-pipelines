#!/usr/bin/env nextflow

params.cpus = 1
params.reference = ""
params.outputLocation = ""

process bcftoolsStats {

    label 'bcftools'
    publishDir "${params.outputLocation}/Analysis/Variants/$sampleName/Mutect2", mode: 'copy'

    input:
    path vcfFile

    tag "${sampleName}"

    output:
    path "${sampleName}.bcftools_stats.txt"

    tag "${sampleName}"

    script:
    sampleName = vcfFile.getSimpleName()
    """
    bcftools stats -s - $vcfFile > ${sampleName}.bcftools_stats.txt
    """
}
