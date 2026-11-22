#!/usr/bin/env nextflow

params.cpus = 1
params.outputLocation = ""
params.snpeffDB = ""

process snpeffEffectHuman {
    
    label 'snpeff'
    publishDir "${params.outputLocation}/Analysis/Variants/Annotation/$sampleName", mode: 'copy'

    input:
    path vcfFile

    output:
    path "${sampleName}.snpEffAnn.vcf"

    tag "${sampleName}"

    script:
    sampleName = vcfFile.getSimpleName()

    """
    snpeff -Xmx8g eff ${params.snpeffDB} $vcfFile > ${sampleName}.snpEffAnn.vcf

    """
}

