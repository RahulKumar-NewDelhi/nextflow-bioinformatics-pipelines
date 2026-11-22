#!/usr/bin/env nextflow

params.echo = true
params.outputLocation = ""


process trimgaloreTrimReads {
    
    label 'trimgalore'
    echo params.echo
    publishDir "${params.outputLocation}/Analysis/Processed_Reads/$sampleName", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    val length
    val quality

    tag "${sampleName}"

    output:
    tuple val(meta), path("*.fq")
    path "*fastqc.html"
    path "*trimming_report.txt"

    tag "${sampleName}"

    script:
    sampleName = meta.id
    if (!meta.pair_end) 
        """
        trim_galore $reads --length $length --quality $quality --fastqc -j 4
        """
    else 
        """
        trim_galore $reads --length $length --quality $quality --paired --fastqc -j 4
        """
}
