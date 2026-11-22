#!/usr/bin/env nextflow

params.outputLocation = ""

process generateReport {
    
    label 'rmarkdown'
    publishDir "${params.outputLocation}/Reports", mode: 'copy'

    input:
    tuple val(sampleName), path(statsFile)

    tag "${sampleName}"

    output:
    path "${sampleName}.html"

    tag "${sampleName}"

    script:
    sampleName = statsFile.getSimpleName()
    """
    cp -v $projectDir/data/sample_report.rmd .
    Rscript $projectDir/bin/render_report.r $sampleName sample_report.rmd
    """
}
