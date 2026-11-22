#!/usr/bin/env nextflow

params.cpus = 1
params.reference = ""
params.outputLocation = ""

process samtoolsSamToBam {
    
    label 'samtools'

    input:
    path(samFile)

    output:
    path("${sampleName}.bam")

    tag "${sampleName}"

    tag "${sampleName}"

    script:
    sampleName = samFile.getSimpleName()
    """
    samtools view -S -b $samFile -o ${sampleName}.bam -@ $task.cpus
    """
}

process samtoolsSortBam {

    label 'samtools'
    cpus params.cpus
    publishDir "${params.outputLocation}/Analysis/Alignment/$sampleName", mode: 'copy'

    input:
    path(bamFile)

    tag "${sampleName}"

    output:
    path("${sampleName}.sorted.bam")

    tag "${sampleName}"

    script:
    sampleName = bamFile.getSimpleName()
    """
    samtools sort $bamFile -o ${sampleName}.sorted.bam -@ $task.cpus
    """
}

process samtoolsIndexFasta {
    
    label 'samtools'
    cpus params.cpus

    input:
    path ref

    output:
    tuple path(ref), path("${ref}.fai")

    script:
    """
    samtools faidx $ref
    """
}

process samtoolsIndexBam {
    
    label 'samtools'
    cpus params.cpus
    publishDir "${params.outputLocation}/Analysis/Alignment/$sampleName", mode: 'copy'

    input:
    path bamFile

    tag "${sampleName}"

    output:
    tuple path(bamFile), path("${sampleName}*.bam.bai")

    tag "${sampleName}"

    script:
    sampleName = bamFile.getSimpleName()
    """
    samtools index $bamFile
    """
}

process samtoolsFlagStat {

    label 'samtools'
    publishDir "${params.outputLocation}/Analysis/Alignment/$sampleName", mode: 'copy'

    input:
    path bamFile

    tag "${sampleName}"

    output:
    path "${sampleName}.flagstat.tsv"

    tag "${sampleName}"

    script:

    sampleName = bamFile.getSimpleName()

    """
    samtools flagstat $bamFile -O tsv > ${sampleName}.flagstat.tsv
    """
}

