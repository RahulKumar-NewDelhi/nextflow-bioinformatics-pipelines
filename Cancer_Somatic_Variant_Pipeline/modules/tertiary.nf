#!/usr/bin/env nextflow

// Tertiary analysis against the shared SQLite database (see db/schema.sql).
//   tertiaryLoadSample  loads a SnpEff-annotated VCF into the database
//   tertiaryReport      writes the gene / disease-association report for that run
//
// The database is a large shared file, so it is passed by path and NOT staged into the work
// directory: params.tertiaryDB must be somewhere the container can see (the docker profile mounts
// $HOME, /mnt and /media). Build it once with the db/load_*.py scripts.

params.outputLocation = ""
params.tertiaryDB = ""
params.tumorType = ""
params.minDepth = 10
params.minAlleleFraction = 0.05

process tertiaryLoadSample {

    label 'python3'
    maxForks 1          // one writer at a time: the database is a single SQLite file
    tag "${sampleName}"

    input:
    path vcfFile
    tuple path(ref), path(refIndex)

    output:
    tuple val(sampleName), path("${sampleName}.run_id.txt")

    script:
    sampleName = vcfFile.getSimpleName()
    def tumor = params.tumorType ? "--tumor-type '${params.tumorType}'" : ""
    """
    python3 $projectDir/db/load_sample.py \
        --db ${params.tertiaryDB} \
        --sample $sampleName \
        --source-dir ${params.outputLocation} \
        --vcf $vcfFile \
        --reference $ref \
        --run-id-file ${sampleName}.run_id.txt $tumor
    """
}

process tertiaryReport {

    label 'python3'
    publishDir "${params.outputLocation}/Analysis/Tertiary/$sampleName", mode: 'copy'
    tag "${sampleName}"

    input:
    tuple val(sampleName), path(runIdFile)

    output:
    path "${sampleName}.priority_variants.tsv"
    path "${sampleName}.gene_disease_associations.tsv"
    path "${sampleName}.tertiary_report.html"

    script:
    """
    python3 $projectDir/db/report.py \
        --db ${params.tertiaryDB} \
        --sample $sampleName \
        --run-id \$(cat $runIdFile) \
        --min-dp ${params.minDepth} \
        --min-af ${params.minAlleleFraction} \
        --outdir .
    """
}
