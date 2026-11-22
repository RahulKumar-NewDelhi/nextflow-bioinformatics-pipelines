#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.platform = ""
params.tool = ""
params.sampleSheet = ""
params.cpus = 1
params.outputLocation = ""

params.length = 20
params.quality = 30
params.reference = ""
params.bed = ""
params.snpeffDB = ""
params.clinvar = ""
params.clinvarDB = ""
params.tertiaryDB = ""
params.tumorType = ""

include { validateSampleSheet; createSamplesChannel } from './modules/utils'
include { SOMATIC_ANALYSIS } from './subworkflows/somatic'
include { NANOPORE_SOMATIC_ANALYSIS } from './subworkflows/nanopore_somatic'

workflow {
    isValidSamplesheet = validateSampleSheet(params.sampleSheet)
    if (isValidSamplesheet) {
        ch_samples = Channel.from(createSamplesChannel(params.sampleSheet))
        if (params.tool != "somatic")
            error "Pipeline is not available: only tool = 'somatic' is supported"
        if (params.platform == "nanopore") {
            NANOPORE_SOMATIC_ANALYSIS(ch_samples)
        } else if (params.platform == "illumina" || params.platform == "ion-torrent") {
            SOMATIC_ANALYSIS(ch_samples)
        } else {
            error "Pipeline is not available for platform '${params.platform}'"
        }
    } else {
        error "ERROR: Not a valid input in sampleSheet. Input should contain fastq files."
    }
}
