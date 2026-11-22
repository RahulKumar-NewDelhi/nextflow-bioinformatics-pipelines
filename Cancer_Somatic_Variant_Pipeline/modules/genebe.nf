#!/usr/bin/env nextflow

params.outputLocation = ""
params.genebe_username = ""
params.genebe_apikey = ""
params.cpus = 1

process runGeneBe {
    label 'genebe'
    publishDir "${params.outputLocation}/Analysis/Variants/GeneBe", mode: 'copy'

    input:
    path input_vcf
    path tbi_file

    output:
    path "${input_vcf.simpleName}_genebe.vcf.gz"

    script:
    """
    genebe annotate \\
      --input $input_vcf \\
      --output ${input_vcf.simpleName}_genebe.vcf.gz \\
      --username ${params.genebe_username} \\
      --api_key ${params.genebe_apikey} \\
      --genome hg19 \\
      --batch_size 1000 \\
      --progress
    """
}
