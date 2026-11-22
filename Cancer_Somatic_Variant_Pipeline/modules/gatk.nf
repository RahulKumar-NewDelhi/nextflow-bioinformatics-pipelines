#!/usr/bin/env nextflow

params.cpus = 1
params.reference = ""
params.platform = ""
params.outputLocation = ""
params.bed = ""
params.dbsnp = ""

process gatkMarkDuplicates {

    label 'gatk'
    cpus params.cpus
    publishDir "${params.outputLocation}/Analysis/Alignment/$sampleName", mode: 'copy'

    input:
    path bamFile

    tag "${sampleName}"

    output:
    path "${sampleName}.markdup.bam"

    tag "${sampleName}"

    script:
    sampleName = bamFile.getSimpleName()
    """
    gatk --java-options "-XX:ParallelGCThreads=1" MarkDuplicates \
    --REMOVE_SEQUENCING_DUPLICATES true \
    -I $bamFile -O ${sampleName}.markdup.bam -M ${sampleName}.metrics.txt
    """
}

process gatkAddReadGroups {

    label 'gatk'
    publishDir "${params.outputLocation}/Analysis/Alignment/$sampleName", mode: 'copy'

    input:
    path bamFile

    tag "${sampleName}"

    output:
    path "${sampleName}.markdup.rg.bam"

    tag "${sampleName}"

    script:
    sampleName = bamFile.getSimpleName()
    """
    READ_ID=`gatk ViewSam -I $bamFile  --ALIGNMENT_STATUS All --PF_STATUS All --RECORDS_ONLY 2> /dev/null|head -n 1|cut -f1`
    ID=`echo \$READ_ID | tr ':' '\t' | cut -f1`
    PU=`echo \$READ_ID | tr ':' '\t' | cut -f3,4 | tr '\t' '.'`
    gatk AddOrReplaceReadGroups \
    -I $bamFile -O ${sampleName}.markdup.rg.bam -ID $ID -LB LIB-${sampleName} -PL ${params.platform} -PU $PU -SM ${sampleName}
    """
}

process gatkBaseRecalibratorRef {

    label 'gatk'
    cpus params.cpus
    publishDir "${params.outputLocation}/Analysis/Alignment/$sampleName", mode: 'copy'

    input:
    path bamFile
    tuple path(ref), path(refIndex), path(refDict)
    path dbsnp

    tag "${sampleName}"

    output:
    path "${sampleName}.markdup.bqsr.report"

    tag "${sampleName}"

    script:
    sampleName = bamFile.getSimpleName()
    """
    gatk IndexFeatureFile -I $dbsnp
    gatk --java-options "-Xms4G -Xmx4G" BaseRecalibrator \
    -I $bamFile \
    -R $ref \
    --known-sites $dbsnp \
    -O ${sampleName}.markdup.bqsr.report
    """
}

process gatkApplyBQSRRef {

    label 'gatk'
    cpus params.cpus
    publishDir "${params.outputLocation}/Analysis/Alignment/$sampleName", mode: 'copy'

    input:
    path bamFile
    path bqsrReport
    tuple path(ref), path(refIndex), path(refDict)

    tag "${sampleName}"

    output:
    path "${sampleName}.markdup.bqsr.bam"

    tag "${sampleName}"

    script:
    sampleName = bamFile.getSimpleName()
    """
    gatk --java-options "-Xms4G -Xmx4G" ApplyBQSR \
    -I $bamFile \
    -R $ref \
    --bqsr-recal-file $bqsrReport \
    -O ${sampleName}.markdup.bqsr.bam
    """
}

process gatkHaplotypeCaller {

    label 'gatk'
    publishDir "${params.outputLocation}/Analysis/Variants/$sampleName/HaplotypeCaller", mode: 'copy'

    input:
    tuple path(bamFile), path(bamIndex)
    tuple path(ref), path(refIndex), path(ref_dict)

    tag "${sampleName}"

    output:
    path "${sampleName}.g.vcf.gz"

    tag "${sampleName}"

    script:
    sampleName = bamFile.getSimpleName()

    """
    gatk HaplotypeCaller \
    -R $ref \
    -I $bamFile \
    -O ${sampleName}.g.vcf.gz \
    --native-pair-hmm-threads ${params.cpus}

    """
}

process gatkMutect2 {

    label 'gatk'
    publishDir "${params.outputLocation}/Analysis/Variants/$sampleName/Mutect2", mode: 'copy'

    input:
    tuple path(bamFile), path(bamIndex)
    tuple path(ref), path(refIndex), path(ref_dict)

    tag "${sampleName}"

    output:
    tuple path("${sampleName}.vcf.gz"), path("${sampleName}.vcf.gz.stats")

    tag "${sampleName}"

    script:
    sampleName = bamFile.getSimpleName()

    """
    gatk Mutect2 \
    -R $ref \
    -I $bamFile \
    -O ${sampleName}.vcf.gz
    """
}

process gatkMutect2TN {

    label 'gatk'
    publishDir "${params.outputLocation}/Analysis/Variants/$patientId/Mutect2", mode: 'copy'

    input:
    tuple val(patientId), val(tumorId), path(tumorBam), path(tumorBai), val(normalId), path(normalBam), path(normalBai)
    tuple path(ref), path(refIndex), path(ref_dict)
    path germlineResource
    path pon

    tag "${patientId}"

    output:
    tuple path("${patientId}.vcf.gz"), path("${patientId}.vcf.gz.stats")

    tag "${patientId}"

    script:
    havePon = pon.name != "NO_FILE"
    ponIndexCmd = havePon ? "gatk IndexFeatureFile -I $pon" : ""
    ponArg = havePon ? "--panel-of-normals $pon" : ""
    """
    gatk IndexFeatureFile -I $germlineResource
    $ponIndexCmd
    gatk Mutect2 \
    -R $ref \
    -I $tumorBam \
    -I $normalBam \
    -normal $normalId \
    --germline-resource $germlineResource \
    $ponArg \
    -O ${patientId}.vcf.gz
    """
}

process gatkFilterMutectCalls {

    label 'gatk'
    publishDir "${params.outputLocation}/Analysis/Variants/$sampleName/Mutect2", mode: 'copy'

    input:
    tuple path(vcfFile), path(statsFile)
    tuple path(ref), path(refIndex), path(ref_dict)

    tag "${sampleName}"

    output:
    path "${sampleName}.filtered.vcf.gz"

    tag "${sampleName}"

    script:
    sampleName = vcfFile.getSimpleName()

    """
    gatk IndexFeatureFile -I $vcfFile
    gatk FilterMutectCalls \
    -R $ref \
    -V $vcfFile \
    -O ${sampleName}.filtered.vcf.gz
    """
}

process gatkVariantsToTableMutect2 {

    label 'gatk'
    publishDir "${params.outputLocation}/Analysis/Variants/$sampleName/Mutect2", mode: 'copy'

    input:
    path vcfFile
    tuple path(ref), path(refIndex), path(ref_dict)

    tag "${sampleName}"

    output:
    path "${sampleName}.varfilter.pass.tsv"

    tag "${sampleName}"

    script:
    sampleName = vcfFile.getSimpleName()
    """
    gatk IndexFeatureFile -I $vcfFile
    gatk VariantsToTable \
    -F CHROM -F POS -F FILTER -F TYPE -GF AD -GF DP -GF GT \
    --show-filtered \
    -R $ref \
    -V $vcfFile \
    -O ${sampleName}.varfilter.pass.tsv
    """
}

process gatkCreateSequenceDictionary {
    
    label 'gatk'
    cpus params.cpus

    input:
    tuple path(ref), path(refIndex)

    tag "${sampleName}"

    output:
    tuple path(ref), path(refIndex), path("${sampleName}.dict")

    tag "${sampleName}"

    script:
    sampleName = ref.getBaseName()
    """
    gatk CreateSequenceDictionary -R $ref -O ${sampleName}.dict
    """
}
