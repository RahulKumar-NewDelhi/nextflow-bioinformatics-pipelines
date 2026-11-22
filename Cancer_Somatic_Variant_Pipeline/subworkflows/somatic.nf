#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.length = 20
params.quality = 30
params.reference = ""
params.dbsnp = ""
params.snpeffDB = ""
params.tertiaryDB = ""

// Matched tumor-normal Mutect2 (additional mode alongside the default
// tumor-only path). Needs the sample sheet's optional patient_id/sample_type
// columns (see modules/utils.nf createTumorNormalPairs) to pair samples.
params.matchedNormal = false
params.germlineResource = ""   // required if matchedNormal is true, e.g. af-only-gnomad.hg38.vcf.gz
params.pon = ""                 // optional panel-of-normals vcf.gz; omit to run without one

// UMI consensus mode (fgbio). When enabled, MarkDuplicates and BQSR are replaced by
// UMI grouping + consensus calling; Mutect2 then runs on the consensus BAM.
// Defaults describe an inline duplex UMI kit: 5 bp UMI + 2 bp spacer at the start of R1 and R2.
params.umi = false
params.umiReadStructure = "5M2S+T 5M2S+T"   // fgbio read structures, one per FASTQ
params.umiDuplex = true                      // true: CallDuplexConsensusReads, false: CallMolecularConsensusReads
params.umiEdits = 1                          // UMI mismatches allowed when grouping
params.umiMinMapQ = 20                       // min mapping quality for reads used in grouping
params.umiConsensusMinReads = 1              // min reads to build a consensus (before filtering)
params.umiFilterMinReads = ""                // FilterConsensusReads --min-reads; empty = "2 1 1" (duplex) or "3" (simplex)
params.umiMinBaseQuality = 30                // consensus bases below this are masked to N
params.umiMaxBaseErrorRate = 0.2             // mask a consensus base if more than this fraction of its reads disagree (fgbio default 0.1)

params.cpus = 1
params.outputLocation = ""

include { rawCounts; processedCounts; readStatistics} from "../modules/counts"
include { trimgaloreTrimReads } from "../modules/trimgalore"
include { bwaIndex; bwaAlignment } from "../modules/bwa"
include { samtoolsSamToBam; samtoolsSortBam; samtoolsFlagStat; samtoolsIndexFasta; samtoolsIndexBam } from "../modules/samtools"
include { gatkMarkDuplicates; gatkAddReadGroups; gatkCreateSequenceDictionary; gatkHaplotypeCaller; gatkBaseRecalibratorRef; gatkApplyBQSRRef; gatkMutect2; gatkMutect2TN; gatkFilterMutectCalls; gatkVariantsToTableMutect2 } from '../modules/gatk'
include { bcftoolsStats } from '../modules/bcftools'
include { snpeffEffectHuman } from '../modules/snpeff'
include { tertiaryLoadSample; tertiaryReport } from '../modules/tertiary'
include { createTumorNormalPairs } from '../modules/utils'


// ---------------------------------------------------------------------------
// UMI processes (fgbio best-practice consensus pipeline)
// ---------------------------------------------------------------------------

// FASTQ -> unmapped BAM, UMI moved from the read sequence into the RX tag
process umiFastqToBam {

	container 'quay.io/biocontainers/fgbio:2.5.21--hdfd78af_0'
	cpus params.cpus

	input:
	tuple val(meta), path(reads)

	output:
	path "${sampleName}.unmapped.bam"

	tag "${sampleName}"

	script:
	sampleName = meta.id
	"""
	fgbio -Xmx4g --compression 1 --async-io FastqToBam \
	--input $reads \
	--read-structures ${params.umiReadStructure} \
	--sample ${sampleName} \
	--library ${sampleName} \
	--read-group-id ${sampleName} \
	--platform illumina \
	--output ${sampleName}.unmapped.bam
	"""
}

// Align raw reads; read order is preserved (-K) so ZipperBams can merge the tags back
process umiAlignRaw {

	label 'samtools'
	cpus params.cpus

	input:
	path unmappedBam
	path bwaIndex

	output:
	tuple path(unmappedBam), path("${sampleName}.raw.mapped.bam")

	tag "${sampleName}"

	script:
	sampleName = unmappedBam.getSimpleName()
	refPrefix = file(params.reference).getName()
	"""
	set -o pipefail
	samtools fastq $unmappedBam \
	| bwa mem -p -K 150000000 -Y -t $task.cpus $refPrefix - \
	| samtools view -b -o ${sampleName}.raw.mapped.bam -
	"""
}

// Restore UMI tags onto the aligned reads, then group reads by UMI (replaces MarkDuplicates)
process umiGroupReads {

	container 'quay.io/biocontainers/fgbio:2.5.21--hdfd78af_0'
	cpus params.cpus
	publishDir "${params.outputLocation}/Analysis/Alignment/$sampleName/UMI", mode: 'copy', pattern: "*.txt"

	input:
	tuple path(unmappedBam), path(mappedBam)
	tuple path(ref), path(refIndex), path(ref_dict)

	output:
	path "${sampleName}.grouped.bam", emit: bam
	path "${sampleName}.*.txt", emit: metrics

	tag "${sampleName}"

	script:
	sampleName = unmappedBam.getSimpleName()
	strategy = params.umiDuplex ? "paired" : "adjacency"
	duplexMetrics = params.umiDuplex ? "fgbio -Xmx4g CollectDuplexSeqMetrics --input ${sampleName}.grouped.bam --output ${sampleName}.duplex_seq_metrics" : ""
	"""
	fgbio -Xmx4g --compression 1 ZipperBams \
	--unmapped $unmappedBam \
	--input $mappedBam \
	--ref $ref \
	--output ${sampleName}.zipped.bam

	fgbio -Xmx4g --compression 1 GroupReadsByUmi \
	--input ${sampleName}.zipped.bam \
	--output ${sampleName}.grouped.bam \
	--strategy $strategy \
	--edits ${params.umiEdits} \
	--min-map-q ${params.umiMinMapQ} \
	--family-size-histogram ${sampleName}.family_sizes.txt

	$duplexMetrics
	rm ${sampleName}.zipped.bam
	"""
}

// One consensus read per UMI family (unmapped BAM)
process umiCallConsensus {

	container 'quay.io/biocontainers/fgbio:2.5.21--hdfd78af_0'
	cpus params.cpus

	input:
	path groupedBam

	output:
	path "${sampleName}.cons.unmapped.bam"

	tag "${sampleName}"

	script:
	sampleName = groupedBam.getSimpleName()
	caller = params.umiDuplex ? "CallDuplexConsensusReads" : "CallMolecularConsensusReads"
	"""
	fgbio -Xmx4g --compression 1 $caller \
	--input $groupedBam \
	--output ${sampleName}.cons.unmapped.bam \
	--min-reads ${params.umiConsensusMinReads} \
	--read-group-id ${sampleName} \
	--threads $task.cpus
	"""
}

// Align the consensus reads
process umiAlignConsensus {

	label 'samtools'
	cpus params.cpus

	input:
	path consensusBam
	path bwaIndex

	output:
	tuple path(consensusBam), path("${sampleName}.cons.mapped.bam")

	tag "${sampleName}"

	script:
	sampleName = consensusBam.getSimpleName()
	refPrefix = file(params.reference).getName()
	"""
	set -o pipefail
	samtools fastq $consensusBam \
	| bwa mem -p -K 150000000 -Y -t $task.cpus $refPrefix - \
	| samtools view -b -o ${sampleName}.cons.mapped.bam -
	"""
}

// Restore consensus tags, filter/mask weak consensus reads, clip overlapping mates
process umiFilterConsensus {

	container 'quay.io/biocontainers/fgbio:2.5.21--hdfd78af_0'
	cpus params.cpus
	publishDir "${params.outputLocation}/Analysis/Alignment/$sampleName", mode: 'copy'

	input:
	tuple path(consensusBam), path(mappedBam)
	tuple path(ref), path(refIndex), path(ref_dict)

	output:
	path "${sampleName}.consensus.bam"

	tag "${sampleName}"

	script:
	sampleName = consensusBam.getSimpleName()
	minReads = params.umiFilterMinReads ?: (params.umiDuplex ? "2 1 1" : "3")
	"""
	fgbio -Xmx4g --compression 1 ZipperBams \
	--unmapped $consensusBam \
	--input $mappedBam \
	--ref $ref \
	--tags-to-reverse Consensus \
	--tags-to-revcomp Consensus \
	--output ${sampleName}.cons.zipped.bam

	fgbio -Xmx4g --compression 1 FilterConsensusReads \
	--input ${sampleName}.cons.zipped.bam \
	--output ${sampleName}.cons.filtered.bam \
	--ref $ref \
	--min-reads $minReads \
	--min-base-quality ${params.umiMinBaseQuality} \
	--max-base-error-rate ${params.umiMaxBaseErrorRate}

	fgbio -Xmx4g ClipBam \
	--input ${sampleName}.cons.filtered.bam \
	--output ${sampleName}.consensus.bam \
	--ref $ref \
	--clipping-mode Hard \
	--clip-overlapping-reads true \
	--sort-order Coordinate

	rm ${sampleName}.cons.zipped.bam ${sampleName}.cons.filtered.bam
	"""
}


workflow SOMATIC_ANALYSIS {

    take: sampleData
    main:

    if (params.matchedNormal && !params.germlineResource) {
        error "params.matchedNormal is true but params.germlineResource is not set (path to an af-only-gnomad vcf.gz)."
    }

    rawCounts(sampleData)
    trimgaloreTrimReads(sampleData, params.length, params.quality)
    processedCounts(trimgaloreTrimReads.out[0])
    readStatistics(rawCounts.out.collect(), processedCounts.out.collect())

    refPath = file(params.reference).getParent()
	refFile = file(params.reference)
	refName = refFile.getBaseName()
	bwaIndexExt = ["amb", "ann", "bwt", "pac", "sa"]
	def checkIndex = []
	for (int i=0; i<bwaIndexExt.size(); i++) {
	    indexFile = refFile + "." + bwaIndexExt[i]
	    if (indexFile.exists()) {
		    checkIndex.add(true)
		}
	}
	if (checkIndex.size != 5) {
		bwaIndex(params.reference)
		bwaIndexCh = bwaIndex.out
	} else {
		def index_files = [:]
		refPath.eachFile { it ->
		    if (it.name.startsWith(refName) && bwaIndexExt.contains(file(it).getExtension())) {
			    index_files[it.name] = it
		    }
	    }
	    index_ch = Channel.from(index_files.values())
		bwaIndexCh = index_ch.collect()
	}

	samtoolsIndexFasta(params.reference)
	gatkCreateSequenceDictionary(samtoolsIndexFasta.out)
	refCh = gatkCreateSequenceDictionary.out

	if (params.umi) {
		// UMIs are read from the untrimmed FASTQs (trimming would shift inline UMIs);
		// trimgalore above still runs for the read statistics / FastQC report.
		// Consensus base qualities are recomputed by fgbio, so no MarkDuplicates or BQSR.
		umiFastqToBam(sampleData)
		umiAlignRaw(umiFastqToBam.out, bwaIndexCh)
		umiGroupReads(umiAlignRaw.out, refCh)
		umiCallConsensus(umiGroupReads.out.bam)
		umiAlignConsensus(umiCallConsensus.out, bwaIndexCh)
		umiFilterConsensus(umiAlignConsensus.out, refCh)
		samtoolsFlagStat(umiFilterConsensus.out)
		samtoolsIndexBam(umiFilterConsensus.out)
	} else {
		bwaAlignment(trimgaloreTrimReads.out[0], bwaIndexCh)
		samtoolsSamToBam(bwaAlignment.out)
		samtoolsSortBam(samtoolsSamToBam.out)
		gatkMarkDuplicates(samtoolsSortBam.out)
		samtoolsFlagStat(gatkMarkDuplicates.out)

		gatkBaseRecalibratorRef(gatkMarkDuplicates.out, refCh, params.dbsnp)
		gatkApplyBQSRRef(gatkMarkDuplicates.out, gatkBaseRecalibratorRef.out, refCh)

		samtoolsIndexBam(gatkApplyBQSRRef.out)
	}

	if (params.matchedNormal) {
		// every sample (tumor + normal) still goes through alignment/BQSR above,
		// individually. Pairing only happens here, right before Mutect2.
		pairs = Channel.from(createTumorNormalPairs(params.sampleSheet))

		bamById = samtoolsIndexBam.out.map { bam, bai -> tuple(bam.getSimpleName(), bam, bai) }

		tnInput = pairs
			.map { patientId, tumorId, normalId -> tuple(tumorId, patientId, normalId) }
			.join(bamById)                                     // -> tumorId, patientId, normalId, tumorBam, tumorBai
			.map { tumorId, patientId, normalId, tumorBam, tumorBai ->
				tuple(normalId, patientId, tumorId, tumorBam, tumorBai)
			}
			.join(bamById)                                     // -> normalId, patientId, tumorId, tumorBam, tumorBai, normalBam, normalBai
			.map { normalId, patientId, tumorId, tumorBam, tumorBai, normalBam, normalBai ->
				tuple(patientId, tumorId, tumorBam, tumorBai, normalId, normalBam, normalBai)
			}

		germlineResourceCh = file(params.germlineResource)
		ponCh = params.pon ? file(params.pon) : file("${projectDir}/assets/NO_FILE")

		gatkMutect2TN(tnInput, refCh, germlineResourceCh, ponCh)
		mutectVcf = gatkMutect2TN.out
	} else {
		gatkMutect2(samtoolsIndexBam.out, refCh)
		mutectVcf = gatkMutect2.out
	}

	gatkFilterMutectCalls(mutectVcf, refCh)
	gatkVariantsToTableMutect2(gatkFilterMutectCalls.out, refCh)
	bcftoolsStats(gatkFilterMutectCalls.out)
	snpeffEffectHuman(gatkFilterMutectCalls.out)

	if (params.tertiaryDB) {
		tertiaryLoadSample(snpeffEffectHuman.out, samtoolsIndexFasta.out)
		tertiaryReport(tertiaryLoadSample.out)
	}


}
