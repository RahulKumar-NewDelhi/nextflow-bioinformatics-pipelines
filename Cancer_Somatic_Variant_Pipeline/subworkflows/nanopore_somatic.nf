#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Nanopore counterpart of SOMATIC_ANALYSIS for UMI-tagged AMPLICON panels
// (UMI consensus method of genepi/umi-pipeline-nf, Karst et al. 2021 / Amstler et al. 2024).
// Each BED line must be one amplicon: reads are assigned to an amplicon only if they span
// >= params.ontMinOverlap of it, so hybrid-capture data will not work here.
//
// Offline mode only (no live MinKNOW watching). Tuning parameters live in modules/umi_ont.nf.

params.reference = ""
params.bed = ""
params.snpeffDB = ""
params.tertiaryDB = ""

params.ontSubsample = false
params.ontChunkSize = 50000                 // reads per mapping job
params.ontMinReadsPerTarget = 1000          // amplicons with fewer on-target reads are skipped
params.ontClusterIdentity = 0.90            // vsearch identity for raw UMI clustering
params.ontMinReadsPerCluster = 20           // smaller UMI families are dropped
params.ontMedakaModel = "r1041_e82_400bps_hac_g615"   // must match the basecalling model; name (downloaded once) or a local *_model_pt.tar.gz
params.ontReferenceBasedPolishing = false   // true: faster medaka consensus against the amplicon, no second UMI round

params.cpus = 1
params.outputLocation = ""

include { rawCounts; processedCounts; readStatistics } from '../modules/counts'
include { nanoplotReadStats } from '../modules/nanostat'
include { samtoolsIndexFasta; samtoolsFlagStat } from '../modules/samtools'
include { bedtoolsCoverage } from '../modules/bedtools'
include { snpeffEffectHuman } from '../modules/snpeff'
include { tertiaryLoadSample; tertiaryReport } from '../modules/tertiary'
include {
	ontParseBed; ontTargetReference; ontMinimap2Index; ontFilterReads; ontSubsample; ontMapReads;
	ontSplitReads; ontSummarizeFilterStats;
	ontDetectUmi; ontDetectUmi as ontDetectUmiConsensus; ontSummarizeUmiStats;
	ontCluster; ontCluster as ontClusterConsensus;
	ontReformatFilterCluster; ontClusterStats; ontGlueClusters;
	ontMedakaModel; ontPolishCluster; ontPolishClusterRef; ontMergeMaskConsensus; ontReformatConsensus;
	ontMapConsensus; ontMergeSampleBam; ontLofreq; ontVcfStats
} from '../modules/umi_ont'


workflow NANOPORE_SOMATIC_ANALYSIS {

	take: sampleData
	main:

	if (!params.bed) {
		error "NANOPORE_SOMATIC_ANALYSIS needs params.bed: one line per UMI amplicon (chrom, start, end, name)."
	}

	def umiBin = "${projectDir}/bin/umi_ont"
	filterReadsPy = file("${umiBin}/filter_reads.py", checkIfExists: true)
	summarizeFilterPy = file("${umiBin}/summarize_split_read_stats.py", checkIfExists: true)
	extractUmisPy = file("${umiBin}/extract_umis.py", checkIfExists: true)
	summarizeUmiPy = file("${umiBin}/summarize_umi_stats.py", checkIfExists: true)
	parseClustersPy = file("${umiBin}/parse_clusters.py", checkIfExists: true)
	clusterReportPy = file("${umiBin}/cluster_report.py", checkIfExists: true)
	parseClusterBamPy = file("${umiBin}/parse_cluster_alignment.py", checkIfExists: true)
	reformatConsensusPy = file("${umiBin}/reformat_consensus.py", checkIfExists: true)

	// Read QC. No porechop: the UMIs sit inside the adapter region that it would trim away.
	rawCounts(sampleData)
	nanoplotReadStats(sampleData)

	ontFilterReads(sampleData.map { meta, reads -> tuple(meta.id, reads) })
	processedCounts(ontFilterReads.out.map { sampleName, reads -> tuple([id: sampleName], reads) })
	readStatistics(rawCounts.out.collect(), processedCounts.out.collect())

	if (params.ontSubsample) {
		ontSubsample(ontFilterReads.out)
		filteredReads = ontSubsample.out.reads
	} else {
		filteredReads = ontFilterReads.out
	}

	samtoolsIndexFasta(params.reference)
	ontMinimap2Index(params.reference)
	ontParseBed(params.bed)
	targetBeds = ontParseBed.out.targets.flatten().map { bed -> tuple(bed.getBaseName(), bed) }

	// Map read chunks, then assign full-length reads to each amplicon
	ontMapReads(filteredReads.splitFastq(by: params.ontChunkSize, file: true), ontMinimap2Index.out)
	ontSplitReads(ontMapReads.out.combine(targetBeds), filterReadsPy)
	ontSummarizeFilterStats(ontSplitReads.out.stats.groupTuple(by: [0, 1]), summarizeFilterPy)

	targetReads = ontSplitReads.out.reads
		.groupTuple(by: [0, 1])
		.filter { sampleName, target, fastqs ->
			def n = fastqs.sum { it.countFastq() }
			if (n <= params.ontMinReadsPerTarget) {
				log.warn "${sampleName}:${target} skipped, only ${n} reads span the amplicon (ontMinReadsPerTarget = ${params.ontMinReadsPerTarget})"
			}
			n > params.ontMinReadsPerTarget
		}
		.transpose(by: 2)

	// UMI extraction and clustering into molecule families
	ontDetectUmi(targetReads, 'raw', extractUmisPy)
	ontSummarizeUmiStats(ontDetectUmi.out.stats.groupTuple(by: [0, 1]), 'raw', summarizeUmiPy)
	ontCluster(ontDetectUmi.out.umis.groupTuple(by: [0, 1]), params.ontClusterIdentity)

	families = ontCluster.out.clusters.flatMap { sampleName, target, clusters ->
		def kept = (clusters instanceof List ? clusters : [clusters]).findAll { it.countFasta() >= params.ontMinReadsPerCluster }
		kept ? [tuple(sampleName, target, kept)] : []
	}

	ontReformatFilterCluster(families, parseClustersPy)
	ontClusterStats(ontReformatFilterCluster.out.stats, clusterReportPy)
	ontGlueClusters(ontReformatFilterCluster.out.smolecules)

	// One item per polishing file; the group key carries the file count per sample/amplicon so
	// each amplicon's consensus is merged as soon as its last file is polished.
	gluedFiles = ontGlueClusters.out.flatMap { sampleName, target, glued ->
		def files = glued instanceof List ? glued : [glued]
		files.collect { f -> tuple(groupKey("${sampleName}:${target}", files.size()), sampleName, target, f) }
	}

	if (file(params.ontMedakaModel).exists()) {
		medakaModel = file(params.ontMedakaModel)
	} else {
		ontMedakaModel(params.ontMedakaModel)
		medakaModel = ontMedakaModel.out
	}

	if (params.ontReferenceBasedPolishing) {
		ontTargetReference(targetBeds, samtoolsIndexFasta.out)
		polishInput = gluedFiles
			.map { key, sampleName, target, f -> tuple(target, key, sampleName, f) }
			.combine(ontTargetReference.out, by: 0)
			.map { target, key, sampleName, f, targetRef -> tuple(key, sampleName, target, f, targetRef) }
		ontPolishClusterRef(polishInput, parseClusterBamPy, medakaModel)
		polished = ontPolishClusterRef.out
	} else {
		ontPolishCluster(gluedFiles, medakaModel)
		polished = ontPolishCluster.out
	}

	ontMergeMaskConsensus(polished.groupTuple(by: 0).map { key, samples, targets, fastqs -> tuple(samples[0], targets[0], fastqs) })

	if (params.ontReferenceBasedPolishing) {
		finalConsensus = ontMergeMaskConsensus.out
	} else {
		// Second UMI round on the consensus reads (99% identity) to collapse split families
		ontDetectUmiConsensus(ontMergeMaskConsensus.out, 'consensus', extractUmisPy)
		ontClusterConsensus(ontDetectUmiConsensus.out.umis, 0.99)
		ontReformatConsensus(ontClusterConsensus.out.consensus, reformatConsensusPy)
		finalConsensus = ontReformatConsensus.out
	}

	// One consensus BAM per sample (all amplicons)
	ontMapConsensus(finalConsensus, ontMinimap2Index.out)
	ontMergeSampleBam(ontMapConsensus.out.groupTuple(by: 0))
	samtoolsFlagStat(ontMergeSampleBam.out.map { bam, bai -> bam })
	bedtoolsCoverage(ontParseBed.out.panel, ontMergeSampleBam.out)

	ontLofreq(ontMergeSampleBam.out, samtoolsIndexFasta.out, ontParseBed.out.panel)
	ontVcfStats(ontLofreq.out)
	snpeffEffectHuman(ontLofreq.out)

	if (params.tertiaryDB) {
		tertiaryLoadSample(snpeffEffectHuman.out, samtoolsIndexFasta.out)
		tertiaryReport(tertiaryLoadSample.out)
	}
}
