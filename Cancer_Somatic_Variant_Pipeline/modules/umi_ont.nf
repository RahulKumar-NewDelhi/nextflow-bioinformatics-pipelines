#!/usr/bin/env nextflow

// Nanopore UMI amplicon consensus processes, adapted from genepi/umi-pipeline-nf v1.0.1 (MPL-2.0).
// The python helpers are the unmodified upstream scripts in bin/umi_ont/.
// Every process runs in the umi-pipeline-nf container (label 'umi_ont': catfishq, seqtk,
// minimap2, samtools, vsearch, medaka, lofreq, bcftools, pysam).

params.cpus = 1
params.outputLocation = ""

// Read filtering / subsampling (catfishq, seqtk)
params.ontMinReadLength = 0
params.ontMinQscore = 0
params.ontSubsampleReads = 100000
params.ontSubsampleSeed = 11

// Mapping (minimap2)
params.ontMinimap2Kmer = 13

// Read-to-amplicon assignment (filter_reads.py)
params.ontMinOverlap = 0.95          // read must align over >= 95% of the amplicon
params.ontAdapterLength = 100        // bases at each read end searched for the UMI
params.ontIncludeSecondary = false

// UMI extraction (extract_umis.py) and clustering (vsearch)
params.ontFwdUmi = "TTTVVVVTTVVVVTTVVVVTTVVVVTTT"
params.ontRevUmi = "AAABBBBAABBBBAABBBBAABBBBAAA"
params.ontFwdContext = "GTATCGTGTAGAGACTGCGTAGG"
params.ontRevContext = "AGTGATCGAGTCAGTGCGAGTG"
params.ontUseContext = false
params.ontUmiErrors = 2
params.ontUmiMinLength = 40          // combined fwd+rev UMI length window for vsearch
params.ontUmiMaxLength = 60

// Cluster filtering (parse_clusters.py)
params.ontMaxDistUmi = 2
params.ontMaxReadsPerCluster = 60
params.ontClusterFilter = "quality"
params.ontBalanceStrands = true

// Polishing (medaka)
params.ontClustersPerPolishingFile = 100
params.ontPolishMinLength = 50
params.ontPolishMinDepth = 2
params.ontMedakaBatchSize = 50000
params.ontMinConsensusQuality = 40
params.ontMaskingStrategy = "softmask"   // "hardmask" turns low-quality consensus bases into N

// Variant calling (LoFreq)
params.ontMinDepth = 5
params.ontMinVaf = 0.005

def umiOntOut(sample, target) { "${params.outputLocation}/Analysis/UMI_Nanopore/${sample}/${target}" }


// Per-amplicon BED files (one per BED line) plus a cleaned panel BED.
// Header lines are dropped; a missing name column becomes chrom_start_end.
process ontParseBed {

	label 'umi_ont'

	input:
	path bed

	output:
	path "targets/*.bed", emit: targets
	path "panel.bed", emit: panel

	script:
	"""
	mkdir targets
	awk 'BEGIN { OFS = "\\t" }
	/^(track|browser|#)/ || NF < 3 { next }
	{
		name = (NF >= 4 && \$4 != "") ? \$4 : \$1 "_" \$2 "_" \$3
		gsub(/[^A-Za-z0-9._-]/, "_", name)
		print \$1, \$2, \$3, name > ("targets/" name ".bed")
		print \$1, \$2, \$3, name > "panel.bed"
	}' $bed
	"""
}

// Small reference per amplicon (+/- 200 bp), used only by reference-based polishing so that
// medaka's mini_align never has to index the whole genome.
process ontTargetReference {

	label 'umi_ont'
	tag "${target}"

	input:
	tuple val(target), path(bed)
	tuple path(ref), path(refIndex)

	output:
	tuple val(target), path("${target}.fa")

	script:
	"""
	read chrom start end name < <(head -n 1 $bed)
	start=\$(( start > 200 ? start - 200 + 1 : 1 ))
	samtools faidx $ref "\${chrom}:\${start}-\$(( end + 200 ))" > ${target}.fa
	"""
}

process ontMinimap2Index {

	label 'umi_ont'
	cpus params.cpus

	input:
	path ref

	output:
	path "${ref.getBaseName()}.mmi"

	script:
	"""
	minimap2 -x map-ont -k ${params.ontMinimap2Kmer} -t $task.cpus -d ${ref.getBaseName()}.mmi $ref
	"""
}

// The container's medaka ships without models and cannot write its own model cache when run as
// the host user, so the model is fetched once and passed to medaka as a file.
process ontMedakaModel {

	label 'umi_ont'

	input:
	val model

	output:
	path "*_model_pt.tar.gz"

	script:
	"""
	wget -q https://github.com/nanoporetech/medaka/raw/master/medaka/data/${model}_model_pt.tar.gz
	"""
}

process ontFilterReads {

	label 'umi_ont'
	tag "${sampleName}"

	input:
	tuple val(sampleName), path(reads)

	output:
	tuple val(sampleName), path("${sampleName}.filtered.fastq")

	script:
	"""
	catfishq \
	--min-length ${params.ontMinReadLength} \
	--min-qscore ${params.ontMinQscore} \
	$reads > ${sampleName}.filtered.fastq
	"""
}

process ontSubsample {

	label 'umi_ont'
	tag "${sampleName}"
	publishDir "${params.outputLocation}/Analysis/UMI_Nanopore/${sampleName}", mode: 'copy', pattern: "*.tsv"

	input:
	tuple val(sampleName), path(reads)

	output:
	tuple val(sampleName), path("${sampleName}.subsampled.fastq"), emit: reads
	path "${sampleName}.subsampling.tsv"

	script:
	"""
	echo -e "seed\\tsubsampling_readnumber\\n${params.ontSubsampleSeed}\\t${params.ontSubsampleReads}" > ${sampleName}.subsampling.tsv
	seqtk sample -s ${params.ontSubsampleSeed} $reads ${params.ontSubsampleReads} > ${sampleName}.subsampled.fastq
	"""
}

// First-pass alignment of one read chunk, used only to assign reads to amplicons
process ontMapReads {

	label 'umi_ont'
	cpus params.cpus
	tag "${sampleName}"

	input:
	tuple val(sampleName), path(reads)
	path mmi

	output:
	tuple val(sampleName), path("${reads.getBaseName()}.bam"), path("${reads.getBaseName()}.bam.bai")

	script:
	"""
	set -o pipefail
	minimap2 -ax map-ont --MD -t $task.cpus $mmi $reads \
	| samtools sort -@ $task.cpus -o ${reads.getBaseName()}.bam -
	samtools index ${reads.getBaseName()}.bam
	"""
}

// Keep reads that span the amplicon (not too short / long / chimeric), as FASTQ
process ontSplitReads {

	label 'umi_ont'
	tag "${sampleName}:${target}"

	input:
	tuple val(sampleName), path(bam), path(bai), val(target), path(bed)
	path script

	output:
	tuple val(sampleName), val(target), path("*_filtered.fastq"), optional: true, emit: reads
	tuple val(sampleName), val(target), path("*.tsv"), optional: true, emit: stats

	script:
	def secondary = params.ontIncludeSecondary ? "--include_secondary_reads" : ""
	"""
	python $script \
	--min_overlap ${params.ontMinOverlap} \
	--output_format fastq \
	--adapter_length ${params.ontAdapterLength} \
	--output_filename ${sampleName}_${target}_${bam.getBaseName()} \
	$secondary \
	--tsv \
	-o ./ $bed $bam
	"""
}

process ontSummarizeFilterStats {

	label 'umi_ont'
	tag "${sampleName}:${target}"
	publishDir "${umiOntOut(sampleName, target)}/stats", mode: 'copy'

	input:
	tuple val(sampleName), val(target), path(stats)
	path script

	output:
	path "${sampleName}_${target}_filtered_reads.tsv"

	script:
	"""
	python $script ${stats.join(' ')} \
	--sample $sampleName \
	--target $target \
	-o ${sampleName}_${target}_filtered_reads.tsv
	"""
}

// Cut the fwd/rev UMIs out of each read end; the read itself is kept in the header
process ontDetectUmi {

	label 'umi_ont'
	tag "${sampleName}:${target}"

	input:
	tuple val(sampleName), val(target), path(reads)
	val type
	path script

	output:
	tuple val(sampleName), val(target), path("${reads.getBaseName()}_umis.fastq"), emit: umis
	tuple val(sampleName), val(target), path("${reads.getBaseName()}_umis.tsv"), emit: stats

	script:
	def context = params.ontUseContext ? "--use_context" : ""
	"""
	python $script \
	--fwd_umi ${params.ontFwdUmi} \
	--rev_umi ${params.ontRevUmi} \
	--fwd_primer ${params.ontFwdContext} \
	--rev_primer ${params.ontRevContext} \
	--max_error ${params.ontUmiErrors} \
	--adapter_length ${params.ontAdapterLength} \
	--output_format fastq \
	--output_filename ${reads.getBaseName()}_umis \
	$context \
	--tsv \
	-o ./ $reads
	"""
}

process ontSummarizeUmiStats {

	label 'umi_ont'
	tag "${sampleName}:${target}"
	publishDir "${umiOntOut(sampleName, target)}/stats/${type}", mode: 'copy'

	input:
	tuple val(sampleName), val(target), path(stats)
	val type
	path script

	output:
	path "${sampleName}_${target}_umi_summary.tsv"

	script:
	"""
	python $script ${stats.join(' ')} \
	--sample $sampleName \
	--target $target \
	-o ${sampleName}_${target}_umi_summary.tsv
	"""
}

// vsearch clustering of the combined fwd+rev UMI sequences.
// clusters: one FASTA per UMI family; consensus: one centroid per family
process ontCluster {

	label 'umi_ont'
	cpus params.cpus
	tag "${sampleName}:${target}"

	input:
	tuple val(sampleName), val(target), path(umis)
	val identity

	output:
	tuple val(sampleName), val(target), path("cluster*"), optional: true, emit: clusters
	tuple val(sampleName), val(target), path("consensus.fasta"), optional: true, emit: consensus

	script:
	"""
	cat $umis \
	| vsearch \
	--clusterout_id \
	--clusters cluster \
	--consout consensus.fasta \
	--minseqlength ${params.ontUmiMinLength} \
	--maxseqlength ${params.ontUmiMaxLength} \
	--threads $task.cpus \
	--cluster_fast - \
	--clusterout_sort \
	--gapopen 0E/5I \
	--gapext 0E/2I \
	--mismatch -8 \
	--match 6 \
	--iddef 0 \
	--minwordmatches 0 \
	--qmask none \
	--id $identity
	"""
}

// Split admixed UMI clusters and keep the best min..max reads per family (strand balanced)
process ontReformatFilterCluster {

	label 'umi_ont'
	tag "${sampleName}:${target}"
	publishDir "${umiOntOut(sampleName, target)}/stats/raw", mode: 'copy', pattern: "*.tsv"

	input:
	tuple val(sampleName), val(target), path(clusters)
	path script

	output:
	tuple val(sampleName), val(target), path("smolecule*"), optional: true, emit: smolecules
	tuple val(sampleName), val(target), path("*.tsv"), optional: true, emit: stats

	script:
	def balance = params.ontBalanceStrands ? "--balance_strands" : ""
	"""
	python $script \
	--filter_strategy ${params.ontClusterFilter} \
	--min_reads_per_clusters ${params.ontMinReadsPerCluster} \
	--max_reads_per_clusters ${params.ontMaxReadsPerCluster} \
	--max_dist_umi ${params.ontMaxDistUmi} \
	--cluster $clusters \
	--output_format fastq \
	$balance \
	--tsv \
	-o .
	"""
}

process ontClusterStats {

	label 'umi_ont'
	tag "${sampleName}:${target}"
	publishDir "${umiOntOut(sampleName, target)}/stats/raw", mode: 'copy'

	input:
	tuple val(sampleName), val(target), path(stats)
	path script

	output:
	path "${sampleName}_${target}_cluster_stats.tsv"
	path "${sampleName}_${target}_cluster_report.pdf"

	script:
	"""
	python $script \
	--input $stats \
	--output-pdf ${sampleName}_${target}_cluster_report.pdf \
	--output-tsv ${sampleName}_${target}_cluster_stats.tsv \
	--min-reads ${params.ontMinReadsPerCluster} \
	--threshold ${params.ontMinReadsPerCluster - 2}
	"""
}

// Pack families into files of ontClustersPerPolishingFile clusters each (one medaka job per file)
process ontGlueClusters {

	label 'umi_ont'
	tag "${sampleName}:${target}"

	input:
	tuple val(sampleName), val(target), path(smolecules)

	output:
	tuple val(sampleName), val(target), path("${sampleName}_${target}_glued_*.fastq")

	script:
	"""
	file_index=1
	cluster_count=0
	for f in ${smolecules.join(' ')}; do
		if [ "\$cluster_count" -ge ${params.ontClustersPerPolishingFile} ]; then
			file_index=\$((file_index + 1))
			cluster_count=0
		fi
		cluster_count=\$((cluster_count + 1))
		awk -v cnum=\$cluster_count '{ if ((NR - 1) % 4 == 0) print \$0 "_" cnum; else print }' "\$f" \
			>> ${sampleName}_${target}_glued_\${file_index}.fastq
	done
	"""
}

// POA (spoa) + medaka consensus, one read per UMI family
process ontPolishCluster {

	label 'umi_ont'
	cpus 2
	errorStrategy 'retry'
	maxRetries 2
	tag "${sampleName}:${target}"

	input:
	tuple val(key), val(sampleName), val(target), path(glued)
	path model

	output:
	tuple val(key), val(sampleName), val(target), path("${glued.getBaseName()}_consensus.fastq")

	script:
	"""
	medaka smolecule \
	--threads $task.cpus \
	--batch_size ${params.ontClustersPerPolishingFile} \
	--length ${params.ontPolishMinLength} \
	--depth ${params.ontPolishMinDepth} \
	--model $model \
	--method spoa . \
	--qualities \
	$glued 2> ${glued.getBaseName()}_smolecule.log

	mv consensus.fastq ${glued.getBaseName()}_consensus.fastq
	"""
}

// Reference-based alternative to ontPolishCluster: align to the amplicon, medaka inference + sequence.
// Faster and lighter, but the consensus no longer carries the UMIs (so no second UMI round).
process ontPolishClusterRef {

	label 'umi_ont'
	cpus 2
	errorStrategy 'retry'
	maxRetries 2
	tag "${sampleName}:${target}"

	input:
	tuple val(key), val(sampleName), val(target), path(glued), path(targetRef)
	path script
	path model

	output:
	tuple val(key), val(sampleName), val(target), path("${name}_consensus.fastq")

	script:
	name = glued.getBaseName()
	"""
	mini_align -d map-ont -m -t $task.cpus -r $targetRef -i $glued -p $name

	python3 $script \
	--bam ${name}.bam \
	--reference $targetRef \
	--output_reference reference_parsed.fasta \
	--output_bam ${name}_parsed.bam

	medaka inference \
	--check_output \
	--save_features \
	--full_precision \
	--batch_size ${params.ontMedakaBatchSize} \
	--model $model \
	${name}_parsed.bam ${name}.hdf 2> ${name}.medaka.log

	medaka sequence --qualities ${name}.hdf reference_parsed.fasta ${name}_consensus.fastq
	"""
}

// All consensus reads of one sample/amplicon, low-quality bases masked
process ontMergeMaskConsensus {

	label 'umi_ont'
	tag "${sampleName}:${target}"
	publishDir "${umiOntOut(sampleName, target)}/consensus", mode: 'copy'

	input:
	tuple val(sampleName), val(target), path(consensus)

	output:
	tuple val(sampleName), val(target), path("${sampleName}_${target}.consensus.fastq")

	script:
	def hardmask = params.ontMaskingStrategy == "hardmask" ? "--hardmask" : ""
	"""
	cat $consensus > merged.fastq
	vsearch \
	--fastx_mask merged.fastq \
	--fastq_qmax 90 \
	--fastq_qmin ${params.ontMinConsensusQuality} \
	$hardmask \
	--fastqout ${sampleName}_${target}.consensus.fastq
	"""
}

// Second UMI round: families split across clusters in round 1 are merged, one read per molecule
process ontReformatConsensus {

	label 'umi_ont'
	tag "${sampleName}:${target}"
	publishDir "${umiOntOut(sampleName, target)}/consensus", mode: 'copy'

	input:
	tuple val(sampleName), val(target), path(consensusFasta)
	path script

	output:
	tuple val(sampleName), val(target), path("${sampleName}_${target}.final.fastq")

	script:
	"""
	python $script --consensus_fasta $consensusFasta -o .
	mv final.fastq ${sampleName}_${target}.final.fastq
	"""
}

// Consensus reads -> genome. Read names are prefixed with the amplicon so they stay unique
// after the per-sample merge; the read group lets GATK/freebayes see the sample name.
process ontMapConsensus {

	label 'umi_ont'
	cpus params.cpus
	tag "${sampleName}:${target}"

	input:
	tuple val(sampleName), val(target), path(consensus)
	path mmi

	output:
	tuple val(sampleName), path("${sampleName}_${target}.bam")

	script:
	"""
	set -o pipefail
	awk -v t=$target 'NR % 4 == 1 { sub(/^@/, "@" t "_") } { print }' $consensus \
	| minimap2 -ax map-ont --MD -t $task.cpus -R "@RG\\tID:${sampleName}\\tSM:${sampleName}\\tPL:ONT" $mmi - \
	| samtools sort -@ $task.cpus -o ${sampleName}_${target}.bam -
	"""
}

process ontMergeSampleBam {

	label 'umi_ont'
	cpus params.cpus
	tag "${sampleName}"
	publishDir "${params.outputLocation}/Analysis/Alignment/${sampleName}", mode: 'copy'

	input:
	tuple val(sampleName), path(bams)

	output:
	tuple path("${sampleName}.consensus.bam"), path("${sampleName}.consensus.bam.bai")

	script:
	"""
	samtools merge -c -p -f -@ $task.cpus -o ${sampleName}.consensus.bam $bams
	samtools index ${sampleName}.consensus.bam
	"""
}

// LoFreq low-frequency caller on the consensus BAM. Strand-bias filtering is off (as upstream):
// amplicon consensus reads are not strand balanced. The VCF gets a GT:AD:DP:AF sample column
// so the tertiary loader reads it like a Mutect2 VCF.
process ontLofreq {

	label 'umi_ont'
	cpus params.cpus
	tag "${sampleName}"
	publishDir "${params.outputLocation}/Analysis/Variants/${sampleName}/LoFreq", mode: 'copy'

	input:
	tuple path(bam), path(bai)
	tuple path(ref), path(refIndex)
	path bed

	output:
	path "${sampleName}.vcf"

	script:
	sampleName = bam.getSimpleName()
	"""
	lofreq indelqual --dindel -f $ref -o indelqual.bam $bam
	samtools index indelqual.bam

	lofreq call-parallel --pp-threads $task.cpus \
	-f $ref \
	-l $bed \
	--call-indels \
	--min-cov ${params.ontMinDepth} \
	--no-default-filter \
	-o raw.vcf \
	indelqual.bam

	lofreq filter --no-defaults \
	--cov-min ${params.ontMinDepth} \
	--af-min ${params.ontMinVaf} \
	-i raw.vcf -o filtered.vcf

	python3 $projectDir/bin/lofreq_add_sample.py filtered.vcf $sampleName > ${sampleName}.vcf
	"""
}

process ontVcfStats {

	label 'umi_ont'
	tag "${sampleName}"
	publishDir "${params.outputLocation}/Analysis/Variants/${sampleName}/LoFreq", mode: 'copy'

	input:
	path vcf

	output:
	path "${sampleName}.bcftools_stats.txt"

	script:
	sampleName = vcf.getSimpleName()
	"""
	bcftools stats -s - $vcf > ${sampleName}.bcftools_stats.txt
	"""
}
