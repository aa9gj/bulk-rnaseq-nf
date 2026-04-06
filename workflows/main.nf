#!/usr/bin/env nextflow

/*
 * bulk-rnaseq-nf: A Nextflow pipeline for bulk RNA-seq analysis
 *
 * Pipeline steps:
 *   1. Concatenate multi-lane FASTQ files
 *   2. Quality control and adapter trimming (Trim Galore)
 *   3. Alignment to reference genome (HISAT2)
 *   4. BAM sorting and indexing (SAMtools)
 *   5. RSeQC quality control suite
 *   6. Transcript assembly (StringTie - first pass)
 *   7. Merge transcript assemblies (StringTie merge)
 *   8. Quantification (StringTie - second pass)
 *   9. Prepare count matrices (prepDE.py)
 *  10. MultiQC report generation
 */

nextflow.enable.dsl = 2

// Import all process modules
include { CONCAT_FASTQ } from '../modules/preprocess/concat_fastq.nf'
include { TRIM_GALORE } from '../modules/preprocess/trim_galore.nf'
include { HISAT2_INDEX } from '../modules/align/hisat2_index.nf'
include { HISAT2_ALIGN } from '../modules/align/hisat2_align.nf'
include { SAMTOOLS_SORT } from '../modules/align/samtools_sort.nf'
include { GTF_TO_BED } from '../modules/qc/gtf_to_bed.nf'
include { RSEQC_BAM_STAT } from '../modules/qc/rseqc_bam_stat.nf'
include { RSEQC_INFER_EXPERIMENT } from '../modules/qc/rseqc_infer_experiment.nf'
include { RSEQC_READ_DISTRIBUTION } from '../modules/qc/rseqc_read_distribution.nf'
include { RSEQC_INNER_DISTANCE } from '../modules/qc/rseqc_inner_distance.nf'
include { RSEQC_GENE_BODY_COVERAGE } from '../modules/qc/rseqc_gene_body_coverage.nf'
include { RSEQC_TIN } from '../modules/qc/rseqc_tin.nf'
include { STRINGTIE_FIRST } from '../modules/transcript_assembly/stringtie_first.nf'
include { STRINGTIE_MERGE } from '../modules/transcript_assembly/stringtie_merge.nf'
include { STRINGTIE_SECOND } from '../modules/transcript_assembly/stringtie_second.nf'
include { PREPDE } from '../modules/transcript_assembly/prepde.nf'
include { MULTIQC } from '../modules/qc/multiqc.nf'

/*
 * Validate required parameters
 */
def validateParams() {
    if (!params.samples) {
        error "ERROR: 'samples' parameter is required. Please provide sample information in params.yaml"
    }
    if (!params.gtf_annotation) {
        error "ERROR: 'gtf_annotation' parameter is required"
    }
    if (!params.hisat2_index && !params.genome_fasta) {
        error "ERROR: Either 'hisat2_index' or 'genome_fasta' must be provided"
    }
}

/*
 * Main workflow
 */
workflow {
    // Validate parameters before running
    validateParams()

    /*
     * STEP 1: Create sample channel from params.yaml
     * Parse the samples map and create tuples of (sample_id, R1_files, R2_files)
     */
    samples_ch = Channel
        .fromList(params.samples.collect { sample_id, reads ->
            tuple(
                sample_id,
                reads.R1.tokenize(',').collect { it.trim() },
                reads.R2.tokenize(',').collect { it.trim() }
            )
        })

    /*
     * STEP 2: Concatenate multi-lane FASTQ files
     * For samples sequenced across multiple lanes
     */
   //  CONCAT_FASTQ(samples_ch)

    /*
     * STEP 3: Quality control and adapter trimming
     * Uses Trim Galore (wrapper for Cutadapt + FastQC)
     */
    TRIM_GALORE(samples_ch)

    /*
     * STEP 4: Prepare HISAT2 index
     * Either use provided index or generate from genome FASTA
     */
    if (params.hisat2_index) {
        hisat2_index_ch = Channel.value(file(params.hisat2_index))
    } else {
        HISAT2_INDEX(file(params.genome_fasta))
        hisat2_index_ch = HISAT2_INDEX.out.index
    }

    /*
     * STEP 5: Align reads to reference genome using HISAT2
     * Combine trimmed reads with the index
     */
    HISAT2_ALIGN(
        TRIM_GALORE.out.trimmed.combine(hisat2_index_ch)
    )

    /*
     * STEP 6: Sort and index BAM files
     */
    SAMTOOLS_SORT(HISAT2_ALIGN.out.bam)

    /*
     * STEP 7: RSeQC Quality Control Suite
     *
     * Convert GTF to BED12 format (required by RSeQC), then run
     * all RSeQC modules in parallel on sorted BAMs.
     *
     * These QC checks help diagnose PCA outliers by identifying:
     *   - Strandedness mismatches (infer_experiment)
     *   - RNA degradation (TIN scores, gene body coverage)
     *   - Abnormal read distributions
     *   - Insert size anomalies
     */
    if (!params.skip_rseqc) {
        // Convert GTF annotation to BED12 for RSeQC
        gtf_annotation_file = file(params.gtf_annotation)
        GTF_TO_BED(gtf_annotation_file)

        // Prepare input: sorted BAM + BED gene model
        sorted_bam_with_bed = SAMTOOLS_SORT.out.sorted_bam.combine(GTF_TO_BED.out.bed)

        // Run all RSeQC modules in parallel
        // These have no dependencies on each other and can execute concurrently
        RSEQC_BAM_STAT(SAMTOOLS_SORT.out.sorted_bam)
        RSEQC_INFER_EXPERIMENT(sorted_bam_with_bed)
        RSEQC_READ_DISTRIBUTION(sorted_bam_with_bed)
        RSEQC_INNER_DISTANCE(sorted_bam_with_bed)
        RSEQC_GENE_BODY_COVERAGE(sorted_bam_with_bed)
        RSEQC_TIN(sorted_bam_with_bed)
    }

    /*
     * STEP 8: First-pass transcript assembly with StringTie
     * Uses reference annotation to guide assembly
     */
    gtf_annotation_ch = Channel.value(file(params.gtf_annotation))
    STRINGTIE_FIRST(
        SAMTOOLS_SORT.out.sorted_bam.combine(gtf_annotation_ch)
    )

    /*
     * STEP 9: Merge all transcript assemblies
     * Creates a unified set of transcripts across all samples
     */
    STRINGTIE_MERGE(
        STRINGTIE_FIRST.out.gtf.map { sample_id, gtf -> gtf }.collect()
    )

    /*
     * STEP 10: Second-pass quantification with StringTie
     * Re-estimates abundances using the merged annotation
     */
    STRINGTIE_SECOND(
        SAMTOOLS_SORT.out.sorted_bam.combine(STRINGTIE_MERGE.out.merged_gtf)
    )

    /*
     * STEP 11: Generate count matrices for differential expression
     */
    PREPDE(
        STRINGTIE_SECOND.out.gtf.map { sample_id, gtf -> gtf }.collect()
    )

    /*
     * STEP 12: Generate MultiQC report
     * Collect QC outputs from Trim Galore, HISAT2, SAMtools, StringTie, and RSeQC
     */
    qc_files_ch = TRIM_GALORE.out.qc_reports
        .mix(TRIM_GALORE.out.logs)
        .mix(HISAT2_ALIGN.out.logs)
        .mix(SAMTOOLS_SORT.out.stats)
        .mix(STRINGTIE_FIRST.out.stats)

    // Add RSeQC outputs to MultiQC if RSeQC was run
    if (!params.skip_rseqc) {
        qc_files_ch = qc_files_ch
            .mix(RSEQC_BAM_STAT.out.stats)
            .mix(RSEQC_INFER_EXPERIMENT.out.strandedness)
            .mix(RSEQC_READ_DISTRIBUTION.out.distribution)
            .mix(RSEQC_INNER_DISTANCE.out.distance)
            .mix(RSEQC_GENE_BODY_COVERAGE.out.coverage)
            .mix(RSEQC_TIN.out.tin_scores)
            .mix(RSEQC_TIN.out.summary)
    }

    MULTIQC(qc_files_ch.collect())
}

/*
 * Workflow completion handler
 */
workflow.onComplete {
    log.info """
    ============================================
    Pipeline completed!
    ============================================
    Status:     ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration:   ${workflow.duration}
    Results:    ${params.outdir ?: 'results'}
    ============================================
    """.stripIndent()
}
