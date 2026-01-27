#!/usr/bin/env nextflow

/*
 * bulk-rnaseq-nf: A Nextflow pipeline for bulk RNA-seq analysis
 *
 * Pipeline steps:
 *   1. Concatenate multi-lane FASTQ files
 *   2. Quality control and adapter trimming (Trim Galore)
 *   3. Alignment to reference genome (HISAT2)
 *   4. BAM sorting and indexing (SAMtools)
 *   5. Transcript assembly (StringTie - first pass)
 *   6. Merge transcript assemblies (StringTie merge)
 *   7. Quantification (StringTie - second pass)
 *   8. Prepare count matrices (prepDE.py)
 *   9. MultiQC report generation
 */

nextflow.enable.dsl = 2

// Import all process modules
include { CONCAT_FASTQ } from '../modules/preprocess/concat_fastq.nf'
include { TRIM_GALORE } from '../modules/preprocess/trim_galore.nf'
include { HISAT2_INDEX } from '../modules/align/hisat2_index.nf'
include { HISAT2_ALIGN } from '../modules/align/hisat2_align.nf'
include { SAMTOOLS_SORT } from '../modules/align/samtools_sort.nf'
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
    CONCAT_FASTQ(samples_ch)

    /*
     * STEP 3: Quality control and adapter trimming
     * Uses Trim Galore (wrapper for Cutadapt + FastQC)
     */
    TRIM_GALORE(CONCAT_FASTQ.out.reads)

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
     * STEP 7: First-pass transcript assembly with StringTie
     * Uses reference annotation to guide assembly
     */
    gtf_annotation_ch = Channel.value(file(params.gtf_annotation))
    STRINGTIE_FIRST(
        SAMTOOLS_SORT.out.sorted_bam.combine(gtf_annotation_ch)
    )

    /*
     * STEP 8: Merge all transcript assemblies
     * Creates a unified set of transcripts across all samples
     */
    STRINGTIE_MERGE(
        STRINGTIE_FIRST.out.gtf.map { sample_id, gtf -> gtf }.collect()
    )

    /*
     * STEP 9: Second-pass quantification with StringTie
     * Re-estimates abundances using the merged annotation
     */
    STRINGTIE_SECOND(
        SAMTOOLS_SORT.out.sorted_bam.combine(STRINGTIE_MERGE.out.merged_gtf)
    )

    /*
     * STEP 10: Generate count matrices for differential expression
     */
    PREPDE(
        STRINGTIE_SECOND.out.gtf.map { sample_id, gtf -> gtf }.collect()
    )

    /*
     * STEP 11: Generate MultiQC report
     * Collect QC outputs from various steps
     */
    qc_files_ch = TRIM_GALORE.out.qc_reports
        .mix(TRIM_GALORE.out.logs)
        .mix(HISAT2_ALIGN.out.logs)
        .mix(SAMTOOLS_SORT.out.stats)
        .mix(STRINGTIE_FIRST.out.stats)
        .collect()

    MULTIQC(qc_files_ch)
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
