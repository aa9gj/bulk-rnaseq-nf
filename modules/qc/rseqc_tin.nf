/*
 * RSEQC_TIN
 *
 * Calculate the Transcript Integrity Number (TIN) for each sample.
 * TIN is a reliable metric for RNA quality that is directly computed
 * from the BAM file. It is analogous to RIN (RNA Integrity Number)
 * but can be calculated post-sequencing.
 *
 * THIS IS THE MOST IMPORTANT RSEQC MODULE FOR DIAGNOSING PCA OUTLIERS.
 *
 * Interpretation:
 *   - TIN > 70: High quality RNA
 *   - TIN 50-70: Moderate degradation
 *   - TIN < 50: Severe degradation (likely PCA outlier)
 *
 * Samples with low TIN scores should be flagged for removal or
 * included as a covariate in differential expression analysis.
 *
 * Input:
 *   - sample_id: Sample identifier
 *   - bam: Sorted BAM file
 *   - bai: BAM index file
 *   - bed: BED12 gene model file
 *
 * Output:
 *   - tin_scores: Per-transcript TIN scores
 *   - summary: Summary TIN statistics (median TIN across all transcripts)
 *
 * Tools: RSeQC tin.py (https://rseqc.sourceforge.net/)
 */

process RSEQC_TIN {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/rseqc/tin", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(bed)

    output:
    path("*.tin.xls"), emit: tin_scores
    path("*.summary.txt"), emit: summary

    script:
    """
    tin.py \\
        -i ${bam} \\
        -r ${bed}
    """
}
