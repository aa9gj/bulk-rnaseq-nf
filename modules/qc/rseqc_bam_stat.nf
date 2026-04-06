/*
 * RSEQC_BAM_STAT
 *
 * Summarize mapping quality, duplication rate, and other basic
 * alignment statistics. Useful for identifying low-quality samples
 * that may appear as PCA outliers.
 *
 * Input:
 *   - sample_id: Sample identifier
 *   - bam: Sorted BAM file
 *   - bai: BAM index file
 *
 * Output:
 *   - stats: BAM statistics text file
 *
 * Tools: RSeQC bam_stat.py (https://rseqc.sourceforge.net/)
 */

process RSEQC_BAM_STAT {
    tag "${sample_id}"
    label 'process_low'

    publishDir "${params.outdir}/rseqc/bam_stat", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path("${sample_id}.bam_stat.txt"), emit: stats

    script:
    """
    bam_stat.py \\
        -i ${bam} \\
        > ${sample_id}.bam_stat.txt
    """
}
