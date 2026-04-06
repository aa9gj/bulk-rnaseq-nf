/*
 * RSEQC_INNER_DISTANCE
 *
 * Calculate the inner distance (insert size) between paired-end reads.
 * Abnormal insert size distributions can indicate library preparation
 * problems and may contribute to PCA outlier behavior.
 *
 * Input:
 *   - sample_id: Sample identifier
 *   - bam: Sorted BAM file
 *   - bai: BAM index file
 *   - bed: BED12 gene model file
 *
 * Output:
 *   - distance: Inner distance statistics
 *   - pdf: Inner distance distribution plot
 *   - freq: Frequency table
 *
 * Tools: RSeQC inner_distance.py (https://rseqc.sourceforge.net/)
 */

process RSEQC_INNER_DISTANCE {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/rseqc/inner_distance", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(bed)

    output:
    path("${sample_id}.inner_distance*.txt"), emit: distance
    path("${sample_id}.inner_distance*.pdf"), emit: pdf, optional: true
    path("${sample_id}.inner_distance*.{r,R}"), emit: rscript, optional: true

    script:
    """
    inner_distance.py \\
        -i ${bam} \\
        -r ${bed} \\
        -o ${sample_id}.inner_distance
    """
}
