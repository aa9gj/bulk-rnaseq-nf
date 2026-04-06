/*
 * RSEQC_READ_DISTRIBUTION
 *
 * Calculate how mapped reads are distributed across genomic features:
 * CDS exons, 5'UTR, 3'UTR, introns, and intergenic regions.
 * Samples with abnormal distributions (e.g., high intergenic or
 * intronic fractions) often appear as outliers in PCA.
 *
 * Input:
 *   - sample_id: Sample identifier
 *   - bam: Sorted BAM file
 *   - bai: BAM index file
 *   - bed: BED12 gene model file
 *
 * Output:
 *   - distribution: Read distribution statistics
 *
 * Tools: RSeQC read_distribution.py (https://rseqc.sourceforge.net/)
 */

process RSEQC_READ_DISTRIBUTION {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/rseqc/read_distribution", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(bed)

    output:
    path("${sample_id}.read_distribution.txt"), emit: distribution

    script:
    """
    read_distribution.py \\
        -i ${bam} \\
        -r ${bed} \\
        > ${sample_id}.read_distribution.txt
    """
}
