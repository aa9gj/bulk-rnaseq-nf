/*
 * RSEQC_INFER_EXPERIMENT
 *
 * Infer the strandedness of the RNA-seq library. This is critical
 * for accurate quantification. Incorrect strandedness settings are
 * a common cause of samples appearing as outliers in PCA plots,
 * as reads get assigned to the wrong genes.
 *
 * Output interpretation:
 *   - "1++,1--,2+-,2-+": Sense (stranded, e.g. dUTP second-strand)
 *   - "1+-,1-+,2++,2--": Antisense (stranded, first-strand)
 *   - ~50/50 split:       Unstranded library
 *
 * Input:
 *   - sample_id: Sample identifier
 *   - bam: Sorted BAM file
 *   - bai: BAM index file
 *   - bed: BED12 gene model file
 *
 * Output:
 *   - strandedness: Strandedness inference results
 *
 * Tools: RSeQC infer_experiment.py (https://rseqc.sourceforge.net/)
 */

process RSEQC_INFER_EXPERIMENT {
    tag "${sample_id}"
    label 'process_low'

    publishDir "${params.outdir}/rseqc/infer_experiment", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(bed)

    output:
    path("${sample_id}.infer_experiment.txt"), emit: strandedness

    script:
    """
    infer_experiment.py \\
        -i ${bam} \\
        -r ${bed} \\
        -s 200000 \\
        > ${sample_id}.infer_experiment.txt
    """
}
