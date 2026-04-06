/*
 * RSEQC_GENE_BODY_COVERAGE
 *
 * Check the uniformity of read coverage across gene bodies.
 * A skew toward the 3' end indicates RNA degradation, which is
 * a major cause of PCA outliers. Degraded samples will have
 * systematically different expression profiles.
 *
 * Interpretation:
 *   - Uniform coverage: Good quality RNA
 *   - 3' bias: RNA degradation (common PCA outlier cause)
 *   - 5' bias: Possible ribosomal RNA contamination
 *
 * Input:
 *   - sample_id: Sample identifier
 *   - bam: Sorted BAM file
 *   - bai: BAM index file
 *   - bed: BED12 gene model file
 *
 * Output:
 *   - coverage: Coverage statistics text file
 *   - pdf: Gene body coverage plot
 *
 * Tools: RSeQC geneBody_coverage.py (https://rseqc.sourceforge.net/)
 */

process RSEQC_GENE_BODY_COVERAGE {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/rseqc/gene_body_coverage", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(bed)

    output:
    path("${sample_id}.geneBodyCoverage.txt"), emit: coverage
    path("${sample_id}.geneBodyCoverage.curves.pdf"), emit: pdf, optional: true
    path("${sample_id}.geneBodyCoverage.{r,R}"), emit: rscript, optional: true

    script:
    """
    geneBody_coverage.py \\
        -i ${bam} \\
        -r ${bed} \\
        -o ${sample_id}
    """
}
