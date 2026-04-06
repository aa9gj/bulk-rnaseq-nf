/*
 * GTF_TO_BED
 *
 * Convert GTF annotation file to BED12 format.
 * RSeQC tools require a BED12 gene model file. This process converts
 * the GTF annotation automatically so users don't need to provide
 * a separate BED file.
 *
 * Input:
 *   - gtf: GTF annotation file
 *
 * Output:
 *   - bed: BED12 gene model file for RSeQC
 *
 * Tools: gtf_to_bed12.py (bundled in bin/)
 */

process GTF_TO_BED {
    tag "gtf_to_bed"
    label 'process_low'

    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path(gtf)

    output:
    path("gene_model.bed"), emit: bed

    script:
    """
    gtf_to_bed12.py -i ${gtf} -o gene_model.bed
    """
}
