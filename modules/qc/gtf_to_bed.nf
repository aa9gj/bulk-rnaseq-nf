/*
 * GTF_TO_BED
 *
 * Convert GTF annotation file to BED12 format.
 * RSeQC tools require a BED12 gene model file. This process converts
 * the GTF annotation using the standard UCSC utilities:
 *   gtfToGenePred -> genePredToBed
 *
 * Input:
 *   - gtf: GTF annotation file
 *
 * Output:
 *   - bed: BED12 gene model file for RSeQC
 *
 * Tools: UCSC gtfToGenePred and genePredToBed
 *        (conda install -c bioconda ucsc-gtftogenepred ucsc-genepredtobed)
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
    gtfToGenePred ${gtf} gene_model.genePred
    genePredToBed gene_model.genePred gene_model.bed
    """
}
