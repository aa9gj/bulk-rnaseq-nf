/*
 * FEATURECOUNTS
 *
 * Count reads mapping to genomic features using featureCounts (Subread).
 * featureCounts is fast and memory-efficient, making it ideal for
 * large-scale RNA-seq experiments.
 *
 * Input:
 *   - bam_files: Collection of sorted BAM files
 *   - gtf: GTF annotation file
 *
 * Output:
 *   - counts: Gene count matrix
 *   - summary: Assignment summary statistics
 *
 * Tools: featureCounts (http://subread.sourceforge.net/)
 */

process FEATURECOUNTS {
    tag "featureCounts"
    label 'process_medium'

    publishDir "${params.outdir}/featurecounts", mode: 'copy'

    input:
    path(bam_files)
    path(gtf)

    output:
    path("featurecounts.txt"), emit: counts
    path("featurecounts.txt.summary"), emit: summary

    script:
    def threads = task.cpus ?: params.threads ?: 8
    // Strandedness: 0 = unstranded, 1 = stranded, 2 = reversely stranded
    def strandedness = params.strandedness ?: 0
    """
    featureCounts \\
        -T ${threads} \\
        -p \\
        -t exon \\
        -g gene_id \\
        -s ${strandedness} \\
        -a ${gtf} \\
        -o featurecounts.txt \\
        ${bam_files}
    """
}
