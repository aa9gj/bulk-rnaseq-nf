/*
 * STRINGTIE_FIRST
 *
 * First-pass transcript assembly with StringTie.
 * Assembles transcripts for each sample using the reference annotation
 * as a guide. This is part of the StringTie two-pass approach for
 * more accurate transcript quantification.
 *
 * Input:
 *   - sample_id: Sample identifier
 *   - bam: Sorted BAM file from alignment
 *   - bai: BAM index file
 *   - gtf_annotation: Reference GTF annotation file
 *
 * Output:
 *   - gtf: Assembled transcripts in GTF format
 *   - stats: Assembly statistics for QC
 *
 * Tools: StringTie (https://ccb.jhu.edu/software/stringtie/)
 */

process STRINGTIE_FIRST {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/stringtie/first_pass", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(gtf_annotation)

    output:
    tuple val(sample_id), path("${sample_id}.gtf"), emit: gtf
    path("${sample_id}.stats.txt"), emit: stats

    script:
    def threads = task.cpus ?: params.threads ?: 8
    """
    # Run StringTie first-pass assembly
    stringtie \\
        ${bam} \\
        -G ${gtf_annotation} \\
        -o ${sample_id}.gtf \\
        -p ${threads} \\
        -A ${sample_id}.stats.txt \\
        -l ${sample_id}
    """
}
