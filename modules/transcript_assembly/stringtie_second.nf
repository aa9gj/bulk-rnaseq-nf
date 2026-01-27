/*
 * STRINGTIE_SECOND
 *
 * Second-pass quantification with StringTie.
 * Re-estimates transcript abundances using the merged annotation
 * from all samples. The -e flag restricts quantification to the
 * transcripts in the merged GTF, and -B generates Ballgown tables.
 *
 * Input:
 *   - sample_id: Sample identifier
 *   - bam: Sorted BAM file
 *   - bai: BAM index file
 *   - merged_gtf: Merged transcript annotation from STRINGTIE_MERGE
 *
 * Output:
 *   - gtf: Quantified transcripts with abundance estimates
 *   - ballgown: Ballgown table files for downstream analysis
 *
 * Tools: StringTie (https://ccb.jhu.edu/software/stringtie/)
 */

process STRINGTIE_SECOND {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/stringtie/quantification", mode: 'copy'
    publishDir "${params.outdir}/ballgown/${sample_id}", mode: 'copy', pattern: "*.ctab"

    input:
    tuple val(sample_id), path(bam), path(bai), path(merged_gtf)

    output:
    tuple val(sample_id), path("${sample_id}_quant.gtf"), emit: gtf
    path("*.ctab"), emit: ballgown

    script:
    def threads = task.cpus ?: params.threads ?: 8
    """
    # Run StringTie second-pass quantification
    # -e: Only estimate abundance of transcripts in the reference
    # -B: Generate Ballgown table files for downstream analysis
    stringtie \\
        ${bam} \\
        -G ${merged_gtf} \\
        -o ${sample_id}_quant.gtf \\
        -p ${threads} \\
        -e \\
        -B
    """
}
