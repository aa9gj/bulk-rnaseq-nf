/*
 * SAMTOOLS_SORT
 *
 * Sort and index BAM files using SAMtools.
 * Sorting is required for downstream analysis tools like StringTie.
 *
 * Input:
 *   - sample_id: Sample identifier
 *   - bam: Unsorted BAM file from HISAT2 alignment
 *
 * Output:
 *   - sorted_bam: Coordinate-sorted BAM file
 *   - bai: BAM index file (.bai)
 *
 * Tools: SAMtools (http://www.htslib.org/)
 */

process SAMTOOLS_SORT {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/aligned", mode: 'copy', pattern: "*.{bam,bai}"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: sorted_bam
    path("${sample_id}.flagstat.txt"), emit: stats

    script:
    def threads = task.cpus ?: params.threads ?: 4
    """
    # Sort BAM file by coordinate
    samtools sort \\
        -@ ${threads} \\
        -o ${sample_id}.sorted.bam \\
        ${bam}

    # Index the sorted BAM file
    samtools index \\
        -@ ${threads} \\
        ${sample_id}.sorted.bam

    # Generate alignment statistics
    samtools flagstat \\
        ${sample_id}.sorted.bam > ${sample_id}.flagstat.txt
    """
}
