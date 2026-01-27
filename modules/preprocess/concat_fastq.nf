/*
 * CONCAT_FASTQ
 *
 * Concatenate FASTQ files from samples sequenced across multiple lanes.
 * This is common when samples are split across multiple sequencing lanes
 * to increase coverage.
 *
 * Input:
 *   - sample_id: Sample identifier
 *   - reads_R1: List of R1 (forward) FASTQ files to concatenate
 *   - reads_R2: List of R2 (reverse) FASTQ files to concatenate
 *
 * Output:
 *   - Tuple of (sample_id, concatenated R1, concatenated R2)
 *
 * Note: Uses 'cat' for gzipped files which preserves the gzip format
 *       when concatenating multiple .gz files.
 */

process CONCAT_FASTQ {
    tag "${sample_id}"
    label 'process_low'

    input:
    tuple val(sample_id), val(reads_R1), val(reads_R2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq.gz"), path("${sample_id}_R2.fastq.gz"), emit: reads

    script:
    // Convert file paths to actual file objects
    def r1_files = reads_R1.collect { file(it) }.join(' ')
    def r2_files = reads_R2.collect { file(it) }.join(' ')
    """
    # Concatenate R1 (forward) reads
    cat ${r1_files} > ${sample_id}_R1.fastq.gz

    # Concatenate R2 (reverse) reads
    cat ${r2_files} > ${sample_id}_R2.fastq.gz
    """
}
