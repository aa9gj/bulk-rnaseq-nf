/*
 * TRIM_GALORE
 *
 * Perform adapter trimming and quality control using Trim Galore.
 * Trim Galore is a wrapper around Cutadapt and FastQC that automates
 * adapter detection and quality trimming for paired-end reads.
 *
 * Input:
 *   - sample_id: Sample identifier
 *   - R1: Forward reads (R1) FASTQ file
 *   - R2: Reverse reads (R2) FASTQ file
 *
 * Output:
 *   - trimmed: Tuple of (sample_id, trimmed R1, trimmed R2)
 *   - qc_reports: FastQC reports for MultiQC
 *
 * Tools: Trim Galore (https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
 */

process TRIM_GALORE {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: "*_val_{1,2}.fq.gz"
    publishDir "${params.outdir}/fastqc", mode: 'copy', pattern: "*_fastqc.{html,zip}"
    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: "*_trimming_report.txt"

    input:
    tuple val(sample_id), path(R1), path(R2)

    output:
    tuple val(sample_id), path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed
    path("*_fastqc.{html,zip}"), emit: qc_reports
    path("*_trimming_report.txt"), emit: logs

    script:
    def cores = task.cpus ?: params.threads ?: 4
    """
    # Run Trim Galore with FastQC on paired-end reads
    trim_galore \\
        --paired \\
        --fastqc \\
        --cores ${cores} \\
        --gzip \\
        ${R1} ${R2}
    """
}
