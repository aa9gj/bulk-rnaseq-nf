/*
 * MULTIQC
 *
 * Aggregate quality control reports from multiple tools into a
 * single interactive HTML report using MultiQC.
 *
 * Input:
 *   - qc_files: Collection of QC files from various pipeline steps
 *               (FastQC, HISAT2 logs, StringTie stats, etc.)
 *
 * Output:
 *   - report: Interactive HTML report
 *   - data: MultiQC data directory for programmatic access
 *
 * Tools: MultiQC (https://multiqc.info/)
 */

process MULTIQC {
    tag "multiqc"
    label 'process_low'

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path(qc_files)

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data"), emit: data

    script:
    """
    # Run MultiQC on all collected QC files
    multiqc \\
        --force \\
        --title "bulk-rnaseq-nf QC Report" \\
        --filename multiqc_report.html \\
        .
    """
}
