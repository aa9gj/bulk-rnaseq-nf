process MULTIQC {
    input:
    path(qc_files)

    output:
    path("multiqc_report.html")

    publishDir "results/multiqc", mode: 'copy'

    script:
    """
    multiqc ${qc_files} -o .
    """
}
