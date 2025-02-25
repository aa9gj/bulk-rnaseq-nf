process PREPDE {
    input:
    path(gtf_files)

    output:
    path("gene_count_matrix.csv"), path("transcript_count_matrix.csv")

    script:
    """
    python3 PrepDE.py -i ${gtf_files}
    """
}
