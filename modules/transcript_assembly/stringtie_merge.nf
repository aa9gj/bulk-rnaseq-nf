process STRINGTIE_MERGE {
    input:
    path(gtf_list)

    output:
    path("merged.gtf")

    script:
    """
    stringtie --merge -G ${params.gtf_annotation} -o merged.gtf ${gtf_list}
    """
}
