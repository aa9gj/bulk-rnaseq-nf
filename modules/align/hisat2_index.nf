process HISAT2_INDEX {
    input:
    path(genome_fasta)

    output:
    path("hisat2_index")

    when:
    params.hisat2_index == ''

    script:
    """
    hisat2-build ${genome_fasta} hisat2_index
    """
}
