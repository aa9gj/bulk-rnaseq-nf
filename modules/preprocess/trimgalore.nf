process TRIM_GALORE {
    input:
    tuple val(sample_id), path(R1), path(R2)

    output:
    tuple val(sample_id), path("*.fq.gz")

    script:
    """
    trim_galore --paired --fastqc -o trimmed/ ${R1} ${R2}
    mv trimmed/*.fq.gz .
    """
}
