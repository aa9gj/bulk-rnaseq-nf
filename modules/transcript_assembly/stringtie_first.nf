process STRINGTIE_FIRST {
    input:
    tuple val(sample_id), path(bam), path(gtf_annotation)

    output:
    tuple val(sample_id), path("${sample_id}.gtf")

    script:
    """
    stringtie ${bam} -G ${gtf_annotation} -o ${sample_id}.gtf -p 8
    """
}
