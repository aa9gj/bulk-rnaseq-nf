process HISAT2_ALIGN {
    input:
    tuple val(sample_id), path(R1), path(R2), path(index)

    output:
    tuple val(sample_id), path("*.bam")

    script:
    """
    hisat2 -x ${index} -1 ${R1} -2 ${R2} | samtools view -bS - > ${sample_id}.bam
    """
}
