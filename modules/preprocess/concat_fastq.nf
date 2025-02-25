process CONCAT_FASTQ {
    input:
    tuple val(sample_id), path(reads_R1), path(reads_R2)

    output:
    tuple val(sample_id), path("*.fastq.gz")

    script:
    """
    cat ${reads_R1} > ${sample_id}_R1.fastq.gz
    cat ${reads_R2} > ${sample_id}_R2.fastq.gz
    """
}

