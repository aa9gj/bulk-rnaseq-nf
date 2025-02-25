process STRINGTIE_SECOND {
    input:
    tuple val(sample_id), path(bam), path(merged_gtf)

    output:
    tuple val(sample_id), path("${sample_id}_quant.gtf")

    script:
    """
    stringtie ${bam} -G ${merged_gtf} -o ${sample_id}_quant.gtf -p 8 -e -B
    """
}
