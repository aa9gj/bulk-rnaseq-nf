/*
 * HISAT2_ALIGN
 *
 * Align paired-end RNA-seq reads to a reference genome using HISAT2.
 * HISAT2 is a fast and sensitive alignment program for mapping
 * next-generation sequencing reads to a population of human genomes.
 *
 * Input:
 *   - sample_id: Sample identifier
 *   - R1: Trimmed forward reads (R1)
 *   - R2: Trimmed reverse reads (R2)
 *   - index: Path to HISAT2 index
 *
 * Output:
 *   - bam: Aligned reads in BAM format (unsorted)
 *   - logs: Alignment summary statistics
 *
 * Tools: HISAT2 (https://ccb.jhu.edu/software/hisat2/)
 *        SAMtools (http://www.htslib.org/)
 */

process HISAT2_ALIGN {
    tag "${sample_id}"
    label 'process_high'

    publishDir "${params.outdir}/hisat2", mode: 'copy', pattern: "*.log"

    input:
    tuple val(sample_id), path(R1), path(R2), path(index)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    path("${sample_id}.hisat2.log"), emit: logs

    script:
    def threads = task.cpus ?: params.threads ?: 8
    // Determine index prefix (handle both directory and file prefix formats)
    def index_prefix = index.isDirectory() ? "${index}/genome" : index
    """
    # Align reads with HISAT2 and convert to BAM
    hisat2 \\
        -p ${threads} \\
        --dta \\
        --rna-strandness RF \\
        --summary-file ${sample_id}.hisat2.log \\
        -x ${index_prefix} \\
        -1 ${R1} \\
        -2 ${R2} \\
        | samtools view -@ ${threads} -bS - > ${sample_id}.bam
    """
}
