/*
 * HISAT2_INDEX
 *
 * Build a HISAT2 index from a reference genome FASTA file.
 * This step is only required if a pre-built index is not available.
 *
 * Input:
 *   - genome_fasta: Reference genome in FASTA format
 *
 * Output:
 *   - index: Directory containing HISAT2 index files
 *
 * Note: Building an index is computationally intensive and should
 *       only be done once per reference genome. Pre-built indices
 *       for common genomes are available from the HISAT2 website.
 *
 * Tools: HISAT2 (https://ccb.jhu.edu/software/hisat2/)
 */

process HISAT2_INDEX {
    tag "${genome_fasta.simpleName}"
    label 'process_high'

    publishDir "${params.outdir}/index", mode: 'copy'

    input:
    path(genome_fasta)

    output:
    path("hisat2_index"), emit: index

    script:
    def threads = task.cpus ?: params.threads ?: 8
    """
    # Create output directory
    mkdir -p hisat2_index

    # Build HISAT2 index
    hisat2-build \\
        -p ${threads} \\
        ${genome_fasta} \\
        hisat2_index/genome
    """
}
