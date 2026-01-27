/*
 * STRINGTIE_MERGE
 *
 * Merge transcript assemblies from all samples into a unified annotation.
 * This creates a non-redundant set of transcripts observed across all
 * samples, which is then used for consistent quantification.
 *
 * Input:
 *   - gtf_list: Collection of GTF files from first-pass assembly
 *
 * Output:
 *   - merged_gtf: Merged transcript annotation
 *
 * Tools: StringTie (https://ccb.jhu.edu/software/stringtie/)
 */

process STRINGTIE_MERGE {
    tag "merge"
    label 'process_medium'

    publishDir "${params.outdir}/stringtie/merged", mode: 'copy'

    input:
    path(gtf_list)

    output:
    path("merged.gtf"), emit: merged_gtf

    script:
    def threads = task.cpus ?: params.threads ?: 8
    """
    # Create a list file of all GTFs for StringTie merge
    ls *.gtf > gtf_list.txt

    # Merge all transcript assemblies
    stringtie \\
        --merge \\
        -G ${params.gtf_annotation} \\
        -o merged.gtf \\
        -p ${threads} \\
        gtf_list.txt
    """
}
