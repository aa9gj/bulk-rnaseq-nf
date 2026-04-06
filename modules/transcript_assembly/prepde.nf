/*
 * PREPDE
 *
 * Generate count matrices from StringTie output using prepDE.py.
 * Creates gene-level and transcript-level count matrices suitable
 * for differential expression analysis with DESeq2, edgeR, or limma.
 *
 * Input:
 *   - gtf_files: Collection of quantified GTF files from STRINGTIE_SECOND
 *
 * Output:
 *   - gene_counts: Gene-level count matrix (CSV)
 *   - transcript_counts: Transcript-level count matrix (CSV)
 *
 * Note: Uses the real prepDE.py that ships with the StringTie package.
 *
 * Tools: prepDE.py (bundled with StringTie)
 */

process PREPDE {
    tag "count_matrix"
    label 'process_low'

    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path(gtf_files)

    output:
    path("gene_count_matrix.csv"), emit: gene_counts
    path("transcript_count_matrix.csv"), emit: transcript_counts
    path("sample_list.txt"), emit: sample_list

    script:
    """
    # Create sample list file for prepDE.py
    # Format: sample_name<tab>path_to_gtf
    for gtf in *.gtf; do
        sample=\$(basename \${gtf} _quant.gtf)
        echo -e "\${sample}\\t\${gtf}"
    done > sample_list.txt

    # Generate count matrices
    prepDE.py3 \\
        -i sample_list.txt \\
        -g gene_count_matrix.csv \\
        -t transcript_count_matrix.csv
    """
}
