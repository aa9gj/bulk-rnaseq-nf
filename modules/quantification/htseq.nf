/*
 * HTSEQ_COUNT
 *
 * Count reads mapping to genomic features using HTSeq-count.
 * HTSeq is widely used and produces results compatible with
 * DESeq2 and edgeR.
 *
 * Note: HTSeq processes one BAM at a time, so this runs per-sample.
 *
 * Input:
 *   - sample_id: Sample identifier
 *   - bam: Sorted BAM file
 *   - bai: BAM index file
 *   - gtf: GTF annotation file
 *
 * Output:
 *   - counts: Per-sample gene counts
 *
 * Tools: HTSeq (https://htseq.readthedocs.io/)
 */

process HTSEQ_COUNT {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/htseq", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(gtf)

    output:
    path("${sample_id}.htseq.txt"), emit: counts

    script:
    // Strandedness: yes, no, reverse
    def strand = params.strandedness == 1 ? 'yes' : (params.strandedness == 2 ? 'reverse' : 'no')
    """
    htseq-count \\
        -f bam \\
        -r pos \\
        -s ${strand} \\
        -t exon \\
        -i gene_id \\
        ${bam} \\
        ${gtf} > ${sample_id}.htseq.txt
    """
}

/*
 * HTSEQ_MERGE
 *
 * Merge individual HTSeq count files into a single count matrix.
 *
 * Input:
 *   - count_files: Collection of per-sample HTSeq count files
 *
 * Output:
 *   - merged_counts: Combined count matrix
 */

process HTSEQ_MERGE {
    tag "merge_htseq"
    label 'process_low'

    publishDir "${params.outdir}/htseq", mode: 'copy'

    input:
    path(count_files)

    output:
    path("htseq_counts.txt"), emit: merged_counts

    script:
    """
    # Get gene IDs from first file (excluding summary lines starting with __)
    head -1 *.htseq.txt | cut -f1 | grep -v "^__" > genes.txt || true

    # Create header
    echo -n "gene_id" > htseq_counts.txt
    for f in *.htseq.txt; do
        sample=\$(basename \$f .htseq.txt)
        echo -ne "\\t\${sample}" >> htseq_counts.txt
    done
    echo "" >> htseq_counts.txt

    # Merge counts (exclude summary lines)
    paste *.htseq.txt | grep -v "^__" | awk '{
        printf \$1
        for (i=2; i<=NF; i+=2) printf "\\t"\$i
        print ""
    }' >> htseq_counts.txt
    """
}
