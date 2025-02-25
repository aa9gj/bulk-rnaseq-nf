include { CONCAT_FASTQ } from './modules/preprocess/concat_fastq.nf'
include { TRIM_GALORE } from './modules/preprocess/trim_galore.nf'
include { HISAT2_INDEX } from './modules/align/hisat2_index.nf'
include { HISAT2_ALIGN } from './modules/align/hisat2_align.nf'
include { STRINGTIE_FIRST } from './modules/transcript_assembly/stringtie_first.nf'
include { STRINGTIE_MERGE } from './modules/transcript_assembly/stringtie_merge.nf'
include { STRINGTIE_SECOND } from './modules/transcript_assembly/stringtie_second.nf'
include { PREPDE } from './modules/transcript_assembly/prepde.nf'

workflow {
    samples_ch = Channel.fromPath(params.yaml)
        .map { key, val -> tuple(key, val.R1.tokenize(','), val.R2.tokenize(',')) }

    concatenated = samples_ch | CONCAT_FASTQ
    trimmed = concatenated | TRIM_GALORE

    hisat2_index = params.hisat2_index ? params.hisat2_index : HISAT2_INDEX(params.genome_fasta)
   
    aligned = trimmed.combine(hisat2_index) | HISAT2_ALIGN
   
    assembled_transcripts = aligned.combine(params.gtf_annotation) | STRINGTIE_FIRST

    merged_gtf = assembled_transcripts.collect() | STRINGTIE_MERGE

    quantified = aligned.combine(merged_gtf) | STRINGTIE_SECOND

    prepde = quantified.collect() | PREPDE
}
