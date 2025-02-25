## Intro
**hills_bulk_RNAseq** is a bioinformatics pipeline that can be used to analyse RNA sequencing data obtained from canines and felines. It takes a samplesheet and FASTQ files as input, performs quality control (QC), trimming, alignment, and assembly. The tool produces a gene expression matrix ready for input to DESeq2 in R for differential expression analysis

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter and quality trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
3. [`HiSAT2`](https://ccb.jhu.edu/software/hisat2/index.shtml) -> **NO QUANTIFICATION**
4. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
5. Transcript assembly and quantification ([`StringTie`](https://ccb.jhu.edu/software/stringtie/))
6. Present QC for raw read, alignment, gene biotype, sample similarity, and strand-specificity checks ([`MultiQC`](http://multiqc.info/), [`R`](https://www.r-project.org/))

## Pipeline structure
~
│── conf/
│── modules/
│   ├── preprocess/
│   │   ├── concat_fastq.nf
│   │   ├── trim_galore.nf
│   ├── align/
│   │   ├── hisat2_index.nf
│   │   ├── hisat2_align.nf
│   ├── transcript_assembly/
│   │   ├── stringtie_first.nf   # First StringTie run
│   │   ├── stringtie_merge.nf   # Merge transcript assemblies
│   │   ├── stringtie_second.nf  # Second StringTie run (quantification)
│   │   ├── prepde.nf            # PrepDE.py for count matrix
│── workflows/
│   ├── main.nf
│── bin/                        # Any helper scripts
│── params.yaml                  # YAML-based input configuration
│── nextflow.config              # Pipeline-wide configuration
