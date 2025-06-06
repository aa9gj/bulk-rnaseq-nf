## Introduction
**bulk-rnaseq-nf** is a bioinformatics pipeline that can be used to analyse RNA sequencing data. It takes a samplesheet and FASTQ files as input, performs lane concatenation, quality control (QC), trimming, alignment, assembly, quantification, and prepares data for input into packages (e.g. DESeq2) for differential expression analysis.

1. Lane concatenation for samples sequences on multiple lanes
2. Adapter trimming, and read QC ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
3. [`HiSAT2`](https://ccb.jhu.edu/software/hisat2/index.shtml) index generation if not readily available 
4. [`HiSAT2`](https://ccb.jhu.edu/software/hisat2/index.shtml) alignment 
5. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
6. Transcript assembly and quantification ([`StringTie`](https://ccb.jhu.edu/software/stringtie/))
7. Present QC for raw read, alignment, gene biotype, sample similarity, and strand-specificity checks ([`MultiQC`](http://multiqc.info/), [`R`](https://www.r-project.org/))

## Pipeline structure
Each directory and file is structured to facilitate the processing pipeline.

- `conf/`: Configuration files related to the project.
- `modules/`: Contains sub-modules, each serving specific roles like preprocessing, alignment, and transcript assembly.
  - `preprocess/`: Preprocessing scripts, such as concatenating and trimming fastqs.
  - `align/`: Contains scripts for indexing and alignment using HISAT2.
  - `transcript_assembly/`: Scripts for transcript assembly and quantification using StringTie.
  - `qc`: Quality control 
- `workflows/`: Main pipeline scripts.
- `bin/`: Directory for helper scripts.
- `params.yaml`: Configuration file specifying input parameters.
- `nextflow.config`: Pipeline-wide configuration settings.

## Pipeline usage
```sh
nextflow run workflows/main.nf -params-file params.yaml
```
