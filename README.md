# bulk-rnaseq-nf

A Nextflow pipeline for bulk RNA sequencing analysis.

## Overview

**bulk-rnaseq-nf** is a bioinformatics pipeline that performs comprehensive analysis of bulk RNA sequencing data. It takes FASTQ files as input and produces gene/transcript count matrices suitable for differential expression analysis with tools like DESeq2, edgeR, or limma.

## Pipeline Steps

1. **Lane concatenation** - Merge FASTQ files for samples sequenced across multiple lanes
2. **Quality control & trimming** - Adapter trimming and QC with [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
3. **Index generation** - Build [HISAT2](https://ccb.jhu.edu/software/hisat2/) index (optional, if not pre-built)
4. **Alignment** - Map reads to reference genome with [HISAT2](https://ccb.jhu.edu/software/hisat2/)
5. **BAM processing** - Sort and index alignments with [SAMtools](http://www.htslib.org/)
6. **Transcript assembly** - Assemble transcripts with [StringTie](https://ccb.jhu.edu/software/stringtie/) (first pass)
7. **Merge assemblies** - Create unified transcript annotation across samples
8. **Quantification** - Estimate abundances with [StringTie](https://ccb.jhu.edu/software/stringtie/) (second pass)
9. **Count matrices** - Generate gene/transcript counts with prepDE.py
10. **QC report** - Aggregate QC metrics with [MultiQC](https://multiqc.info/)

## Requirements

### System Requirements

- Linux or macOS
- Nextflow >= 21.10.0
- Java >= 11

### Software Dependencies

Install via conda (recommended):

```bash
conda env create -f environment.yml
conda activate bulk-rnaseq
```

Or ensure the following tools are in your PATH:
- trim_galore >= 0.6.7
- fastqc >= 0.11.9
- hisat2 >= 2.2.1
- samtools >= 1.15
- stringtie >= 2.2.1
- multiqc >= 1.12
- python >= 3.8

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/bulk-rnaseq-nf.git
cd bulk-rnaseq-nf

# Install dependencies with conda
conda env create -f environment.yml
conda activate bulk-rnaseq

# Verify installation
nextflow run workflows/main.nf --help
```

## Quick Start

1. **Prepare your parameter file**

   Edit `params.yaml` with your sample information and reference paths:

   ```yaml
   samples:
     sample1:
       R1: "data/sample1_R1.fastq.gz"
       R2: "data/sample1_R2.fastq.gz"
     sample2:
       R1: "data/sample2_R1.fastq.gz"
       R2: "data/sample2_R2.fastq.gz"

   hisat2_index: "/path/to/hisat2/index"
   gtf_annotation: "/path/to/annotation.gtf"
   outdir: "results"
   ```

2. **Run the pipeline**

   ```bash
   nextflow run workflows/main.nf -params-file params.yaml
   ```

## Usage

### Basic Usage

```bash
nextflow run workflows/main.nf -params-file params.yaml
```

### With Profile

```bash
# Using conda environment
nextflow run workflows/main.nf -params-file params.yaml -profile conda

# Using Docker containers
nextflow run workflows/main.nf -params-file params.yaml -profile docker

# Using Singularity containers
nextflow run workflows/main.nf -params-file params.yaml -profile singularity

# On SLURM cluster
nextflow run workflows/main.nf -params-file params.yaml -profile slurm
```

### Resume Failed Run

```bash
nextflow run workflows/main.nf -params-file params.yaml -resume
```

## Parameters

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `samples` | Sample definitions with R1/R2 FASTQ paths |
| `gtf_annotation` | Reference GTF annotation file |
| `hisat2_index` OR `genome_fasta` | Pre-built HISAT2 index or genome FASTA to build one |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `outdir` | `results` | Output directory |
| `threads` | `4` | Default threads per process |
| `max_cpus` | `16` | Maximum CPUs |
| `max_memory` | `64.GB` | Maximum memory |
| `max_time` | `48.h` | Maximum time per process |

## Output Structure

```
results/
├── trimmed/              # Trimmed FASTQ files
├── fastqc/               # FastQC reports
├── hisat2/               # Alignment logs
├── aligned/              # Sorted BAM files and indices
├── stringtie/
│   ├── first_pass/       # Initial transcript assemblies
│   ├── merged/           # Merged annotation
│   └── quantification/   # Quantified transcripts
├── ballgown/             # Ballgown table files
├── counts/               # Count matrices for DE analysis
│   ├── gene_count_matrix.csv
│   └── transcript_count_matrix.csv
├── multiqc/              # MultiQC report
└── pipeline_info/        # Execution reports
    ├── timeline.html
    ├── report.html
    ├── trace.txt
    └── dag.svg
```

## Project Structure

```
bulk-rnaseq-nf/
├── workflows/
│   └── main.nf           # Main pipeline workflow
├── modules/
│   ├── preprocess/       # Concatenation and trimming
│   ├── align/            # HISAT2 indexing and alignment
│   ├── transcript_assembly/  # StringTie assembly and quantification
│   └── qc/               # MultiQC reporting
├── bin/
│   └── prepDE.py         # Count matrix generation script
├── nextflow.config       # Pipeline configuration
├── params.yaml           # Parameter template
├── environment.yml       # Conda environment
├── CITATIONS.md          # Tool citations
└── README.md             # This file
```

## Troubleshooting

### Common Issues

**Pipeline fails with "No such file" error**
- Ensure all FASTQ paths in `params.yaml` are correct and files exist
- Use absolute paths for reliability

**Out of memory errors**
- Reduce `max_memory` in params or use a profile with more resources
- For large genomes, HISAT2 indexing requires ~200GB RAM

**HISAT2 alignment fails**
- Ensure HISAT2 index was built with the same version
- Check that index files are not corrupted

### Getting Help

1. Check the [Nextflow documentation](https://www.nextflow.io/docs/latest/)
2. Review error messages in `.nextflow.log`
3. Check process logs in `work/` directory

## Citations

See [CITATIONS.md](CITATIONS.md) for the complete list of tools and their citations.

## License

This project is licensed under the MIT License.
