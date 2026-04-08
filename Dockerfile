# Dockerfile for bulk-rnaseq-nf pipeline
# Build: docker build -t bulk-rnaseq-nf .
# Run with Nextflow: nextflow run workflows/main.nf -params-file params.yaml -profile docker

FROM mambaorg/micromamba:1.5-jammy

LABEL maintainer="RNA-seq Pipeline Authors"
LABEL description="Docker container for bulk-rnaseq-nf Nextflow pipeline"
LABEL version="1.0.0"

# Copy environment definition
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# Install all dependencies via micromamba
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Ensure conda environment is activated
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Verify key tools are available
RUN hisat2 --version && \
    samtools --version | head -1 && \
    stringtie --version && \
    trim_galore --version && \
    fastqc --version && \
    multiqc --version && \
    featureCounts -v 2>&1 | head -1 && \
    python --version

WORKDIR /data
