#!/usr/bin/env python3
"""
prepDE.py - Prepare DESeq2 input from StringTie output

This script generates gene and transcript count matrices from StringTie
output files that can be used for differential expression analysis
with DESeq2, edgeR, or limma.

Original script from StringTie package:
https://github.com/gpertea/stringtie

Usage:
    prepDE.py -i sample_list.txt -g gene_counts.csv -t transcript_counts.csv

Input file format (tab-separated):
    sample_name<tab>path_to_gtf_or_ballgown_dir
"""

import re
import csv
import sys
import argparse
from collections import defaultdict
from operator import itemgetter

def get_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Prepare DESeq2 input from StringTie output'
    )
    parser.add_argument(
        '-i', '--input',
        metavar='input',
        help='Input file listing sample-name and path to the GTF/Ballgown file',
        required=True
    )
    parser.add_argument(
        '-g', '--gene',
        metavar='gene_count_matrix',
        default='gene_count_matrix.csv',
        help='Output gene count matrix (default: gene_count_matrix.csv)'
    )
    parser.add_argument(
        '-t', '--transcript',
        metavar='transcript_count_matrix',
        default='transcript_count_matrix.csv',
        help='Output transcript count matrix (default: transcript_count_matrix.csv)'
    )
    parser.add_argument(
        '-l', '--length',
        metavar='length',
        default=75,
        type=int,
        help='Average read length (default: 75)'
    )
    parser.add_argument(
        '-p', '--pattern',
        metavar='pattern',
        default='.',
        help='File path pattern for Ballgown ctab files (default: .)'
    )
    parser.add_argument(
        '-e', '--legend',
        metavar='legend',
        default='legend.csv',
        help='Legend file for ID to gene name mapping (optional)'
    )
    parser.add_argument(
        '-s', '--string',
        metavar='string',
        default='.',
        help='If -i is a directory, only consider files containing this string (default: .)'
    )
    parser.add_argument(
        '-c', '--cluster',
        action='store_true',
        help='Enable gene/transcript clustering'
    )
    parser.add_argument(
        '--version',
        action='version',
        version='prepDE.py 1.0.0'
    )
    return parser.parse_args()


def parse_attributes(attr_str):
    """Parse GTF attribute string into a dictionary."""
    attrs = {}
    for attr in attr_str.strip().rstrip(';').split(';'):
        attr = attr.strip()
        if not attr:
            continue
        # Handle both quoted and unquoted values
        match = re.match(r'(\S+)\s+"?([^"]*)"?', attr)
        if match:
            key, value = match.groups()
            attrs[key] = value.strip('"')
    return attrs


def get_counts_from_gtf(gtf_file, read_length):
    """Extract gene and transcript counts from a GTF file."""
    gene_counts = defaultdict(float)
    transcript_counts = defaultdict(float)
    gene_names = {}

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature = parts[2]
            attrs = parse_attributes(parts[8])

            if feature == 'transcript':
                transcript_id = attrs.get('transcript_id', '')
                gene_id = attrs.get('gene_id', '')
                gene_name = attrs.get('gene_name', gene_id)

                # Get coverage and length for count estimation
                cov = float(attrs.get('cov', 0))
                fpkm = float(attrs.get('FPKM', 0))

                # Estimate count from coverage
                # coverage = reads * read_length / transcript_length
                # So count ≈ coverage * length / read_length
                start = int(parts[3])
                end = int(parts[4])
                length = end - start + 1

                count = int(round(cov * length / read_length))

                transcript_counts[transcript_id] = count
                gene_counts[gene_id] += count
                gene_names[gene_id] = gene_name

    return gene_counts, transcript_counts, gene_names


def get_counts_from_ctab(sample_dir, read_length, pattern='.'):
    """Extract counts from Ballgown ctab files."""
    import os

    gene_counts = defaultdict(float)
    transcript_counts = defaultdict(float)
    gene_names = {}

    # Find t_data.ctab file
    t_data_file = None
    for root, dirs, files in os.walk(sample_dir):
        for f in files:
            if f == 't_data.ctab' and pattern in root:
                t_data_file = os.path.join(root, f)
                break
        if t_data_file:
            break

    if not t_data_file:
        print(f"Warning: No t_data.ctab found in {sample_dir}", file=sys.stderr)
        return gene_counts, transcript_counts, gene_names

    with open(t_data_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            transcript_id = row.get('t_name', '')
            gene_id = row.get('gene_id', '')
            gene_name = row.get('gene_name', gene_id)
            cov = float(row.get('cov', 0))
            length = int(row.get('length', 1))

            count = int(round(cov * length / read_length))

            transcript_counts[transcript_id] = count
            gene_counts[gene_id] += count
            gene_names[gene_id] = gene_name

    return gene_counts, transcript_counts, gene_names


def main():
    args = get_args()

    samples = []
    all_gene_counts = {}
    all_transcript_counts = {}
    all_genes = set()
    all_transcripts = set()
    gene_name_map = {}

    # Read sample list
    with open(args.input, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                sample_name, path = parts[0], parts[1]
            else:
                # Assume single column is just the path
                path = parts[0]
                sample_name = path.replace('.gtf', '').replace('_quant', '')

            samples.append(sample_name)

            # Determine if path is GTF or directory (Ballgown)
            import os
            if os.path.isfile(path) and path.endswith('.gtf'):
                gene_counts, transcript_counts, gene_names = get_counts_from_gtf(
                    path, args.length
                )
            elif os.path.isdir(path):
                gene_counts, transcript_counts, gene_names = get_counts_from_ctab(
                    path, args.length, args.pattern
                )
            else:
                print(f"Warning: Cannot process {path}", file=sys.stderr)
                continue

            all_gene_counts[sample_name] = gene_counts
            all_transcript_counts[sample_name] = transcript_counts
            all_genes.update(gene_counts.keys())
            all_transcripts.update(transcript_counts.keys())
            gene_name_map.update(gene_names)

    # Sort genes and transcripts
    all_genes = sorted(all_genes)
    all_transcripts = sorted(all_transcripts)

    # Write gene count matrix
    with open(args.gene, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['gene_id'] + samples)
        for gene in all_genes:
            row = [gene]
            for sample in samples:
                row.append(int(all_gene_counts.get(sample, {}).get(gene, 0)))
            writer.writerow(row)

    print(f"Gene count matrix written to {args.gene}", file=sys.stderr)

    # Write transcript count matrix
    with open(args.transcript, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['transcript_id'] + samples)
        for transcript in all_transcripts:
            row = [transcript]
            for sample in samples:
                row.append(int(all_transcript_counts.get(sample, {}).get(transcript, 0)))
            writer.writerow(row)

    print(f"Transcript count matrix written to {args.transcript}", file=sys.stderr)


if __name__ == '__main__':
    main()
