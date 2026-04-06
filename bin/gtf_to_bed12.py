#!/usr/bin/env python3
"""
gtf_to_bed12.py - Convert GTF annotation to BED12 format

RSeQC tools require a BED12 gene model file. This script converts
a standard GTF file (Ensembl/GENCODE format) to BED12 format.

Usage:
    gtf_to_bed12.py -i annotation.gtf -o annotation.bed
"""

import argparse
import sys
from collections import defaultdict


def parse_attributes(attr_str):
    """Parse GTF attribute column into a dictionary."""
    attrs = {}
    for attr in attr_str.strip().rstrip(';').split(';'):
        attr = attr.strip()
        if not attr:
            continue
        parts = attr.split(' ', 1)
        if len(parts) == 2:
            key = parts[0]
            value = parts[1].strip('"').strip()
            attrs[key] = value
    return attrs


def gtf_to_bed12(gtf_file, bed_file):
    """Convert GTF to BED12 format."""
    # Collect exons per transcript
    transcripts = {}
    transcript_info = {}

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom = fields[0]
            feature = fields[2]
            start = int(fields[3]) - 1  # Convert to 0-based
            end = int(fields[4])
            strand = fields[6]
            attrs = parse_attributes(fields[8])

            transcript_id = attrs.get('transcript_id', '')
            if not transcript_id:
                continue

            if feature == 'transcript':
                transcript_info[transcript_id] = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'gene_name': attrs.get('gene_name', attrs.get('gene_id', transcript_id)),
                }
            elif feature == 'exon':
                if transcript_id not in transcripts:
                    transcripts[transcript_id] = []
                transcripts[transcript_id].append((start, end))

                # If we haven't seen a transcript feature, build info from exons
                if transcript_id not in transcript_info:
                    transcript_info[transcript_id] = {
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'gene_name': attrs.get('gene_name', attrs.get('gene_id', transcript_id)),
                    }
                else:
                    # Update boundaries
                    info = transcript_info[transcript_id]
                    info['start'] = min(info['start'], start)
                    info['end'] = max(info['end'], end)

    # Write BED12 output
    with open(bed_file, 'w') as out:
        for tid, exons in sorted(transcripts.items()):
            if tid not in transcript_info:
                continue

            info = transcript_info[tid]
            chrom = info['chrom']
            tx_start = info['start']
            tx_end = info['end']
            strand = info['strand']
            name = info['gene_name']

            # Sort exons by start position
            exons.sort(key=lambda x: x[0])

            block_count = len(exons)
            block_sizes = ','.join(str(e[1] - e[0]) for e in exons)
            block_starts = ','.join(str(e[0] - tx_start) for e in exons)

            # BED12: chrom, start, end, name, score, strand,
            #        thickStart, thickEnd, rgb, blockCount, blockSizes, blockStarts
            out.write(f"{chrom}\t{tx_start}\t{tx_end}\t{name}\t0\t{strand}\t"
                      f"{tx_start}\t{tx_end}\t0\t{block_count}\t{block_sizes}\t{block_starts}\n")

    print(f"Converted {len(transcripts)} transcripts to BED12 format", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='Convert GTF to BED12 format for RSeQC')
    parser.add_argument('-i', '--input', required=True, help='Input GTF file')
    parser.add_argument('-o', '--output', required=True, help='Output BED12 file')
    args = parser.parse_args()

    gtf_to_bed12(args.input, args.output)


if __name__ == '__main__':
    main()
