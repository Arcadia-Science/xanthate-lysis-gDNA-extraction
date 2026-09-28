#!/usr/bin/env bash
# 02_map_reads.sh — map ONT amplicon reads to the combined reference.
#
# Usage:
#   bash 02_map_reads.sh <fastq_dir> <reference.mmi> <out_bam_dir> [threads]
#
# Expects files named <order>_<num>_<reaction>_<prep>_<sample>.fastq.gz
# (e.g. AI-8632_1_A_HF-7_cc124_A3.fastq.gz). Requires minimap2 (tested v2.31,
# map-ont preset) and samtools (tested v1.22.1) on PATH.

set -euo pipefail

FASTQ_DIR=${1:?fastq directory required}
REF=${2:?minimap2 index (.mmi) required}
OUT=${3:?output BAM directory required}
THREADS=${4:-4}

mkdir -p "$OUT"

for f in "$FASTQ_DIR"/*.fastq.gz; do
    b=$(basename "$f" .fastq.gz)
    echo "mapping $b"
    minimap2 -ax map-ont -t "$THREADS" "$REF" "$f" 2>/dev/null \
        | samtools sort -@2 -o "$OUT/$b.bam" -
    samtools index "$OUT/$b.bam"
done

echo "ALL_DONE"
