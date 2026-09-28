#!/usr/bin/env bash
# design/01_download_and_align_wgs.sh — download public WGS reads and align to the reference.
#
# Downloads whole-genome resequencing reads for CC-124 (SRR1734609) and CC-1373
# (SRR1734611) (Flowers et al. 2015, BioProject PRJNA271103) from ENA, aligns them
# to the C. reinhardtii v5.5 reference (GCF_000002595.2; plus organellar references
# NC_005353.1 / NC_001638.1 if present in the FASTA) with BWA-MEM, and computes
# per-window coverage with mosdepth (100 bp windows). Expect ~90x coverage per strain.
#
# Requires: bwa (used: 0.7.19-r1273), samtools (1.22.1), mosdepth (0.3.14), curl.
# Usage:
#   bash design/01_download_and_align_wgs.sh <reference.fna> <workdir> [threads]

set -euo pipefail

REF=${1:?reference FASTA required (see 01_prepare_reference.py)}
WORK=${2:?workdir required}
THREADS=${3:-4}

mkdir -p "$WORK"/{fastq,bam,mosdepth}

# ENA direct FASTQ download (no SRA toolkit needed).
# SRR1734609 = CC-124, SRR1734611 = CC-1373.
declare -A R1=( [CC124]=https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR173/009/SRR1734609/SRR1734609_1.fastq.gz
                [CC1373]=https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR173/001/SRR1734611/SRR1734611_1.fastq.gz )
declare -A R2=( [CC124]=https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR173/009/SRR1734609/SRR1734609_2.fastq.gz
                [CC1373]=https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR173/001/SRR1734611/SRR1734611_2.fastq.gz )

for strain in CC124 CC1373; do
    r1="$WORK/fastq/${strain}_R1.fastq.gz"
    r2="$WORK/fastq/${strain}_R2.fastq.gz"
    [ -f "$r1" ] || curl -fSL --retry 3 -o "$r1" "${R1[$strain]}"
    # _2 may not exist for single-end libraries; tolerate failure
    [ -f "$r2" ] || curl -fSL --retry 3 -o "$r2" "${R2[$strain]}" || rm -f "$r2"

    bam="$WORK/bam/${strain}.bam"
    if [ ! -f "$bam" ]; then
        if [ -f "$r2" ]; then
            bwa mem -t "$THREADS" "$REF" "$r1" "$r2" | samtools sort -@2 -o "$bam" -
        else
            bwa mem -t "$THREADS" "$REF" "$r1" | samtools sort -@2 -o "$bam" -
        fi
        samtools index "$bam"
    fi

    # 100 bp non-overlapping per-window mean coverage
    mosdepth -n -t "$THREADS" --by 100 "$WORK/mosdepth/${strain}" "$bam"
done

echo "DONE: BAMs in $WORK/bam, mosdepth windows in $WORK/mosdepth"
