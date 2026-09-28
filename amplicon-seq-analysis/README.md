# Amplicon-seq strain typing: CC-124 vs CC-1373

Code for the species-typing marker design and amplicon-sequencing analysis
described in "A high-throughput method for genomic DNA extraction from
Chlamydomonas" (xanthate lysis pub). The validated three-marker panel:

| Marker | Locus | CC-124 product | CC-1373 product | Distinguishes by |
|---|---|---|---|---|
| chr9 deletion | NC_057012.1:7,367,489–7,369,276 | 1,937 bp | ~149 bp | size (1,788 bp deletion) |
| rbcL | chloroplast NC_005353.1:123,488–123,876 | 389 bp | 389 bp | organelle presence (QC) |
| cox1 | mitochondrion NC_001638.1:5,729–6,183 | 455 bp | 455 bp | 2 SNPs (5,911 C→T; 5,929 C→A) |

## Layout

- `01_prepare_reference.py` — download the C. reinhardtii v5.5 reference
  (GCF_000002595.2) plus the organellar references (NC_005353.1, NC_001638.1);
  build the minimap2 index and BLAST database.
- `02_map_reads.sh` — map provider FASTQs with minimap2 (`-ax map-ont`),
  sort/index with samtools.
- `03_analyze_reads.py` — per-sample on-target rates (primary reads overlapping
  the intended locus ±200 bp), chr9 genotype from read depth over the deletion
  interior (midpoint NC_057012.1:7,368,383), cox1 genotype from base counts at
  the two SNPs.
- `04_compare_consensus.py` — assign each vendor consensus to its locus by
  BLASTN against the reference, and compare consensuses within a reaction by
  strand-corrected MAFFT pairwise identity.
- `design/` — how the panel was designed from public WGS data:
  - `01_download_and_align_wgs.sh` — ENA download of SRR1734609 (CC-124) and
    SRR1734611 (CC-1373), BWA-MEM alignment, mosdepth 100 bp windows.
  - `02_find_strain_deletions.py` — coverage-gap screen (≥20× CC-124, ≤5×
    CC-1373, ≥10× flanking) + bridging read-pair confirmation.
  - `03_design_typing_primers.py` — SNP-free flanks (bcftools mpileup), Primer3
    design, blastn-short uniqueness, in-silico product sizes.
  - `04_organellar_qc_primers.py` — organellar gap/indel screen (documents why
    presence/absence markers were chosen) + rbcL/cox1 primer design.
  - `05_make_amplicon_genbank.py` — annotated GenBank amplicon files for each
    marker × strain (self-validating: primers are located by exact match).

## Requirements

```
pip install -r requirements.txt
```

Plus CLI tools (versions used): minimap2 2.31, samtools 1.22.1, MAFFT 7.526,
BLAST+ 2.16.0; design phase only: bwa 0.7.19-r1273, mosdepth 0.3.14,
bcftools 1.22, curl.

## Usage

Analysis (starts from provider FASTQs named `<order>_<num>_<reaction>_<prep>_<sample>.fastq.gz`):

```
python 01_prepare_reference.py --outdir reference --email you@example.com
bash 02_map_reads.sh fastqs reference/combined_reference.mmi bams 8
python 03_analyze_reads.py --bam-dir bams --out results
python 04_compare_consensus.py --fasta-dir fastas \
    --blast-db reference/combined_reference_db --out results
```

Design (optional reproduction; downloads ~40 GB of public WGS reads):

```
bash design/01_download_and_align_wgs.sh reference/combined_reference.fna workdir 8
python design/02_find_strain_deletions.py --ref-windows workdir/mosdepth/CC124.regions.bed.gz \
    --query-windows workdir/mosdepth/CC1373.regions.bed.gz \
    --ref-bam workdir/bam/CC124.bam --query-bam workdir/bam/CC1373.bam \
    --mode nuclear --out deletions_CC1373.csv
python design/03_design_typing_primers.py --reference reference/combined_reference.fna \
    --deletions deletions_CC1373.csv --ref-bam workdir/bam/CC124.bam \
    --query-bam workdir/bam/CC1373.bam --blast-db reference/combined_reference_db \
    --out primers.tsv
python design/04_organellar_qc_primers.py --reference reference/combined_reference.fna \
    --ref-bam workdir/bam/CC124.bam --query-bam workdir/bam/CC1373.bam \
    --blast-db reference/combined_reference_db --out organellar_primers.tsv
python design/05_make_amplicon_genbank.py --reference reference/combined_reference.fna \
    --outdir amplicon_genbank
```

## Notes

- The BLAST database contains the v5.5 nuclear genome plus both organellar
  references. The designed amplicon sequences are exact subsequences of these,
  so locus assignment is identical whether or not they are added separately.
- The sequenced run also contained a chr13 candidate marker (reaction C) that
  amplified off-target; it was excluded from the panel. The analysis scripts
  still handle reaction C because the run's data contains those samples — the
  on-target table is what revealed the failure. A chr11 candidate was designed
  but not pursued.
- The analysis pipeline (`01`–`04`) was verified end-to-end against the
  provider data (order AI-8632): it reproduces the reported on-target rates,
  chr9 genotypes, cox1 SNP calls, and consensus comparisons. The `design/`
  scripts document the design method; `05` regenerates GenBank files whose
  organellar sequences are byte-identical to the originals, and the exact
  primer sequences are pinned in its marker table (Primer3 re-runs may return
  equivalent alternative pairs depending on the variant set and version).

## Data availability

- Raw QC data (NanoDrop, Qubit, TapeStation): Zenodo DOI ___ (see pub).
- Amplicon sequencing reads/consensuses (order AI-8632): ___ .
- Protocols: protocols.io (see pub).

## License

MIT (see `LICENSE`).
