#!/usr/bin/env python3
"""
design/03_design_typing_primers.py — design one shared primer pair per deletion.

For each confirmed CC-1373 deletion (from 02_find_strain_deletions.py):
  1. Extract the flanking sequence from the v5.5 reference.
  2. Call variants in the flanks for each strain (bcftools mpileup -q 20 -Q 20 |
     bcftools call -mv --ploidy 1; keep QUAL >= 20 and DP >= 5) and exclude any
     position variant in either strain -> SNP-free windows in BOTH strains.
  3. Design primers with Primer3 in the SNP-free flanks (one shared pair
     amplifies both strains; the CC-1373 product is shorter by the deletion).
  4. Verify genome uniqueness by blastn-short against the full reference
     (1 perfect hit per primer; off-targets <= 15 bp do not prime).
  5. Report in-silico product sizes for CC-124 (reference) and CC-1373
     (reference minus the deletion).

Note on reproduction: the exact pair Primer3 picks depends on the called
variant set and the primer3-py version, so a re-run may return an equivalent
SNP-free pair rather than the historical one. The final panel sequence
(RAF040/041, chr9) is fixed and verified in 05_make_amplicon_genbank.py. The
blastn-short verification step matters: a chr13 candidate pair was excluded
from the final panel after this check flagged its reverse primer as
non-unique (a second perfect site ~12 kb away on chr13), consistent with that
reaction's off-target amplification in the validation run.

Requires: primer3-py, biopython, pandas (pip); bcftools (used: 1.22), blastn
          (blastn-short, 2.16.0+) on PATH.

Usage:
    python design/03_design_typing_primers.py \
        --reference reference/combined_reference.fna \
        --deletions deletions_CC1373.csv \
        --ref-bam workdir/bam/CC124.bam --query-bam workdir/bam/CC1373.bam \
        --blast-db reference/combined_reference_db \
        --flank 500 --out primers.tsv
"""

import argparse
import subprocess
import tempfile
import os

import pandas as pd
import primer3
from Bio import SeqIO
from Bio.Seq import Seq

# Primer3 settings used for the panel (balanced Tm ~60 C, no hairpins /
# self- or pair-complementarity at the 3' end, max poly-X 4).
PRIMER3_GLOBAL_ARGS = {
    "PRIMER_TASK": "generic",
    "PRIMER_PICK_LEFT_PRIMER": 1,
    "PRIMER_PICK_RIGHT_PRIMER": 1,
    "PRIMER_PICK_INTERNAL_OLIGO": 0,
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_MAX_SIZE": 24,
    "PRIMER_OPT_TM": 60.0,
    "PRIMER_MIN_TM": 57.0,
    "PRIMER_MAX_TM": 63.0,
    "PRIMER_MIN_GC": 40.0,
    "PRIMER_MAX_GC": 65.0,
    "PRIMER_MAX_POLY_X": 4,
    "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1,
    "PRIMER_MAX_SELF_ANY_TH": 45.0,
    "PRIMER_MAX_SELF_END_TH": 35.0,
    "PRIMER_MAX_HAIRPIN_TH": 24.0,
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_DNA_CONC": 250.0,
    "PRIMER_NUM_RETURN": 5,
}


def call_flank_variants(bam, ref_fasta, chrom, start, end, qual_min=20, dp_min=5):
    """Variant positions (1-based) in a region for one strain."""
    region = f"{chrom}:{start}-{end}"
    mp = subprocess.run(["bcftools", "mpileup", "-f", ref_fasta, "-r", region,
                         "-q", "20", "-Q", "20", bam],
                        capture_output=True, text=True, check=True)
    call = subprocess.run(["bcftools", "call", "-mv", "--ploidy", "1"],
                          input=mp.stdout, capture_output=True, text=True, check=True)
    positions = set()
    for line in call.stdout.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        qual = float(f[5]) if f[5] != "." else 0.0
        dp = 0
        for kv in f[7].split(";"):
            if kv.startswith("DP="):
                dp = int(kv[3:])
        if qual >= qual_min and dp >= dp_min:
            positions.add(int(f[1]))
    return positions


def snp_free_windows(start, end, variant_positions, min_len=25):
    """Contiguous [start, end) windows (1-based) containing no variant position."""
    windows = []
    w_start = start
    for p in sorted(variant_positions):
        if p - w_start >= min_len:
            windows.append((w_start, p - 1))
        w_start = max(w_start, p + 1)
    if end - w_start + 1 >= min_len:
        windows.append((w_start, end))
    return windows


def blastn_short_uniqueness(primers, db):
    """blastn-short each primer; return {name: (n_perfect_hits, max_offtarget_aln_len)}."""
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tf:
        for name, seq in primers.items():
            tf.write(f">{name}\n{seq}\n")
        path = tf.name
    out = subprocess.run(
        ["blastn", "-task", "blastn-short", "-query", path, "-db", db,
         "-evalue", "1", "-word_size", "7",
         "-outfmt", "6 qseqid sseqid pident length"],
        capture_output=True, text=True, check=True).stdout
    os.unlink(path)
    stats = {}
    for line in out.splitlines():
        q, _, pident, length = line.split("\t")[:4]
        perfect = float(pident) == 100.0 and int(length) == len(primers[q])
        n_perfect, max_off = stats.get(q, (0, 0))
        if perfect:
            n_perfect += 1
        else:
            max_off = max(max_off, int(length))
        stats[q] = (n_perfect, max_off)
    return stats


def insilico_pcr_size(ref_seqs, chrom, f_seq, r_seq):
    """Exact-match both primers on the chromosome; return (start, end, size) or None.

    Case-insensitive: the v5.5 reference is soft-masked (lowercase repeats).
    """
    seq = str(ref_seqs[chrom].seq).upper()
    f_pos = seq.find(f_seq.upper())
    r_pos = seq.find(str(Seq(r_seq.upper()).reverse_complement()))
    if f_pos == -1 or r_pos == -1 or r_pos <= f_pos:
        return None
    return (f_pos + 1, r_pos + len(r_seq), r_pos + len(r_seq) - f_pos)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--deletions", required=True, help="CSV from 02_find_strain_deletions.py")
    ap.add_argument("--ref-bam", required=True)
    ap.add_argument("--query-bam", required=True)
    ap.add_argument("--blast-db", required=True)
    ap.add_argument("--flank", type=int, default=500)
    ap.add_argument("--product-range", default="1500-2500",
                    help="CC-124 product size range, e.g. 1500-2500")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ref_seqs = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))
    deletions = pd.read_csv(args.deletions)
    if "confirmed" in deletions.columns:
        deletions = deletions[deletions.confirmed]
    pmin, pmax = [int(x) for x in args.product_range.split("-")]

    rows = []
    for _, d in deletions.iterrows():
        chrom, ds, de = d.chrom, int(d.start), int(d.end)
        fs, fe = max(1, ds - args.flank), de + args.flank
        print(f"\n=== {chrom}:{ds}-{de} ({de - ds} bp deletion) ===")

        var = (call_flank_variants(args.ref_bam, args.reference, chrom, fs, fe) |
               call_flank_variants(args.query_bam, args.reference, chrom, fs, fe))
        left_windows = snp_free_windows(fs, ds - 1, var)
        right_windows = snp_free_windows(de + 1, fe, var)
        if not left_windows or not right_windows:
            print("  no SNP-free flank windows; skipping")
            continue
        lw = max(left_windows, key=lambda w: w[1] - w[0])
        rw = max(right_windows, key=lambda w: w[1] - w[0])

        template = str(ref_seqs[chrom].seq[fs - 1:fe])  # 0-based slice
        res = primer3.design_primers(
            seq_args={
                "SEQUENCE_ID": f"{chrom}_{ds}",
                "SEQUENCE_TEMPLATE": template,
                # Primer3 OK regions are 0-based relative to template
                "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [
                    lw[0] - fs, lw[1] - lw[0] + 1,
                    rw[0] - fs, rw[1] - rw[0] + 1],
            },
            global_args={**PRIMER3_GLOBAL_ARGS,
                         "PRIMER_PRODUCT_SIZE_RANGE": [[pmin, pmax]]})
        if "PRIMER_LEFT_0_SEQUENCE" not in res:
            print("  Primer3 found no pair; skipping")
            continue

        f_seq = res["PRIMER_LEFT_0_SEQUENCE"]
        r_seq = res["PRIMER_RIGHT_0_SEQUENCE"]
        uniq = blastn_short_uniqueness({"F": f_seq, "R": r_seq}, args.blast_db)
        pcr = insilico_pcr_size(ref_seqs, chrom, f_seq, r_seq)
        cc124_size = pcr[2] if pcr else None
        cc1373_size = cc124_size - (de - ds) if cc124_size else None
        rows.append(dict(chrom=chrom, del_start=ds, del_end=de, deletion_bp=de - ds,
                         primer_F=f_seq, primer_R=r_seq,
                         F_tm=round(res["PRIMER_LEFT_0_TM"], 1),
                         R_tm=round(res["PRIMER_RIGHT_0_TM"], 1),
                         F_perfect_hits=uniq.get("F", (0, 0))[0],
                         R_perfect_hits=uniq.get("R", (0, 0))[0],
                         cc124_product=cc124_size, cc1373_product=cc1373_size))
        print(f"  F {f_seq}  R {r_seq}")
        print(f"  uniqueness (perfect hits): F={uniq.get('F')}, R={uniq.get('R')}")
        print(f"  in-silico product: CC-124={cc124_size} bp, CC-1373={cc1373_size} bp")

    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"\nWritten to {args.out}")


if __name__ == "__main__":
    main()
