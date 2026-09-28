#!/usr/bin/env python3
"""
design/05_make_amplicon_genbank.py — build annotated GenBank amplicon files.

Builds one GenBank file per marker x strain. Amplicon coordinates are found by
locating the primer sequences exactly on the reference (in-silico PCR), so the
files are self-validating. The v5.5 reference is soft-masked (lowercase), so
matching is case-insensitive. Both primer arrangements are tried: F..rc(R) on
the top strand (conventional) or R..rc(F) (amplicon on the minus strand, e.g.
rbcL). Each primer_bind feature is placed where the primer actually matches,
with the correct strand.

Records always contain the reference top-strand slice. Nuclear-marker records
include 500 bp flanks; the CC-1373 records have the deletion interval removed.
Organellar records are the bare amplicon, with the CC-1373 SNPs applied
(cox1: amp positions 183 C->T and 201 C->A; rbcL: amp position 294 A->G, called
from the WGS data at design time — note this rbcL SNP was NOT reproduced in the
later amplicon sequencing, where the two strains' rbcL products were identical).

Nuclear CC-1373 records apply the deletion only. The original design also
applied WGS-called flanking SNPs to the CC-1373 consensus; those per-strain
flanking variants are not needed for the diagnostic size shift and are left as
an optional per-marker SNP list below.

Requires: biopython.

Usage:
    python design/05_make_amplicon_genbank.py \
        --reference reference/combined_reference.fna --outdir amplicon_genbank
"""

import argparse
import os
from datetime import date

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

# Marker definitions. coords = (accession, approx. region start, region end) used
# to bound the primer search; deletion = (start, end) 1-based inclusive on the
# reference; snps = (1-based amplicon position, CC-124 base, CC-1373 base) on the
# reference top strand; organelle/misc_strand annotate the feature table.
MARKERS = [
    # Final validated panel: one nuclear deletion marker (chr9) plus the two
    # organellar QC markers. (A chr13 candidate was excluded after blastn-short
    # verification flagged its reverse primer as non-unique and the reaction
    # amplified off-target; a chr11 candidate was designed but not pursued.)
    dict(name="Marker_chr9", acc="NC_057012.1", region=(7367000, 7369800),
         deletion=(7367489, 7369276),
         F="CAAGAGACCCGTCCACGTGACG", R="AGTTCTGGCCAGGGAGGGTTGG",
         flank=500, snps=[], kind="nuclear deletion marker",
         organelle=None, misc_strand=1),
    dict(name="Cr_rbcL", acc="NC_005353.1", region=(123400, 123950),
         deletion=None,
         F="AGCCGGTGTAAAAGACTACC", R="AAGCAGGTGGAATACGAAGG",
         flank=0, snps=[(294, "A", "G")], kind="chloroplast rbcL QC",
         organelle="plastid:chloroplast", misc_strand=-1),
    dict(name="Cr_cox1", acc="NC_001638.1", region=(5650, 6250),
         deletion=None,
         F="TTCTGGTTCTTTGGTCACCC", R="AAAACGGCACCCATACTCAA",
         flank=0, snps=[(183, "C", "T"), (201, "C", "A")],
         kind="mitochondrial cox1 QC",
         organelle="mitochondrion", misc_strand=1),
]

STRAIN_LABEL = {"CC124": "CC-124 (mt-)", "CC1373": "CC-1373 (mt+)"}


def find_amplicon(ref_seqs, acc, region, f_seq, r_seq):
    """Locate the amplicon on the reference top strand; return (start, end) 1-based.

    Case-insensitive exact match (the reference is soft-masked). Tries the
    conventional arrangement (F 5' of rc(R) on the top strand) first, then the
    flipped arrangement (R 5' of rc(F)) for minus-strand amplicons such as rbcL.
    """
    seq = str(ref_seqs[acc].seq).upper()
    f_seq, r_seq = f_seq.upper(), r_seq.upper()
    lo, hi = region
    f_pos = seq.find(f_seq, lo - 1, hi)
    r_pos = seq.find(str(Seq(r_seq).reverse_complement()), lo - 1, hi)
    if f_pos != -1 and r_pos != -1 and r_pos > f_pos:
        return f_pos + 1, r_pos + len(r_seq)  # 1-based inclusive
    r_pos2 = seq.find(r_seq, lo - 1, hi)
    f_pos2 = seq.find(str(Seq(f_seq).reverse_complement()), lo - 1, hi)
    if r_pos2 != -1 and f_pos2 != -1 and f_pos2 > r_pos2:
        return r_pos2 + 1, f_pos2 + len(f_seq)
    raise ValueError(f"primers not found in {acc}:{lo}-{hi}")


def locate_primer(rec_seq, primer):
    """Return (start0, end0, strand) of an exact primer match on the record."""
    i = rec_seq.find(primer)
    if i >= 0:
        return i, i + len(primer), 1
    rc = str(Seq(primer).reverse_complement())
    j = rec_seq.find(rc)
    if j >= 0:
        return j, j + len(primer), -1
    raise ValueError(f"primer {primer} not found in record")


def apply_snps(seq, snps):
    """Apply (1-based amplicon pos, ref, alt) SNPs to an amplicon sequence."""
    s = list(seq)
    for pos, ref, alt in snps:
        if s[pos - 1].upper() != ref.upper():
            raise ValueError(f"ref mismatch at amplicon pos {pos}: "
                             f"expected {ref}, got {s[pos - 1]}")
        s[pos - 1] = alt
    return "".join(s)


def build_record(seq_str, marker, strain, amp_start, amp_end):
    """Assemble the annotated SeqRecord; primer features are located by exact match."""
    rec = SeqRecord(Seq(seq_str), id=f"{marker['name']}_{strain}",
                    name=f"{marker['name']}_{strain}",
                    description=f"{marker['kind']} amplicon, {strain}")
    rec.annotations = {
        "molecule_type": "DNA", "topology": "linear", "data_file_division": "PLN",
        "date": date.today().strftime("%d-%b-%Y").upper(),
        "organism": "Chlamydomonas reinhardtii",
        "taxonomy": ["Eukaryota", "Viridiplantae", "Chlorophyta", "Chlorophyceae",
                     "Chlamydomonadales", "Chlamydomonadaceae", "Chlamydomonas"],
        "comment": (f"Amplicon from {marker['acc']} ({marker['kind']}), reference "
                    f"coordinates {amp_start}-{amp_end} (1-based). Strain: "
                    f"{STRAIN_LABEL[strain]}. Primers: {marker['name']} F/R."),
    }
    src_qual = {"organism": ["Chlamydomonas reinhardtii"],
                "strain": [STRAIN_LABEL[strain]],
                "mol_type": ["genomic DNA"], "db_xref": ["taxon:3055"]}
    if marker["organelle"]:
        src_qual["organelle"] = [marker["organelle"]]
    rec.features.append(SeqFeature(
        FeatureLocation(0, len(seq_str), strand=1), type="source",
        qualifiers=src_qual))
    rec.features.append(SeqFeature(
        FeatureLocation(0, len(seq_str), strand=marker["misc_strand"]),
        type="misc_feature",
        qualifiers={"note": [f"{marker['kind']} amplicon ({marker['acc']})"]}))
    for label, pseq in (("F", marker["F"]), ("R", marker["R"])):
        s0, e0, strand = locate_primer(seq_str.upper(), pseq.upper())
        rec.features.append(SeqFeature(
            FeatureLocation(s0, e0, strand=strand), type="primer_bind",
            qualifiers={"label": [f"{marker['name']}_{label}"],
                        "note": [f"5'-{pseq}-3'"]}))
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--outdir", default="amplicon_genbank")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    ref_seqs = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))

    for m in MARKERS:
        amp_s, amp_e = find_amplicon(ref_seqs, m["acc"], m["region"], m["F"], m["R"])
        chrom_seq = str(ref_seqs[m["acc"]].seq).upper()
        fl = m["flank"]
        rec_s = max(1, amp_s - fl)
        rec_e = amp_e + fl
        cc124_seq = chrom_seq[rec_s - 1:rec_e]

        rec124 = build_record(cc124_seq, m, "CC124", amp_s, amp_e)

        if m["deletion"]:
            ds, de = m["deletion"]
            # remove the deletion interval (1-based inclusive) from the record
            cut_s, cut_e = ds - rec_s, de - rec_s + 1  # 0-based slice bounds
            cc1373_seq = cc124_seq[:cut_s] + cc124_seq[cut_e:]
            rec1373 = build_record(cc1373_seq, m, "CC1373", amp_s, amp_e)
            rec1373.features.append(SeqFeature(
                FeatureLocation(cut_s, cut_s, strand=1), type="misc_feature",
                qualifiers={"note": [f"CC-1373 deletion junction: {de - ds + 1} bp "
                                     f"removed ({m['acc']}:{ds}-{de})"]}))
        else:
            cc1373_seq = apply_snps(cc124_seq, m["snps"])
            rec1373 = build_record(cc1373_seq, m, "CC1373", amp_s, amp_e)
            for pos, ref, alt in m["snps"]:
                rec1373.features.append(SeqFeature(
                    FeatureLocation(pos - 1, pos, strand=1), type="variation",
                    qualifiers={"note": [f"CC-1373 SNP: {ref}>{alt} "
                                         f"(CC-124={ref}, CC-1373={alt})"]}))

        for strain, rec in (("CC124", rec124), ("CC1373", rec1373)):
            path = os.path.join(args.outdir, f"{m['name']}_{strain}.gb")
            SeqIO.write(rec, path, "genbank")
            print(f"  wrote {path} ({len(rec.seq)} bp)")


if __name__ == "__main__":
    main()
