#!/usr/bin/env python3
"""
04_nanodrop_qubit.py — NanoDrop and Qubit per-sample values.

Per sample in metadata/samples.csv with a NanoDrop reading:

  nd_ng_per_ul      A260 × 50 ng/µL × dilution factor (the dsDNA conversion)
  a260_a280, a260_a230
  qubit_dsdna       Qubit 1× dsDNA HS, original sample concentration (ng/µL)
  nd_qubit_ratio    nd_ng_per_ul / qubit_dsdna
  qubit_rna         Qubit RNA HS (ng/µL) where measured
  qubit_rna_status  'ok', or 'below_range' when the instrument reported
                    "Out of range" with a sample RFU below its top standard
                    (at or below the detection limit; value left empty and
                    excluded from means), or 'above_range'

It also writes each sample's spectrum (220–320 nm), divided by its own
absorbance at 260 nm, to nanodrop_spectra_normalized.csv (unsmoothed; the
plotting script applies a 9-point display smoothing).

Feeds Figure 2 D–F, Figure 3 D–G and Figure 4D.

Requires: numpy, pandas (see requirements.txt).

Usage:
    python 04_nanodrop_qubit.py --data-dir /path/to/gDNA_extractions --out results
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import qclib as q


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    samples = q.load_samples()
    rows, spectra = [], []
    for sid, s in samples[samples["nanodrop_file"] != ""].iterrows():
        nd = q.read_nanodrop(str(q.find_file(args.data_dir, s["prep"], s["nanodrop_file"])))
        rec = nd[s["nanodrop_time"]]
        conc = rec["A260"] * 50 * float(s["nanodrop_dilution"])
        r = dict(sample_id=sid, prep=s["prep"], nanodrop_time=s["nanodrop_time"], a260=rec["A260"],
                 dilution=float(s["nanodrop_dilution"]), nd_ng_per_ul=conc,
                 a260_a280=rec["A260/A280"], a260_a230=rec["A260/A230"])
        if s["qubit_file"]:
            qd = q.qubit(args.data_dir, s["prep"], s["qubit_file"], s["qubit_key"])
            if qd["status"] != "ok":
                raise ValueError(f"{sid}: Qubit dsDNA reading is {qd['status']}")
            r.update(qubit_dsdna=qd["conc"], nd_qubit_ratio=conc / qd["conc"],
                     qubit_source="notebook transcription" if "from-Notion" in s["qubit_file"]
                     else "instrument export")
        if s["qubit_rna_file"]:
            qr = q.qubit(args.data_dir, s["prep"], s["qubit_rna_file"], s["qubit_rna_key"])
            r.update(qubit_rna=qr["conc"], qubit_rna_status=qr["status"], qubit_rna_rfu=qr["rfu"],
                     qubit_rna_std1_rfu=qr["std1_rfu"])
        rows.append(r)
        w, a = q.spectrum_normalized(rec)
        spectra += [dict(sample_id=sid, wavelength_nm=x, abs_norm_a260=y) for x, y in zip(w, a)]

    res = pd.DataFrame(rows)
    res.to_csv(out / "nanodrop_qubit.csv", index=False, float_format="%.4f")
    pd.DataFrame(spectra).to_csv(out / "nanodrop_spectra_normalized.csv", index=False,
                                 float_format="%.5f")
    below = res.get("qubit_rna_status", pd.Series(dtype=str)).eq("below_range").sum()
    print(f"{len(res)} NanoDrop readings -> {out}/nanodrop_qubit.csv "
          f"({int(below)} Qubit RNA reads at or below the detection limit)")


if __name__ == "__main__":
    main()
