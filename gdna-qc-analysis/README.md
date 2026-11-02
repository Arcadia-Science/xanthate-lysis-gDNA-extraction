# Genomic DNA QC analysis: TapeStation, NanoDrop and Qubit

Code for the quality-control analysis and figure data panels in "A
high-throughput method for genomic DNA extraction from Chlamydomonas". 
Every number is computed from raw instrument exports:
TapeStation electropherogram traces sized against each run's own ladder,
NanoDrop spectra, and Qubit readings. No vendor-computed metric (peak tables,
DIN, calibrated concentrations) is used.

## Layout

- `qclib.py` — shared readers (TapeStation, NanoDrop, Qubit), ladder
  calibration, trace metrics and summary helpers.
- `01_tapestation_metrics.py` — size every Genomic DNA ScreenTape lane against
  its run's ladder (14 bands, 100 bp–48.5 kb; PCHIP in log10(size), linear in
  log10(size) beyond the end bands) and compute average DNA size and the
  percentage of DNA at or above 10 kb from signal between 400 bp and 1 Mb.
  Also reports the percentage of that signal above the ladder's top band
  (48.5 kb), which is sized by extrapolation alone, and integrates each lane's
  signal at or above 10 kb and its 100 bp lower marker (50–150 bp).
- `02_hsd5000_low_mw.py` — High Sensitivity D5000 check that nothing is hidden
  under the Genomic DNA ScreenTape's lower marker: share of lane signal,
  excluding both HS D5000 markers (below 60 bp; 9–11 kb), at 70–150 bp and
  60–400 bp.
- `03_lower_marker_comparison.py` — signal at or above 10 kb divided by the
  lower marker, compared between the "RNA-free" and "Single-SPRI" cleanups.
- `04_nanodrop_qubit.py` — NanoDrop concentration (A260 × 50 ng/µL × dilution),
  A260/A280, A260/A230, Qubit dsDNA and RNA, the NanoDrop:Qubit ratio, and
  spectra normalized to each sample's own absorbance at 260 nm.
- `05_summary_stats.py` — every mean ± SD, range and n quoted in "Observed
  performance", "Optimization notes" and the figure legends, checked against
  the values as written in the pub (`metadata/pub_values.csv`).
- `06_plot_figures.py` — data panels for Figures 2, 3 and 4. The published
  figures were laid out in Illustrator from panels like these.
- `metadata/samples.csv` — one row per sample: which files and which
  lane/reading hold its data (see "Sample sheet" below).
- `metadata/sets.csv` — which samples make up each group in each figure.
- `metadata/pub_values.csv` — the values as written in the pub, for `05`.

## Figure panels and preparations

| Set (`sets.csv`) | Figure panels | Preparation(s) | Groups |
|---|---|---|---|
| `hf_track` | Fig 2 A–G | HF-8 | *C. reinhardtii* CC-124, *C. smithii* CC-1373 (n = 6 each) |
| `ht_track` | Fig 3 A–G | HT-12 | "RNA-free" (row D), "Single-SPRI" (row H), n = 12 each |
| `opt_handling` | Fig 4A | HF-1 (before), HF-3 (after) | n = 4, 6 |
| `opt_timecourse` | Fig 4B | HT-6 | 0, 1, 5, 10, 20, 40 min lysis; CC-124, n = 2 each |
| `opt_thaw_ht` | Fig 4C (high throughput) | HT-6 10 min (hot dry); HT-7 and HT-8 (cold wet) | CC-124, n = 2, 7 |
| `opt_thaw_hf` | Fig 4C (high fidelity) | HF-5 (hot dry); HF-7 (cold wet) | CC-124, n = 3, 7 |
| `opt_additives` | Fig 4D | HT-5 | lysis buffer, + proteinase K, + proteinase K + RNase A/T1; 5 matched cultures |

Sample selection, recorded per row in `samples.csv` (`exclude_reason`):

- Figure 4 B and C are restricted to *C. reinhardtii* CC-124. Deep-well lanes
  labelled cc1373 are excluded because that label was unreliable; the three
  HF-7 lanes labelled cc1373 (G4, H7, H12) were typed as CC-124 by chr9 PCR and
  amplicon sequencing and are included.
- Figure 4D drops HT-5 cc1373_1, a swimming flask culture unlike the
  palmelloid *C. smithii* stock used elsewhere (n = 5 cultures: three
  *C. reinhardtii*, two *C. smithii*).

## Data

Raw QC data are on Zenodo ([DOI: 10.5281/zenodo.23811781](https://doi.org/10.5281/zenodo.23811781)). The archive mirrors the lab's
`gDNA_extractions/` folder: one `<prep>_<date>/` folder per preparation with
its TapeStation, NanoDrop and Qubit exports. Files are found by name anywhere
under their preparation folder, so sub-folder names don't matter.

- TapeStation: `<run>_Electropherogram.csv` (Genomic DNA and HS D5000
  ScreenTapes), exported from TapeStation Analysis Software 4.1.1.
- NanoDrop One: the long-format `dsDNA <date> <time>.tsv` export (the wide
  `_table.tsv` export is column-misaligned and isn't used). HF-8 was read on the
  1:10 dilution prepared for Qubit, so its dilution factor is 10.
- Qubit 4: `QubitData_<date>_<time>.csv` instrument exports (HT-12). The
  exports carry no sample names, so `samples.csv` maps each well to its
  reading by acquisition time. For HF-8 and HT-5 no instrument export was
  kept; `<prep>_Qubit_from-Notion.csv` holds the values transcribed from the
  lab notebook.
- Qubit RNA HS "Out of range" readings had sample fluorescence at the 0 ng/µL
  standard (below the top standard), so they are reported as at or below the
  detection limit, counted rather than averaged, and drawn at 0 in Figure 3F.

## Sample sheet

`metadata/samples.csv` columns: `sample_id`, `prep`, `sample`, `species`,
`strain`; `tapestation_file`, `tapestation_well`, `tapestation_desc` (checked
against the export header); `hsd5000_file`, `hsd5000_well`; `nanodrop_file`,
`nanodrop_time` (the measurement's timestamp in the export), `nanodrop_dilution`;
`qubit_file`, `qubit_key` (Test Date for instrument exports, sample name for
transcriptions); `qubit_rna_file`, `qubit_rna_key`; `exclude_reason`, `notes`.

## Requirements

```
pip install -r requirements.txt
```

Tested with Python 3.11, numpy 2.4.4, pandas 3.0.2, scipy 1.17.1, matplotlib
3.10.9 and arcadia-pycolor 0.8.0.

## Usage

```
DATA=/path/to/gDNA_extractions        # unpacked Zenodo archive
python 01_tapestation_metrics.py --data-dir $DATA --out results
python 02_hsd5000_low_mw.py      --data-dir $DATA --out results
python 03_lower_marker_comparison.py --results results
python 04_nanodrop_qubit.py      --data-dir $DATA --out results
python 05_summary_stats.py       --results results
python 06_plot_figures.py        --data-dir $DATA --results results --out figures
```

Outputs in `results/`: `tapestation_metrics.csv`, `tapestation_ladders.csv`,
`hsd5000_low_mw.csv`, `hsd5000_ladders.csv`, `lower_marker_comparison.csv`,
`nanodrop_qubit.csv`, `nanodrop_spectra_normalized.csv`, `summary_stats.csv`
(each pub value, recomputed, with a match flag) and `group_summary.csv`
(n, mean, SD, min and max of every metric for every group).

## Notes

- Means are arithmetic and SDs are sample SDs (n − 1).
- All 11 Genomic DNA ScreenTape ladders resolve their 14 bands (10 at a 25%
  peak-height threshold, 1 at 20%); both HS D5000 ladders resolve their 11.
- Trace plots are drawn against acquisition position with size labels at their
  calibrated positions, because the size calibration isn't log-linear. Traces
  from different runs are resampled onto one reference run by calibrated size.

## License

MIT (see `LICENSE`).
