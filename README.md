# Code Release — Clozapine-Associated Pulmonary Infections in FAERS (2022–2025)

Analysis code for the manuscript:

**Pneumonia and Fatal Outcomes Associated with Clozapine: A Pharmacovigilance Study Using the FDA Adverse Event Reporting System** (manuscript version v11, 2026-08-13)

## Data Provenance

All analyses use **FAERS ASCII quarterly data files** (not the openFDA API).
The quarterly ASCII files for 2022Q1–2025Q4 are publicly available at:

https://www.fda.gov/drugs/questions-and-answers-fdas-adverse-event-reporting-system-faers/fda-adverse-event-reporting-system-faers-latest-quarterly-data-files

- `clozapine_risperidone_demo_2022_2025.txt`, `clozapine_risperidone_drug_2022_2025.txt`, `clozapine_risperidone_reac_2022_2025.txt` (referenced by the R scripts, located in `F:/faersdata/`): `$`-delimited extracts of the DEMO, DRUG, and REAC tables of the 16 quarterly ASCII files, restricted to records whose `primaryid` appears in a DRUG row matching the clozapine or risperidone drug-name patterns (case-insensitive substring match on `drugname`/`prod_ai`). Columns follow the standard FAERS ASCII layouts (25-field DEMO, 21-field DRUG, 4-field REAC).
- `faers_ascii_2025q4/ASCII/OUTC*.txt`: the 16 quarterly OUTC files (5,045,552 rows), used to define fatal outcome.
- `data/olanzapine_ascii_extraction.json`: olanzapine primary-suspect counts extracted from the same ASCII quarters (see `python/extract_olanzapine_ascii.py`).
- `data/clozapine_lung_signal_detection.csv`: authoritative per-PT signal-detection output (Table 2 numerators and the clozapine side of Table S1).

## Fatal-Outcome Definition (IMPORTANT)

Fatal outcome is defined as **`outc_cod` = "DE" (death) in the OUTC table**
(manuscript Section 2.5). An earlier iteration of this analysis used
`i_f_code == "F"` from the DEMO table, which is the **initial/follow-up report
flag**, not a death indicator. All scripts in this release use the corrected
OUTC definition: 309/1,305 cases (23.7%) fatal.

## File Inventory

### R scripts (run in numeric order)

| File | Purpose | Key outputs |
|---|---|---|
| `R/01_signal_detection.R` | Clozapine vs risperidone disproportionality (ROR/PRR/IC) across pulmonary-infection PTs; writes `clozapine_lung_signal_detection.csv` | Table 2 numerators; Figure 2 data |
| `R/02_fatal_analysis.R` | Corrected fatal-outcome analysis (OUTC DE): overall fatality, age/sex/weight comparisons, per-PT and country fatality, multivariable model m4 (C-statistic, Hosmer-Lemeshow), dose subset m3, DE-or-LT sensitivity | Table 3; Figures 4–5 data |
| `R/03_table1_counts.R` | Table 1 cell values (demographics × fatal status) | Table 1 |
| `R/04_dose_response.R` | Dose standardization (mg/d), RCS dose–aspiration modeling, threshold tests, age/sex subgroup analyses | Figure 3; Supplementary Figure S1 data |

Required R packages: `data.table`, `dplyr`, `pROC`, `ResourceSelection` (versions in manuscript Section 2.9). R 4.6.0.

### Python scripts

| File | Purpose |
|---|---|
| `python/extract_olanzapine_ascii.py` | Extract olanzapine primary-suspect counts from FAERS ASCII → `data/olanzapine_ascii_extraction.json` |
| `python/recompute_disproportionality_ascii.py` | Clozapine vs olanzapine disproportionality (Table S1) from ASCII counts (clozapine side read from `data/clozapine_lung_signal_detection.csv`, olanzapine side from `data/olanzapine_ascii_extraction.json`) |
| `python/generate_figures.py` | Publication figures 1, 2, 4, 5 (matplotlib) |

## File Paths

The R scripts read inputs from `F:/faersdata/`; adjust `data_dir`/`qdir`
variables at the top of each script for other machines. The Python scripts
write to `F:/clazpin/manuscript/analysis/figures/` (figures) and
`F:/clo-ola/data/` (olanzapine outputs).

## Manuscript Cross-Reference

- Denominators: clozapine PS 45,749 DRUG records → 44,055 unique reports (1,694 duplicate drug records across 985 reports); risperidone 15,130; olanzapine 12,701; clozapine pulmonary-infection cases 1,305.
- Fatal analysis: 309/1,305 (23.7%); multivariable model m4: age OR 1.90 (95% CI 1.64–2.19), aspiration pneumonia OR 2.40 (1.58–3.66), male sex 0.99 (0.67–1.45), Europe 0.75 (0.52–1.10), year 0.97 (0.90–1.04); C-statistic 0.774 (0.735–0.812); Hosmer-Lemeshow χ² = 11.38, df = 8, p = 0.181.
- Dose–fatal: OR per 100 mg/d = 0.96 (0.84–1.11); dose subset fatality 89/405 (22.0%) vs no-dose 220/900 (24.4%).
