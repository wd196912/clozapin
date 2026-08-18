# Analysis Code — Clozapine-Associated Pulmonary Infections in FAERS (2022–2025)

R analysis code for the manuscript:

**Pneumonia and Fatal Outcomes Associated with Clozapine: A Pharmacovigilance Study Using the FDA Adverse Event Reporting System** (manuscript version v12, composite-endpoint revision, 2026-08-18)

## Repository Contents

This repository contains only the R analysis scripts (contingency tables, ROR/PRR/IC disproportionality, logistic regression, restricted cubic spline dose–response modeling) and this README. **No FAERS data are redistributed here.**

## Data Source

All analyses use the **FAERS ASCII quarterly data files** (not the openFDA API), which are publicly downloadable from the FDA:

https://www.fda.gov/drugs/questions-and-answers-fdas-adverse-event-reporting-system-faers/fda-adverse-event-reporting-system-faers-latest-quarterly-data-files

The R scripts read the following inputs (referenced by the `data_dir` / `qdir` variables at the top of each script):

- `clozapine_risperidone_demo_2022_2025.txt`, `clozapine_risperidone_drug_2022_2025.txt`, `clozapine_risperidone_reac_2022_2025.txt`: `$`-delimited extracts of the DEMO, DRUG, and REAC tables of the 16 quarterly ASCII files (2022Q1–2025Q4), restricted to records whose `primaryid` appears in a DRUG row matching the clozapine or risperidone drug-name patterns (case-insensitive substring match on `drugname`/`prod_ai`). Columns follow the standard FAERS ASCII layouts (25-field DEMO, 21-field DRUG, 4-field REAC).
- `faers_ascii_2025q4/ASCII/OUTC*.txt`: the 16 quarterly OUTC files, used to define fatal outcome.

## Fatal-Outcome Definition (IMPORTANT)

Fatal outcome is defined as **`outc_cod` = "DE" (death) in the OUTC table** (manuscript Section 2.5). An earlier iteration used `i_f_code == "F"` from the DEMO table, which is the initial/follow-up report flag, not a death indicator. All scripts use the corrected OUTC definition: 309/1,305 cases (23.7%) fatal.

## Composite Endpoint and Revised Positivity Rule (v12)

Screening uses **26 MedDRA Preferred Terms (PTs)** — the 17 protocol PTs plus 9 additional PTs identified in the clozapine cohort (pneumonia necrotising, pneumonia anthrax, pulmonary sepsis, acute interstitial pneumonitis, atypical pneumonia, eosinophilic pneumonia, hypersensitivity pneumonitis, organising pneumonia, pleural infection) — reproducing N = 1,305 exactly.

- **Primary analysis — composite pulmonary infection endpoint (10 PTs):** pneumonia, pneumonia aspiration, lower respiratory tract infection, pneumonia bacterial, pneumonia viral, empyema, lung abscess, pneumonia klebsiella, pneumonia influenzal, pneumonia staphylococcal.
- **Clinical subgroups:** (1) pneumonia spectrum (composite minus lower respiratory tract infection; 9 PTs); (2) non-pneumonia lower respiratory tract infection (1 PT); (3) special/other pulmonary infection terms (the remaining 16 screening PTs). Individual PT analyses are exploratory.
- **Revised positivity rule:** a positive signal requires **IC₀₂₅ > 0 with ≥5 events in the comparator group**; ROR/PRR are not reported when the comparator group has fewer than 5 events (IC reported alone).
- **Principal-diagnosis priority chain** assigns each case to a single infection type for Table 1 (aspiration > pneumonia incl. etiologic subtypes > lower respiratory tract infection > COVID-19 pneumonia > upper respiratory tract infection > respiratory tract infection > pneumonitis > empyema > lung abscess > pulmonary tuberculosis > pulmonary sepsis > idiopathic interstitial pneumonia); counts sum exactly to 1,305.

These design elements were developed in response to expert review after preliminary per-term analyses revealed instability from small comparator denominators; the manuscript declares them as post-hoc (§4.3).

## R Scripts (run in numeric order)

| File | Purpose | Key outputs |
|---|---|---|
| `R/01_signal_detection.R` | Clozapine vs risperidone disproportionality (ROR/PRR/IC) across pulmonary-infection PTs | `clozapine_lung_signal_detection.csv` (Table 2 numerators; Figure 2 data) |
| `R/02_fatal_analysis.R` | Fatal-outcome analysis (OUTC DE): overall fatality, age/sex/weight comparisons, per-PT and country fatality, multivariable model m4 (C-statistic, Hosmer-Lemeshow), dose subset m3, DE-or-LT sensitivity | Table 3; Figures 4–5 data |
| `R/03_table1_counts.R` | Table 1 cell values (demographics × fatal status) | Table 1 |
| `R/04_dose_response.R` | Dose standardization (mg/d), RCS dose–aspiration modeling, threshold tests, age/sex subgroup analyses | Figure 3; Supplementary Figure S1 data |
| `R/05_sensitivity_analyses.R` | Sensitivity analyses of the fatal-outcome model m4: exclusion of imputed years, median-year imputation, categorical year, no-age and age missing-indicator models, multiple imputation by chained equations (mice, m = 10), caseid-level analysis; calibration plot | Supplementary Table S3; Supplementary Figure S2 |
| `R/06_composite_endpoint.R` | v12 revision: composite endpoint + 3 subgroups + 17 individual PTs vs risperidone under the revised positivity rule, principal-diagnosis assignment (priority chain), partition verification (1,305/309) and ROR back-solving | `composite_endpoint_signal_results.csv`, `principal_diagnosis_table1.csv` (Table 1 infection rows; Table 2; Figure 2) |

`R/06_composite_endpoint.R` includes a validation gate that reproduces the per-PT counts from `clozapine_lung_signal_detection.csv`; if that file is absent, the gate is skipped with a notice (run `01_signal_detection.R` first to generate it).

Required R packages: `data.table`, `dplyr`, `pROC`, `ResourceSelection`, `mice` (versions in manuscript Section 2.9). R 4.6.0.

## File Paths

Adjust `data_dir` / `qdir` (and any `out_dir` / `output_dir`) at the top of each script for your machine.

## Key Results Cross-Reference

- Denominators: clozapine PS 45,749 DRUG records → 44,055 unique reports; risperidone 15,130; olanzapine 13,691 (report level; from the same ASCII quarters); clozapine pulmonary-infection cases 1,305; risperidone 119.
- Composite endpoint vs risperidone: 1,177 vs 115 events; ROR 3.58 (95% CI 2.96–4.34); PRR 3.51 (2.91–4.25); IC₀₂₅ 0.168 — **Positive**. Pneumonia spectrum: 963 vs 114; ROR 2.94 (2.42–3.58) — **Positive**. Pneumonia: 723 vs 75; ROR 3.35 (2.64–4.25) — **Positive**. Pneumonia aspiration: 209 vs 40; ROR 1.80 (1.28–2.52); IC₀₂₅ −0.112 — not positive. Lower respiratory tract infection: 237 vs 1 — ROR not reported (IC 0.419, IC₀₂₅ 0.138); the raw per-term estimate (ROR 81.83) is retained only as an example of small-denominator instability.
- Fatal analysis: 309/1,305 (23.7%); multivariable model m4: age OR 1.90 (95% CI 1.64–2.19), aspiration pneumonia OR 2.40 (1.58–3.66); C-statistic 0.774 (0.735–0.812); Hosmer-Lemeshow χ² = 11.38, df = 8, p = 0.181.
