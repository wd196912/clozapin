# Code Release — Clozapine-Associated Pulmonary Infections in FAERS (2022–2025)

Analysis code for the manuscript:

**Pneumonia and Fatal Outcomes Associated with Clozapine: A Pharmacovigilance Study Using the FDA Adverse Event Reporting System** (manuscript version v12, composite-endpoint revision, 2026-08-18)

## Data Provenance

All analyses use **FAERS ASCII quarterly data files** (not the openFDA API).
The quarterly ASCII files for 2022Q1–2025Q4 are publicly available at:

https://www.fda.gov/drugs/questions-and-answers-fdas-adverse-event-reporting-system-faers/fda-adverse-event-reporting-system-faers-latest-quarterly-data-files

- `clozapine_risperidone_demo_2022_2025.txt`, `clozapine_risperidone_drug_2022_2025.txt`, `clozapine_risperidone_reac_2022_2025.txt` (referenced by the R scripts, located in `F:/faersdata/`): `$`-delimited extracts of the DEMO, DRUG, and REAC tables of the 16 quarterly ASCII files, restricted to records whose `primaryid` appears in a DRUG row matching the clozapine or risperidone drug-name patterns (case-insensitive substring match on `drugname`/`prod_ai`). Columns follow the standard FAERS ASCII layouts (25-field DEMO, 21-field DRUG, 4-field REAC).
- `faers_ascii_2025q4/ASCII/OUTC*.txt`: the 16 quarterly OUTC files (5,045,552 rows), used to define fatal outcome.
- `data/olanzapine_ascii_extraction.json`: olanzapine primary-suspect counts extracted from the same ASCII quarters (see `python/extract_olanzapine_ascii.py`); contains both the caseid-level (12,701) and primaryid/report-level (13,691) denominators.
- `data/clozapine_lung_signal_detection.csv`: authoritative per-PT signal-detection output of the original per-term analysis (17 PTs, old positivity criteria) — reproduced by `R/06_composite_endpoint.R` as a validation gate.
- `data/composite_endpoint_signal_results.csv`: authoritative output of the revised (composite-endpoint) signal-detection analysis — composite endpoint, three clinical subgroups, and 17 individual PTs vs risperidone, under the revised positivity rule.
- `data/principal_diagnosis_table1.csv`: principal-diagnosis assignment of the 1,305 cases (Table 1 infection-type rows).
- `data/disproportionality_clo_vs_ola_ascii.json`: clozapine vs olanzapine disproportionality (Table S1 source; report-level denominator 13,691).

## Fatal-Outcome Definition (IMPORTANT)

Fatal outcome is defined as **`outc_cod` = "DE" (death) in the OUTC table**
(manuscript Section 2.5). An earlier iteration of this analysis used
`i_f_code == "F"` from the DEMO table, which is the **initial/follow-up report
flag**, not a death indicator. All scripts in this release use the corrected
OUTC definition: 309/1,305 cases (23.7%) fatal.

## Composite Endpoint and Revised Positivity Rule (v12)

Screening uses **26 MedDRA Preferred Terms (PTs)** — the 17 protocol PTs plus
9 additional PTs identified in the clozapine cohort (pneumonia necrotising,
pneumonia anthrax, pulmonary sepsis, acute interstitial pneumonitis, atypical
pneumonia, eosinophilic pneumonia, hypersensitivity pneumonitis, organising
pneumonia, pleural infection) — reproducing N = 1,305 exactly.

- **Primary analysis — composite pulmonary infection endpoint (10 PTs):**
  pneumonia, pneumonia aspiration, lower respiratory tract infection,
  pneumonia bacterial, pneumonia viral, empyema, lung abscess,
  pneumonia klebsiella, pneumonia influenzal, pneumonia staphylococcal.
- **Clinical subgroups:** (1) pneumonia spectrum (composite minus lower
  respiratory tract infection; 9 PTs); (2) non-pneumonia lower respiratory
  tract infection (1 PT); (3) special/other pulmonary infection terms
  (the remaining 16 screening PTs). Individual PT analyses are exploratory.
- **Revised positivity rule:** a positive signal requires **IC₀₂₅ > 0 with
  ≥5 events in the comparator group**; ROR/PRR are not reported when the
  comparator group has fewer than 5 events (IC reported alone).
- **Principal-diagnosis priority chain** assigns each case to a single
  infection type for Table 1 (aspiration > pneumonia incl. etiologic
  subtypes > lower respiratory tract infection > COVID-19 pneumonia > upper
  respiratory tract infection > respiratory tract infection > pneumonitis >
  empyema > lung abscess > pulmonary tuberculosis > pulmonary sepsis >
  idiopathic interstitial pneumonia); counts sum exactly to 1,305.
- **Supplementary Table S4** compares three classification schemes
  (A: original per-PT, old criteria; B: composite endpoint; C: subgroups).

These design elements were developed in response to expert review after
preliminary per-term analyses revealed instability from small comparator
denominators; the manuscript declares them as post-hoc (§4.3).

## File Inventory

### R scripts (run in numeric order)

| File | Purpose | Key outputs |
|---|---|---|
| `R/01_signal_detection.R` | Clozapine vs risperidone disproportionality (ROR/PRR/IC) across pulmonary-infection PTs; writes `clozapine_lung_signal_detection.csv` | Table 2 numerators; Figure 2 data |
| `R/02_fatal_analysis.R` | Corrected fatal-outcome analysis (OUTC DE): overall fatality, age/sex/weight comparisons, per-PT and country fatality, multivariable model m4 (C-statistic, Hosmer-Lemeshow), dose subset m3, DE-or-LT sensitivity | Table 3; Figures 4–5 data |
| `R/03_table1_counts.R` | Table 1 cell values (demographics × fatal status) | Table 1 |
| `R/04_dose_response.R` | Dose standardization (mg/d), RCS dose–aspiration modeling, threshold tests, age/sex subgroup analyses | Figure 3; Supplementary Figure S1 data |
| `R/05_sensitivity_analyses.R` | Sensitivity analyses of the fatal-outcome model m4: exclusion of imputed years, median-year imputation, categorical year, no-age and age missing-indicator models, multiple imputation by chained equations (mice, m = 10) for age and sex, caseid-level analysis, age structure by region, crude-vs-adjusted Europe contrast; calibration plot | Supplementary Table S3; Supplementary Figure S2 |
| `R/06_composite_endpoint.R` | v12 revision: validates all original per-PT counts (gate), computes composite endpoint + 3 subgroups + 17 PTs vs risperidone under the revised positivity rule, assigns principal diagnoses (priority chain), verifies partition sums (1,305/309) and ROR back-solving; writes `composite_endpoint_signal_results.csv` and `principal_diagnosis_table1.csv` | Table 1 infection rows; Table 2; Figure 2; Scheme B/C of Table S4 |

Required R packages: `data.table`, `dplyr`, `pROC`, `ResourceSelection`, `mice` (3.19.0) (versions in manuscript Section 2.9). R 4.6.0.

### Python scripts

| File | Purpose |
|---|---|
| `python/extract_olanzapine_ascii.py` | Extract olanzapine primary-suspect counts for all 26 pulmonary-infection PTs from FAERS ASCII at both caseid and primaryid level → `data/olanzapine_ascii_extraction.json` |
| `python/recompute_disproportionality_ascii.py` | Clozapine vs olanzapine disproportionality (Table S1) under the revised rule; clozapine counts read from `data/composite_endpoint_signal_results.csv`; writes `data/disproportionality_clo_vs_ola_ascii.json` |
| `python/make_table_s1.py` | Generate Supplementary Table S1 docx (report-level olanzapine denominator 13,691) |
| `python/make_tables_1_2.py` | Rebuild Tables 1 and 2 of `tables_1_3.docx` (principal-diagnosis rows; composite endpoint + NR rule) |
| `python/make_table_s4.py` | Generate Supplementary Table S4 docx (three-scheme classification sensitivity) |
| `python/qc_round2_check.py` | QC: re-derives every Table 1 percentage, Table 2 and S1 ROR/CI from the 2 × 2 counts, and checks the rendered manuscript text |
| `python/generate_figures.py` | Publication figures 1, 2, 4, 5 (matplotlib); Figure 2 shows the composite endpoint with the revised positivity rule |

## File Paths

The R scripts read inputs from `F:/faersdata/`; adjust `data_dir`/`qdir`
variables at the top of each script for other machines. The Python scripts
write to `F:/clazpin/manuscript/analysis/figures/` (figures) and
`F:/clo-ola/data/` (olanzapine outputs).

## Manuscript Cross-Reference

- Denominators: clozapine PS 45,749 DRUG records → 44,055 unique reports (1,694 duplicate drug records across 985 reports); risperidone 15,130; olanzapine 13,691 (report level); clozapine pulmonary-infection cases 1,305; risperidone pulmonary-infection cases 119.
- Composite endpoint vs risperidone: 1,177 vs 115 events; ROR 3.58 (95% CI 2.96–4.34); PRR 3.51 (2.91–4.25); IC₀₂₅ 0.168 — **Positive**. Pneumonia spectrum subgroup: 963 vs 114; ROR 2.94 (2.42–3.58); IC₀₂₅ 0.129 — **Positive**. Pneumonia: 723 vs 75; ROR 3.35 (2.64–4.25); IC₀₂₅ 0.127 — **Positive**. Pneumonia aspiration: 209 vs 40; ROR 1.80 (1.28–2.52); IC₀₂₅ −0.112 — not positive. Lower respiratory tract infection: 237 vs 1 — ROR not reported (IC 0.419, IC₀₂₅ 0.138); the raw Scheme-A estimate (ROR 81.83) is retained only in Table S4 and the discussion as an example of small-denominator instability.
- Vs olanzapine (Supplementary Table S1): composite 1177 vs 360; ROR 1.02 (0.90–1.15), not positive; pneumonia 723 vs 197; ROR 1.14 (0.98–1.34); pneumonia aspiration 209 vs 140; ROR 0.46 (0.37–0.57); lower respiratory tract infection 237 vs 19; ROR 3.89 (2.44–6.21), IC₀₂₅ 0.006 — **Positive**.
- Fatal analysis: 309/1,305 (23.7%); multivariable model m4: age OR 1.90 (95% CI 1.64–2.19), aspiration pneumonia OR 2.40 (1.58–3.66), male sex 0.99 (0.67–1.45), Europe 0.75 (0.52–1.10), year 0.97 (0.90–1.04); C-statistic 0.774 (0.735–0.812); Hosmer-Lemeshow χ² = 11.38, df = 8, p = 0.181.
- Dose–fatal: OR per 100 mg/d = 0.96 (0.84–1.11); dose subset fatality 89/405 (22.0%) vs no-dose 220/900 (24.4%).
- Sensitivity analyses (Supplementary Table S3): complete-case m4 on 797 cases (184 fatal events); excluding imputed years 566 (138); median-year imputation 797 (184); MICE m = 10 pooled n = 1,305 (309); model without age 1,260 (303); age missing-indicator 1,260 (303); caseid-level 707 (171). All estimates consistent with the primary analysis. Caseid-level fatality 294/1,162 (25.3%) vs primaryid-level 23.7%.
