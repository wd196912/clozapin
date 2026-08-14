# Compute all Table 1 cells with corrected fatal definition (OUTC DE)
# Outputs CSV for docx patching: F:/clazpin/table1_cells.csv
library(data.table)
library(dplyr)

data_dir <- "F:/faersdata"
qdir <- "F:/faersdata/faers_ascii_2025q4/ASCII"

demo_cols <- c("primaryid", "caseid", "caseversion", "i_f_code",
               "event_dt", "mfr_dt", "init_fda_dt", "fda_dt",
               "rept_cod", "auth_num", "mfr_num", "mfr_sndr",
               "lit_ref", "age", "age_cod", "age_grp", "sex",
               "e_sub", "wt", "wt_cod", "rept_dt", "to_mfr",
               "occp_cod", "reporter_country", "occr_country")
drug_cols <- c("primaryid", "caseid", "drug_seq", "role_cod",
               "drugname", "prod_ai", "val_vbm", "route",
               "dose_vbm", "cum_dose_chr", "cum_dose_unit",
               "dechal", "rechal", "lot_num", "exp_dt",
               "nda_num", "dose_amt", "dose_unit", "dose_form",
               "dose_freq")
reac_cols <- c("primaryid", "caseid", "pt", "drug_rec_act")

demo_raw <- fread(file.path(data_dir, "clozapine_risperidone_demo_2022_2025.txt"),
                  sep = "$", header = FALSE, col.names = demo_cols,
                  na.strings = "", blank.lines.skip = TRUE, strip.white = TRUE)
drug_raw <- fread(file.path(data_dir, "clozapine_risperidone_drug_2022_2025.txt"),
                  sep = "$", header = FALSE, col.names = drug_cols,
                  na.strings = "", blank.lines.skip = TRUE, strip.white = TRUE)
reac_raw <- fread(file.path(data_dir, "clozapine_risperidone_reac_2022_2025.txt"),
                  sep = "$", header = FALSE, col.names = reac_cols,
                  na.strings = "", blank.lines.skip = TRUE, strip.white = TRUE)

outc_files <- list.files(qdir, pattern = "^OUTC.*\\.txt$", full.names = TRUE)
outc <- rbindlist(lapply(outc_files, function(f) {
  fread(f, sep = "$", header = FALSE, select = 1:3,
        col.names = c("primaryid", "caseid", "outc_cod"),
        na.strings = "", blank.lines.skip = TRUE, strip.white = TRUE)
}))
death_ids <- unique(outc[outc_cod == "DE", primaryid])

drug_raw[, drug_type := ifelse(
  grepl("CLOZAPINE|CLOZARIL|FAZACLO|VERSACLOZ|LEPONEX|ZAPONEX|CLOPINE|DENZAPINE",
        drugname, ignore.case = TRUE) &
    !grepl("RISPERIDONE|RISPERDAL|PERSERIS", drugname, ignore.case = TRUE),
  "Clozapine", "Other")]

lung_terms <- c(
  "Pneumonia$", "Pneumonia aspiration$", "Pneumonitis$",
  "Bronchopneumonia$", "Lung infection$", "Pulmonary infection$",
  "Pulmonary sepsis$", "Lower respiratory tract infection$",
  "Respiratory tract infection$", "Atypical pneumonia$",
  "Bacterial pneumonia$", "Viral pneumonia$", "Fungal pneumonia$",
  "Pneumonia bacterial$", "Pneumonia viral$", "Pneumonia fungal$",
  "Pneumonia staphylococcal$", "Pneumonia streptococcal$",
  "Pneumonia pneumococcal$", "Pneumonia klebsiella$",
  "Pneumonia pseudomonas$", "Pneumonia mycoplasmal$",
  "Pneumonia influenzal$", "Pneumonia necrotising$",
  "Pneumocystis jirovecii pneumonia$", "Lung abscess$",
  "Empyema$", "Pleural infection$", "Pulmonary tuberculosis$",
  "COVID-19 pneumonia$", "Coronavirus pneumonia$",
  "Tracheobronchitis$", "Organising pneumonia$",
  "Eosinophilic pneumonia$", "Pneumonia anthrax$",
  "Radiation pneumonitis$", "Pulmonary mycosis$",
  "Pneumonia cytomegaloviral$",
  "Pneumonia respiratory syncytial viral$",
  "Pneumonia adenoviral$", "Pneumonia herpes viral$",
  "Hypersensitivity pneumonitis$",
  "Idiopathic interstitial pneumonia$",
  "Acute interstitial pneumonitis$",
  "Interstitial pneumonitis$", "Pneumonia chemical$",
  "Pneumonia lipid$", "Metapneumovirus pneumonia$",
  "Upper respiratory tract infection$",
  "Lower respiratory tract infection viral$",
  "Respiratory tract infection viral$",
  "Respiratory tract infection bacterial$",
  "Viral upper respiratory tract infection$",
  "Respiratory moniliasis$",
  "Infective exacerbation of bronchiectasis$"
)
lung_pattern <- paste(lung_terms, collapse = "|")
reac_raw[, is_lung := grepl(lung_pattern, pt, ignore.case = TRUE)]

drug_cloz_ps <- drug_raw[drug_type == "Clozapine" & role_cod == "PS"]

cloz_lung <- reac_raw[is_lung == TRUE] %>%
  inner_join(drug_cloz_ps, by = c("primaryid", "caseid"),
             relationship = "many-to-many") %>%
  inner_join(demo_raw, by = c("primaryid", "caseid"),
             relationship = "many-to-many") %>%
  as.data.table()

cloz_lung[, age_num := as.numeric(age)]
cloz_lung[age_cod %in% c("DEC"), age_num := age_num / 12]
cloz_lung[age_cod %in% c("WK","WEEK"), age_num := age_num / 52]
cloz_lung[age_cod %in% c("DY","DAY"), age_num := age_num / 365]
cloz_lung[, event_year := as.integer(substr(event_dt, 1, 4))]
cloz_lung[, event_month := as.integer(substr(event_dt, 5, 6))]
cloz_lung[is.na(event_month) | event_month < 1 | event_month > 12, event_month := 6L]
cloz_lung[is.na(event_year) | event_year < 1900, c("event_year","event_month") := .(2022L, 6L)]
cloz_lung[, wt_num := as.numeric(wt)]
cloz_lung[wt_cod %in% c("LBS","LB","POUNDS","POUND"), wt_kg := wt_num * 0.4536]
cloz_lung[wt_cod %in% c("KG","KGS","KILOGRAMS","KILOGRAM"), wt_kg := wt_num]
cloz_lung[is.na(wt_kg), wt_kg := wt_num]

cloz_lung[, is_fatal := primaryid %in% death_ids]
cloz_lung[, is_aspiration := grepl("aspiration|aspiratio", pt, ignore.case = TRUE)]
cloz_lung[, caseversion_num := as.numeric(caseversion)]
cloz_by_id <- cloz_lung[order(primaryid, -is_fatal, -caseversion_num)][
  , .SD[1], by = primaryid]

N_ALL <- nrow(cloz_by_id)
N_F <- sum(cloz_by_id$is_fatal)
N_NF <- N_ALL - N_F
cat(sprintf("All %d, fatal %d, non-fatal %d\n", N_ALL, N_F, N_NF))

# ---- Overall (All) column ----
age_all <- cloz_by_id$age_num
wt_all <- cloz_by_id$wt_kg
cat(sprintf("AGE_ALL median %.0f (IQR %.0f-%.0f), n=%d\n",
            median(age_all, na.rm=TRUE), quantile(age_all, 0.25, na.rm=TRUE),
            quantile(age_all, 0.75, na.rm=TRUE), sum(!is.na(age_all))))
cat(sprintf("AGE_F median %.0f (IQR %.0f-%.0f), n=%d\n",
            median(cloz_by_id[is_fatal==TRUE, age_num], na.rm=TRUE),
            quantile(cloz_by_id[is_fatal==TRUE, age_num], 0.25, na.rm=TRUE),
            quantile(cloz_by_id[is_fatal==TRUE, age_num], 0.75, na.rm=TRUE),
            sum(!is.na(cloz_by_id[is_fatal==TRUE, age_num]))))
cat(sprintf("AGE_NF median %.0f (IQR %.0f-%.0f), n=%d\n",
            median(cloz_by_id[is_fatal==FALSE, age_num], na.rm=TRUE),
            quantile(cloz_by_id[is_fatal==FALSE, age_num], 0.25, na.rm=TRUE),
            quantile(cloz_by_id[is_fatal==FALSE, age_num], 0.75, na.rm=TRUE),
            sum(!is.na(cloz_by_id[is_fatal==FALSE, age_num]))))
cat(sprintf("AGE Wilcoxon p = %.4g\n",
            wilcox.test(cloz_by_id[is_fatal==TRUE, age_num],
                        cloz_by_id[is_fatal==FALSE, age_num])$p.value))
cat(sprintf("WT_ALL median %.0f (IQR %.0f-%.0f), n=%d\n",
            median(wt_all, na.rm=TRUE), quantile(wt_all, 0.25, na.rm=TRUE),
            quantile(wt_all, 0.75, na.rm=TRUE), sum(!is.na(wt_all))))

# ---- Age strata ----
age_cuts <- c(0, 30, 40, 50, 60, 70, 80, Inf)
strata_names <- c("<30", "30-39", "40-49", "50-59", "60-69", "70-79", ">=80")
for (i in seq_along(strata_names)) {
  lo <- age_cuts[i]; hi <- age_cuts[i+1]
  a <- sum(cloz_by_id$age_num >= lo & cloz_by_id$age_num < hi, na.rm=TRUE)
  f <- sum(cloz_by_id$is_fatal & cloz_by_id$age_num >= lo & cloz_by_id$age_num < hi, na.rm=TRUE)
  nf <- sum(!cloz_by_id$is_fatal & cloz_by_id$age_num >= lo & cloz_by_id$age_num < hi, na.rm=TRUE)
  cat(sprintf("STRATA %-6s all=%d fatal=%d (%.1f%%) nonfatal=%d (%.1f%%)\n",
              strata_names[i], a, f, 100*f/N_F, nf, 100*nf/N_NF))
}

# ---- Sex ----
sex_ok <- cloz_by_id[sex %in% c("M","F")]
m_f <- sum(sex_ok$sex=="M" & sex_ok$is_fatal); m_nf <- sum(sex_ok$sex=="M" & !sex_ok$is_fatal)
f_f <- sum(sex_ok$sex=="F" & sex_ok$is_fatal); f_nf <- sum(sex_ok$sex=="F" & !sex_ok$is_fatal)
cat(sprintf("SEX male fatal %d (%.1f%%) nonfatal %d (%.1f%%)\n",
            m_f, 100*m_f/N_F, m_nf, 100*m_nf/N_NF))
cat(sprintf("SEX female fatal %d (%.1f%%) nonfatal %d (%.1f%%)\n",
            f_f, 100*f_f/N_F, f_nf, 100*f_nf/N_NF))
cat(sprintf("SEX unknown total %d (fatal %d, nonfatal %d)\n",
            sum(is.na(cloz_by_id$sex)), sum(is.na(cloz_by_id$sex) & cloz_by_id$is_fatal),
            sum(is.na(cloz_by_id$sex) & !cloz_by_id$is_fatal)))
sex_2x2 <- table(sex_ok$sex, sex_ok$is_fatal)
cat(sprintf("SEX Fisher p = %.3f; chisq p = %.3f\n",
            fisher.test(sex_2x2)$p.value, chisq.test(sex_2x2)$p.value))

# ---- Region ----
cloz_by_id[, region := "Other"]
cloz_by_id[reporter_country %in% c("GB","FR","DE","IT","ES","NL","BE","CH","AT","SE",
                                    "NO","DK","FI","IE","PT","GR","PL","CZ","RO","HU",
                                    "BG","HR","SK","LT","LV","EE","LU","MT","CY"), region := "Europe"]
cloz_by_id[reporter_country %in% c("US","CA"), region := "North America"]
cloz_by_id[reporter_country %in% c("AU","NZ"), region := "Oceania"]
cloz_by_id[reporter_country %in% c("JP","CN","KR","TW","HK","IN","TH","MY","SG",
                                    "ID","PH","VN","PK","BD"), region := "Asia"]
cloz_by_id[reporter_country %in% c("BR","AR","CL","CO","MX","PE","VE","UY","EC",
                                    "BO","PY"), region := "Latin America"]
for (r in c("Europe","North America","Oceania","Asia","Latin America","Other")) {
  tot <- sum(cloz_by_id$region == r)
  ff <- sum(cloz_by_id$region == r & cloz_by_id$is_fatal)
  nf <- sum(cloz_by_id$region == r & !cloz_by_id$is_fatal)
  cat(sprintf("REGION %-14s all=%d fatal=%d (%.1f%%) nonfatal=%d (%.1f%%)\n",
              r, tot, ff, 100*ff/N_F, nf, 100*nf/N_NF))
}

# ---- Infection PTs (one row per primaryid after dedup) ----
for (p in c("Pneumonia", "Lower respiratory tract infection",
            "Pneumonia aspiration", "Respiratory tract infection",
            "COVID-19 pneumonia")) {
  s <- cloz_by_id[pt == p]
  ff <- sum(s$is_fatal); nf <- nrow(s) - ff
  cat(sprintf("INF %-38s all=%d fatal=%d (%.1f%%) nonfatal=%d (%.1f%%)\n",
              p, nrow(s), ff, 100*ff/N_F, nf, 100*nf/N_NF))
}
