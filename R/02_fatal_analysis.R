# Recompute ALL fatal-outcome results with the correct definition:
# fatal = outc_cod == "DE" in OUTC (per manuscript section 2.5)
# Replaces the erroneous i_f_code == "F" (initial/follow-up flag).
library(data.table)
library(dplyr)
library(pROC)
library(ResourceSelection)

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

# ---- OUTC death ids ----
outc_files <- list.files(qdir, pattern = "^OUTC.*\\.txt$", full.names = TRUE)
outc <- rbindlist(lapply(outc_files, function(f) {
  fread(f, sep = "$", header = FALSE, select = 1:3,
        col.names = c("primaryid", "caseid", "outc_cod"),
        na.strings = "", blank.lines.skip = TRUE, strip.white = TRUE)
}))
death_ids <- unique(outc[outc_cod == "DE", primaryid])
cat(sprintf("OUTC death primaryids (all FAERS): %d\n", length(death_ids)))

# ---- drug classification ----
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

cat(sprintf("Cloz lung: %d rows, %d unique primaryid\n",
            nrow(cloz_lung), uniqueN(cloz_lung$primaryid)))

# ---- cleaning (exact replication) ----
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

# CORRECTED fatal definition
cloz_lung[, is_fatal := primaryid %in% death_ids]
cloz_lung[, is_aspiration := grepl("aspiration|aspiratio", pt, ignore.case = TRUE)]

cloz_lung[, caseversion_num := as.numeric(caseversion)]
cloz_by_id <- cloz_lung[order(primaryid, -is_fatal, -caseversion_num)][
  , .SD[1], by = primaryid]

cat(sprintf("\n=== OVERALL ===\n"))
cat(sprintf("Cases: %d, fatal: %d (%.1f%%), non-fatal: %d (%.1f%%)\n",
            nrow(cloz_by_id), sum(cloz_by_id$is_fatal),
            100*sum(cloz_by_id$is_fatal)/nrow(cloz_by_id),
            sum(!cloz_by_id$is_fatal),
            100*sum(!cloz_by_id$is_fatal)/nrow(cloz_by_id)))

cat(sprintf("\n=== AGE ===\n"))
f_age <- cloz_by_id[is_fatal == TRUE, age_num]
n_age <- cloz_by_id[is_fatal == FALSE, age_num]
cat(sprintf("Fatal: median %.0f (IQR %.0f-%.0f), n=%d\n",
            median(f_age, na.rm=TRUE), quantile(f_age, 0.25, na.rm=TRUE),
            quantile(f_age, 0.75, na.rm=TRUE), sum(!is.na(f_age))))
cat(sprintf("Non-fatal: median %.0f (IQR %.0f-%.0f), n=%d\n",
            median(n_age, na.rm=TRUE), quantile(n_age, 0.25, na.rm=TRUE),
            quantile(n_age, 0.75, na.rm=TRUE), sum(!is.na(n_age))))
cat(sprintf("Wilcoxon p = %.3f; t-test p = %.3f\n",
            wilcox.test(f_age, n_age)$p.value, t.test(f_age, n_age)$p.value))

cat(sprintf("\n=== SEX ===\n"))
sex_tbl <- cloz_by_id[, .N, by = sex]
print(sex_tbl)
m_f <- cloz_by_id[sex %in% c("M","F")]
sex_2x2 <- table(m_f$sex, m_f$is_fatal)
print(sex_2x2)
cat(sprintf("Fatal male: %d/%d (%.1f%%), non-fatal male: %d/%d (%.1f%%)\n",
            sum(m_f$sex=="M" & m_f$is_fatal), sum(m_f$is_fatal),
            100*sum(m_f$sex=="M" & m_f$is_fatal)/sum(m_f$is_fatal),
            sum(m_f$sex=="M" & !m_f$is_fatal), sum(!m_f$is_fatal),
            100*sum(m_f$sex=="M" & !m_f$is_fatal)/sum(!m_f$is_fatal)))
cat(sprintf("Chisq p = %.3f; Fisher p = %.3f\n",
            chisq.test(sex_2x2)$p.value, fisher.test(sex_2x2)$p.value))

cat(sprintf("\n=== WEIGHT ===\n"))
cat(sprintf("Weight available: %d (%.1f%%)\n",
            sum(!is.na(cloz_by_id$wt_kg)), 100*sum(!is.na(cloz_by_id$wt_kg))/nrow(cloz_by_id)))
cat(sprintf("Fatal median wt: %.0f (IQR %.0f-%.0f), n=%d\n",
            median(cloz_by_id[is_fatal==TRUE, wt_kg], na.rm=TRUE),
            quantile(cloz_by_id[is_fatal==TRUE, wt_kg], 0.25, na.rm=TRUE),
            quantile(cloz_by_id[is_fatal==TRUE, wt_kg], 0.75, na.rm=TRUE),
            sum(!is.na(cloz_by_id[is_fatal==TRUE, wt_kg]))))
cat(sprintf("Non-fatal median wt: %.0f (IQR %.0f-%.0f), n=%d\n",
            median(cloz_by_id[is_fatal==FALSE, wt_kg], na.rm=TRUE),
            quantile(cloz_by_id[is_fatal==FALSE, wt_kg], 0.25, na.rm=TRUE),
            quantile(cloz_by_id[is_fatal==FALSE, wt_kg], 0.75, na.rm=TRUE),
            sum(!is.na(cloz_by_id[is_fatal==FALSE, wt_kg]))))
wt_grp <- cloz_by_id[!is.na(wt_kg) & wt_kg >= 90 & wt_kg <= 99]
if (nrow(wt_grp) > 0) {
  k <- sum(wt_grp$is_fatal); n <- nrow(wt_grp)
  bt <- binom.test(k, n)
  cat(sprintf("90-99kg subgroup: %d/%d fatal (%.1f%%, CI %.1f-%.1f)\n",
              k, n, 100*k/n, 100*bt$conf.int[1], 100*bt$conf.int[2]))
}

cat(sprintf("\n=== PER-PT FATALITY (one PT per case) ===\n"))
for (this_pt in c("Pneumonia viral","Respiratory tract infection",
                  "Lower respiratory tract infection","Upper respiratory tract infection",
                  "Pneumonia aspiration","Pneumonia","Pulmonary tuberculosis")) {
  s <- cloz_by_id[pt == this_pt]
  cat(sprintf("%-42s %d/%d (%.1f%%)\n", this_pt, sum(s$is_fatal), nrow(s),
              100*sum(s$is_fatal)/nrow(s)))
}
cat("\nAll PTs (n>=5):\n")
pt_all <- cloz_by_id[, .(Total=.N, Fatal=sum(is_fatal)), by=pt][Total>=5][order(-Total)]
print(pt_all)

cat(sprintf("\n=== COUNTRY (n>=10) ===\n"))
country_fatal <- cloz_by_id[, .(Total=.N, Fatal=sum(is_fatal),
                                Pct=round(100*sum(is_fatal)/.N,1)),
                            by=reporter_country][Total>=10][order(-Pct)]
print(country_fatal)

cat(sprintf("\n=== REGION ===\n"))
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
region_fatal <- cloz_by_id[, .(Total=.N, Fatal=sum(is_fatal),
                               Pct=round(100*sum(is_fatal)/.N,1)), by=region][order(-Pct)]
print(region_fatal)
eur <- cloz_by_id[, .(Total=.N, Fatal=sum(is_fatal)), by=.(region=="Europe")]
print(eur)
eur_2x2 <- table(cloz_by_id$region == "Europe", cloz_by_id$is_fatal)
cat(sprintf("Europe vs other chisq p = %.4f\n", chisq.test(eur_2x2)$p.value))

cat(sprintf("\n=== YEARLY ===\n"))
yearly <- cloz_by_id[event_year >= 2015 & event_year <= 2025, .(
  Total=.N, Fatal=sum(is_fatal), Pct=round(100*sum(is_fatal)/.N,1)), by=event_year][order(event_year)]
print(yearly)
yr_glm <- glm(cbind(Fatal, Total-Fatal) ~ event_year, data=yearly, family=binomial)
cat(sprintf("Yearly trend OR per year = %.3f (95%% CI %.3f-%.3f), p = %.3f\n",
            exp(coef(yr_glm)["event_year"]),
            exp(confint.default(yr_glm)["event_year",1]),
            exp(confint.default(yr_glm)["event_year",2]),
            coef(summary(yr_glm))["event_year","Pr(>|z|)"]))

cat(sprintf("\n=== MULTIVARIABLE m4 ===\n"))
model_data <- cloz_by_id[, .(primaryid, is_fatal, age_num, sex, is_aspiration, region, event_year)]
model_data[, sex_bin := ifelse(sex == "M", 1, 0)]
model_data[, aspiration := as.integer(is_aspiration)]
model_data[, region_eur := ifelse(region == "Europe", 1, 0)]
m4 <- glm(is_fatal ~ age_num + sex_bin + aspiration + region_eur + event_year,
          data = model_data, family = binomial)
cc_idx <- as.integer(rownames(model.frame(m4)))
cat(sprintf("Complete cases: %d / %d (%.1f%%)\n", length(cc_idx), nrow(model_data),
            100*length(cc_idx)/nrow(model_data)))
cat(sprintf("Fatal among complete cases: %d (%.1f%%)\n",
            sum(model_data$is_fatal[cc_idx]), 100*sum(model_data$is_fatal[cc_idx])/length(cc_idx)))
or_df <- data.frame(
  Predictor = c("Age (per 10yr)", "Male sex", "Aspiration pneumonia",
                "Europe region", "Year (per year)"),
  OR = c(exp(coef(m4)["age_num"] * 10), exp(coef(m4)["sex_bin"]),
         exp(coef(m4)["aspiration"]), exp(coef(m4)["region_eur"]),
         exp(coef(m4)["event_year"])),
  Lower = c(exp(confint.default(m4)["age_num", 1] * 10),
            exp(confint.default(m4)["sex_bin", 1]),
            exp(confint.default(m4)["aspiration", 1]),
            exp(confint.default(m4)["region_eur", 1]),
            exp(confint.default(m4)["event_year", 1])),
  Upper = c(exp(confint.default(m4)["age_num", 2] * 10),
            exp(confint.default(m4)["sex_bin", 2]),
            exp(confint.default(m4)["aspiration", 2]),
            exp(confint.default(m4)["region_eur", 2]),
            exp(confint.default(m4)["event_year", 2]))
)
or_df$P <- coef(summary(m4))[c("age_num","sex_bin","aspiration","region_eur","event_year"), "Pr(>|z|)"]
print(or_df, row.names = FALSE, digits = 4)

pred_prob <- predict(m4, type = "response")
roc_obj <- roc(model_data$is_fatal[cc_idx], pred_prob)
cat(sprintf("C-statistic = %.3f (95%% CI %.3f-%.3f)\n",
            auc(roc_obj), ci(roc_obj)[1], ci(roc_obj)[3]))
hl <- hoslem.test(model_data$is_fatal[cc_idx], pred_prob, g = 10)
cat(sprintf("Hosmer-Lemeshow: Chi-sq = %.2f, df = %d, p = %.3f\n",
            hl$statistic, hl$parameter, hl$p.value))

# E-value for Europe OR
or_eur <- exp(coef(m4)["region_eur"])
ci_eur <- exp(confint.default(m4)["region_eur", ])
e_pt <- or_eur + sqrt(or_eur * (or_eur - 1))
e_ci <- ci_eur[1] + sqrt(ci_eur[1] * (ci_eur[1] - 1))
cat(sprintf("Europe OR = %.3f (%.3f-%.3f), E-value: point %.2f, CI-limit %.2f\n",
            or_eur, ci_eur[1], ci_eur[2], e_pt, e_ci))

cat(sprintf("\n=== DOSE SUBSET (fatal) ===\n"))
# dose cleaning (replicate dose script)
dose_raw <- copy(cloz_lung)
dose_raw[, dose_amt_num := as.numeric(dose_amt)]
dose_raw[, dose_unit_clean := toupper(gsub("[^A-Z]", "", dose_unit))]
dose_raw[dose_unit_clean %in% c("MG","MGS","MILLIGRAMS","MILLIGRAM"), dose_mg := dose_amt_num]
dose_raw[dose_unit_clean %in% c("G","GM","GRAMS","GRAM"), dose_mg := dose_amt_num * 1000]
dose_raw[dose_unit_clean %in% c("UG","MCG","MICROGRAMS"), dose_mg := dose_amt_num / 1000]
dose_raw[dose_unit_clean %in% c("ML","MILLILITERS","MILLILITER"), dose_mg := dose_amt_num * 100]
dose_raw[dose_unit_clean == "NG", dose_mg := dose_amt_num / 1e6]
dose_raw[dose_unit_clean %in% c("IU","UNITS","UNIT"), dose_mg := NA]
dose_raw[, dose_freq_clean := toupper(gsub("[^A-Z0-9]", "", dose_freq))]
dose_raw[, daily_mult := 1]
dose_raw[grepl("BID|BID$|TWICE.*DAILY|TWODAILY|2XD|2X.*DAILY|TWICEADAY|EVERY12H|Q12H|BIDLY",
               dose_freq_clean), daily_mult := 2]
dose_raw[grepl("TID|TID$|THREE.*DAILY|3XD|3X.*DAILY|THREETIMESDAILY|EVERY8H|Q8H|TIDLY",
               dose_freq_clean), daily_mult := 3]
dose_raw[grepl("QID|QID$|FOUR.*DAILY|4XD|4X.*DAILY|FOURTIMESDAILY|EVERY6H|Q6H|QIDLY",
               dose_freq_clean), daily_mult := 4]
dose_raw[grepl("QOD|EVERYOTHER|ALTERNATE|EVERY2D|EVERY48H|Q48H",
               dose_freq_clean), daily_mult := 0.5]
dose_raw[grepl("QW|WEEKLY|ONCEAWEEK|EVERYWEEK|Q7D|EVERY7D",
               dose_freq_clean), daily_mult := 1/7]
dose_raw[grepl("BIW|TWICEWEEKLY|TWICEAWEEK|2XWEEKLY",
               dose_freq_clean), daily_mult := 2/7]
dose_raw[grepl("TIW|THREETIMESWEEKLY|3XWEEKLY",
               dose_freq_clean), daily_mult := 3/7]
dose_raw[, daily_dose_mg := dose_mg * daily_mult]
dose_by_id <- dose_raw[!is.na(daily_dose_mg) & daily_dose_mg > 0 & daily_dose_mg <= 2000,
                        .(daily_dose = median(daily_dose_mg, na.rm = TRUE)),
                        by = primaryid]
demo_unique <- cloz_by_id[, .(primaryid, age_num, sex, is_fatal, is_aspiration,
                              event_year, reporter_country)]
dose_full <- dose_by_id %>%
  inner_join(demo_unique, by = "primaryid") %>%
  as.data.table()
dose_full[, sex_bin := ifelse(sex == "M", 1, ifelse(sex == "F", 0, NA))]
dose_full[, fatal := as.integer(is_fatal)]

cat(sprintf("Dose subset: %d cases (%.1f%% of 1305)\n", nrow(dose_full),
            100*nrow(dose_full)/1305))
cat(sprintf("Dose subset fatal: %d/%d (%.1f%%)\n",
            sum(dose_full$fatal), nrow(dose_full), 100*sum(dose_full$fatal)/nrow(dose_full)))
cat(sprintf("No-dose subset fatal: %d/%d (%.1f%%)\n",
            sum(cloz_by_id[!primaryid %in% dose_full$primaryid, is_fatal]),
            nrow(cloz_by_id) - nrow(dose_full),
            100*sum(cloz_by_id[!primaryid %in% dose_full$primaryid, is_fatal])/
              (nrow(cloz_by_id) - nrow(dose_full))))
no_dose_fatal <- cloz_by_id[!primaryid %in% dose_full$primaryid, is_fatal]
dose_fatal <- dose_full$fatal
cat(sprintf("Dose vs no-dose fatality chisq p = %.3f\n",
            chisq.test(table(c(rep("dose", length(dose_fatal)), rep("nodose", length(no_dose_fatal))),
                             c(dose_fatal, no_dose_fatal)))$p.value))

m3 <- glm(fatal ~ daily_dose + age_num + sex_bin, data = dose_full, family = binomial)
or_fatal_100 <- exp(coef(m3)["daily_dose"] * 100)
ci_fatal_100 <- exp(confint.default(m3)["daily_dose", ] * 100)
cat(sprintf("Dose-fatal adjusted OR per 100 mg/d = %.2f (95%% CI %.2f-%.2f), p = %.3f\n",
            or_fatal_100, ci_fatal_100[1], ci_fatal_100[2],
            coef(summary(m3))["daily_dose","Pr(>|z|)"]))

# sensitivity: DE or LT
lt_ids <- unique(outc[outc_cod == "LT", primaryid])
cat(sprintf("\n=== SENSITIVITY: DE or LT ===\n"))
cat(sprintf("DE or LT fatal: %d (%.1f%%)\n",
            sum(cloz_by_id$primaryid %in% c(death_ids, lt_ids)),
            100*sum(cloz_by_id$primaryid %in% c(death_ids, lt_ids))/nrow(cloz_by_id)))
