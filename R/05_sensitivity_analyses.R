# Sensitivity analyses for the clozapine FAERS manuscript (v11 revision round)
# - Year-imputation sensitivities (exclude imputed / categorical / median imputation)
# - Complete-case sensitivities (model without age; age missing-indicator)
# - Multiple imputation by chained equations (mice, m=10) for age and sex
# - Caseid-level (latest version per case) sensitivity
# - Calibration plot (Supplementary Figure S2)
# Fatal definition: OUTC outc_cod == "DE" (manuscript section 2.5)
library(data.table)
library(dplyr)
library(pROC)
library(ResourceSelection)
library(mice)

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
cloz_lung[, event_year_raw := as.integer(substr(event_dt, 1, 4))]
cloz_lung[, event_year := event_year_raw]
cloz_lung[, year_imputed := is.na(event_year_raw) | event_year_raw < 1900]
cloz_lung[year_imputed == TRUE, c("event_year") := 2022L]
cloz_lung[, wt_num := as.numeric(wt)]
cloz_lung[wt_cod %in% c("LBS","LB","POUNDS","POUND"), wt_kg := wt_num * 0.4536]
cloz_lung[wt_cod %in% c("KG","KGS","KILOGRAMS","KILOGRAM"), wt_kg := wt_num]
cloz_lung[is.na(wt_kg), wt_kg := wt_num]

cloz_lung[, is_fatal := primaryid %in% death_ids]
cloz_lung[, is_aspiration := grepl("aspiration|aspiratio", pt, ignore.case = TRUE)]
cloz_lung[, caseversion_num := as.numeric(caseversion)]
cloz_by_id <- cloz_lung[order(primaryid, -is_fatal, -caseversion_num)][
  , .SD[1], by = primaryid]

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

md <- cloz_by_id[, .(primaryid, caseid, is_fatal, age_num, sex, is_aspiration,
                     region, event_year, year_imputed)]
md[, sex_bin := ifelse(sex == "M", 1, ifelse(sex == "F", 0, NA))]
md[, aspiration := as.integer(is_aspiration)]
md[, region_eur := ifelse(region == "Europe", 1, 0)]

or_line <- function(fit, name) {
  cf <- coef(fit); ci <- confint.default(fit); pv <- coef(summary(fit))[, "Pr(>|z|)"]
  cat(sprintf("[%s] n=%d events=%d\n", name, length(fit$y), sum(fit$y)))
  for (v in c("age_num","sex_bin","aspiration","region_eur","event_year")) {
    if (v %in% names(cf)) {
      if (v == "age_num") {
        cat(sprintf("  %-12s OR(1yr) %.3f (%.3f-%.3f); OR(10yr) %.3f (%.3f-%.3f) p=%.3f\n",
                    v, exp(cf[v]), exp(ci[v,1]), exp(ci[v,2]),
                    exp(10*cf[v]), exp(10*ci[v,1]), exp(10*ci[v,2]), pv[v]))
      } else {
        cat(sprintf("  %-12s OR %.3f (%.3f-%.3f) p=%.3f\n", v,
                    exp(cf[v]), exp(ci[v,1]), exp(ci[v,2]), pv[v]))
      }
    }
  }
}

cat("=== S1a: m4 excluding imputed years ===\n")
m4a <- glm(is_fatal ~ age_num + sex_bin + aspiration + region_eur + event_year,
           data = md[year_imputed == FALSE], family = binomial)
or_line(m4a, "excl imputed year")

cat("\n=== S1b: m4 with year categorical (imputed category) ===\n")
md[, year_cat := factor(ifelse(year_imputed, "imputed", as.character(event_year)))]
m4b <- glm(is_fatal ~ age_num + sex_bin + aspiration + region_eur + year_cat,
           data = md, family = binomial)
cat(sprintf("[cat-year] n=%d events=%d\n", length(m4b$y), sum(m4b$y)))
cf <- coef(m4b); ci <- confint.default(m4b); pv <- coef(summary(m4b))[, "Pr(>|z|)"]
for (v in names(cf)) cat(sprintf("  %-20s OR %.3f (%.3f-%.3f) p=%.3f\n",
                                 v, exp(cf[v]), exp(ci[v,1]), exp(ci[v,2]), pv[v]))

cat("\n=== S1c: yearly trend excluding imputed years ===\n")
yearly_r <- md[year_imputed == FALSE, .(Total=.N, Fatal=sum(is_fatal)),
               by=event_year][order(event_year)]
print(yearly_r)
if (nrow(yearly_r) >= 2) {
  yr_g <- glm(cbind(Fatal, Total-Fatal) ~ event_year, data=yearly_r, family=binomial)
  cat(sprintf("Trend excl imputed: OR %.3f (%.3f-%.3f) p=%.3f\n",
              exp(coef(yr_g)["event_year"]),
              exp(confint.default(yr_g)["event_year",1]),
              exp(confint.default(yr_g)["event_year",2]),
              coef(summary(yr_g))["event_year","Pr(>|z|)"]))
}

cat("\n=== S1d: m4 with median-year imputation ===\n")
med_year <- median(md[year_imputed == FALSE, event_year], na.rm = TRUE)
cat(sprintf("Median observed year = %d\n", med_year))
md[, event_year_med := ifelse(year_imputed, med_year, event_year)]
m4d <- glm(is_fatal ~ age_num + sex_bin + aspiration + region_eur + event_year_med,
           data = md, family = binomial)
or_line(m4d, "median-year imputation")
c4d <- coef(m4d); ci4d <- confint.default(m4d)
cat(sprintf("  median-year OR per year %.3f (%.3f-%.3f) p=%.3f\n",
            exp(c4d["event_year_med"]), exp(ci4d["event_year_med",1]),
            exp(ci4d["event_year_med",2]),
            coef(summary(m4d))["event_year_med","Pr(>|z|)"]))

cat("\n=== S1e: yearly trend restricted to observed 2022-2025 ===\n")
yr_obs <- md[year_imputed == FALSE & event_year >= 2022,
             .(Total=.N, Fatal=sum(is_fatal)), by=event_year][order(event_year)]
print(yr_obs)
if (nrow(yr_obs) >= 2) {
  yrg <- glm(cbind(Fatal, Total-Fatal) ~ event_year, data=yr_obs, family=binomial)
  cat(sprintf("Trend 2022-2025 observed only: OR %.3f (%.3f-%.3f) p=%.3f\n",
              exp(coef(yrg)["event_year"]),
              exp(confint.default(yrg)["event_year",1]),
              exp(confint.default(yrg)["event_year",2]),
              coef(summary(yrg))["event_year","Pr(>|z|)"]))
}

cat("\n=== S2a: m4 without age ===\n")
m4_noage <- glm(is_fatal ~ sex_bin + aspiration + region_eur + event_year,
                data = md, family = binomial)
or_line(m4_noage, "no age")

cat("\n=== S2b: m4 with age missing-indicator ===\n")
md[, age_miss := as.integer(is.na(age_num))]
md[, age_fill := ifelse(is.na(age_num), median(age_num, na.rm = TRUE), age_num)]
m4_mi <- glm(is_fatal ~ age_fill + age_miss + sex_bin + aspiration + region_eur + event_year,
             data = md, family = binomial)
or_line(m4_mi, "age missing-indicator")
cmi <- coef(m4_mi); cimi <- confint.default(m4_mi)
cat(sprintf("  age_fill OR per 10yr %.3f (%.3f-%.3f) p=%.3f\n",
            exp(10*cmi["age_fill"]), exp(10*cimi["age_fill",1]), exp(10*cimi["age_fill",2]),
            coef(summary(m4_mi))["age_fill","Pr(>|z|)"]))
cat(sprintf("  age_miss OR %.3f (%.3f-%.3f) p=%.3f\n",
            exp(cmi["age_miss"]), exp(cimi["age_miss",1]), exp(cimi["age_miss",2]),
            coef(summary(m4_mi))["age_miss","Pr(>|z|)"]))

cat("\n=== S3: MICE (m=10) for age and sex, pooled m4 ===\n")
imp_md <- md[, .(is_fatal, age_num, sex_bin, aspiration, region_eur, event_year)]
set.seed(20260814)
imp <- mice(imp_md, m = 10, maxit = 10, method = c("", "pmm", "logreg", "", "", ""),
            printFlag = FALSE)
m4_mice <- with(imp, glm(is_fatal ~ age_num + sex_bin + aspiration + region_eur + event_year,
                         family = binomial))
pooled <- summary(pool(m4_mice), conf.int = TRUE)
cat(sprintf("Pooled n per imputation = %d\n", length(m4_mice$analyses[[1]]$y)))
for (v in c("age_num","sex_bin","aspiration","region_eur","event_year")) {
  row <- pooled[pooled$term == v, ]
  cat(sprintf("  %-12s OR %.3f (%.3f-%.3f) p=%.3f\n", v,
              exp(row$estimate), exp(row$conf.low), exp(row$conf.high), row$p.value))
}

cat("\n=== S4: caseid-level (latest version per case) ===\n")
cloz_by_case <- cloz_lung[order(caseid, -is_fatal, -caseversion_num)][
  , .SD[1], by = caseid]
cat(sprintf("Unique caseids: %d (vs %d unique primaryid)\n",
            nrow(cloz_by_case), nrow(cloz_by_id)))
cat(sprintf("Caseid-level fatality: %d/%d (%.1f%%)\n",
            sum(cloz_by_case$is_fatal), nrow(cloz_by_case),
            100*sum(cloz_by_case$is_fatal)/nrow(cloz_by_case)))
cloz_by_case[, region := "Other"]
cloz_by_case[reporter_country %in% c("GB","FR","DE","IT","ES","NL","BE","CH","AT","SE",
                                     "NO","DK","FI","IE","PT","GR","PL","CZ","RO","HU",
                                     "BG","HR","SK","LT","LV","EE","LU","MT","CY"), region := "Europe"]
cloz_by_case[reporter_country %in% c("US","CA"), region := "North America"]
cloz_by_case[reporter_country %in% c("AU","NZ"), region := "Oceania"]
cloz_by_case[reporter_country %in% c("JP","CN","KR","TW","HK","IN","TH","MY","SG",
                                     "ID","PH","VN","PK","BD"), region := "Asia"]
cloz_by_case[reporter_country %in% c("BR","AR","CL","CO","MX","PE","VE","UY","EC",
                                     "BO","PY"), region := "Latin America"]
mc <- cloz_by_case[, .(is_fatal, age_num, sex, is_aspiration, region, event_year)]
mc[, sex_bin := ifelse(sex == "M", 1, ifelse(sex == "F", 0, NA))]
mc[, aspiration := as.integer(is_aspiration)]
mc[, region_eur := ifelse(region == "Europe", 1, 0)]
m4_case <- glm(is_fatal ~ age_num + sex_bin + aspiration + region_eur + event_year,
               data = mc, family = binomial)
or_line(m4_case, "caseid-level m4")

cat("\n=== S5: age structure by region (for crude-rate confounding) ===\n")
for (r in c("Europe", "North America", "Oceania", "Asia", "Latin America", "Other")) {
  a <- cloz_by_id[region == r & !is.na(age_num), age_num]
  cat(sprintf("  %-14s n_age=%4d median=%.0f (IQR %.0f-%.0f)\n",
              r, length(a), median(a), quantile(a, 0.25), quantile(a, 0.75)))
}

cat("\n=== S6: crude vs adjusted Europe contrast ===\n")
eur_tab <- table(cloz_by_id$region == "Europe", cloz_by_id$is_fatal)
print(eur_tab)
cat(sprintf("Crude Europe OR = %.3f; chisq p = %.4f\n",
            (eur_tab[2,2]*eur_tab[1,1])/(eur_tab[1,2]*eur_tab[2,1]),
            chisq.test(eur_tab)$p.value))

# ---- Calibration plot (Supplementary Figure S2) ----
m4_main <- glm(is_fatal ~ age_num + sex_bin + aspiration + region_eur + event_year,
               data = md, family = binomial)
pred <- predict(m4_main, type = "response")
obs <- m4_main$y
png("F:/clazpin/manuscript/analysis/figures/figS2_calibration.png",
    width = 6.5, height = 6, units = "in", res = 300)
plot(pred, jitter(as.numeric(obs), amount = 0.02), pch = 19, cex = 0.35, col = rgb(0,0,0,0.35),
     xlab = "Predicted probability of fatal outcome",
     ylab = "Observed fatal outcome",
     main = sprintf("Calibration of the multivariable model (complete cases, n = %d)",
                    length(obs)))
abline(0, 1, lty = 2, col = "grey50")
smo <- loess(obs ~ pred, span = 0.8)
px <- seq(min(pred), max(pred), length.out = 100)
lines(px, predict(smo, newdata = data.frame(pred = px)), lwd = 2, col = "steelblue")
legend("topleft", legend = c("Ideal (45° line)", "Loess-smoothed observed"),
       col = c("grey50", "steelblue"), lty = c(2, 1), lwd = 2, bty = "n")
dev.off()
cat("\nCalibration plot written: figS2_calibration.png\n")

cat(sprintf("\nMICE version: %s\n", as.character(packageVersion("mice"))))
