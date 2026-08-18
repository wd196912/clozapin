# ============================================================
# 复合终点主分析 + 主要诊断分配（呼吸科专家第二轮意见实施方案）
# 数据源: FAERS ASCII 2022Q1-2025Q4 预提取文件 (F:/faersdata)
# 输出: composite_endpoint_signal_results.csv
#       principal_diagnosis_table1.csv
#       composite_endpoint_log.txt
# ============================================================

suppressMessages({library(data.table); library(dplyr)})

data_dir <- "F:/faersdata"
out_dir  <- "F:/clazpin"
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

logf <- file(file.path(out_dir, "composite_endpoint_log.txt"), "w")
sink(logf, append = TRUE); sink(logf, append = TRUE, type = "message")
cat("=== 复合终点主分析 (2026-08-18) ===\n\n")

# ---------- 药物分类 ----------
drug_raw[, drug_type := ifelse(
  grepl("CLOZAPINE|CLOZARIL|FAZACLO|VERSACLOZ|LEPONEX|ZAPONEX|CLOPINE|DENZAPINE",
        drugname, ignore.case = TRUE) &
    !grepl("RISPERIDONE|RISPERDAL|PERSERIS", drugname, ignore.case = TRUE),
  "Clozapine", "Other")]
drug_raw[grepl("RISPERIDONE|RISPERDAL|PERSERIS|RISPERDAL CONSTA",
               drugname, ignore.case = TRUE) & drug_type != "Clozapine",
         drug_type := "Risperidone"]
drug_ps <- drug_raw[role_cod == "PS"]
cloz_ids <- unique(drug_ps[drug_type == "Clozapine", primaryid])
risp_ids <- unique(drug_ps[drug_type == "Risperidone", primaryid])

cat("GATE cloz_ids =", length(cloz_ids), "(expected 44055)\n")
cat("GATE risp_ids =", length(risp_ids), "(expected 15130)\n")
stopifnot(length(cloz_ids) == 44055, length(risp_ids) == 15130)

# ---------- 死亡定义 (OUTC DE) ----------
outc_files <- list.files(qdir, pattern = "^OUTC.*\\.txt$", full.names = TRUE)
outc <- rbindlist(lapply(outc_files, function(f) {
  fread(f, sep = "$", header = FALSE, select = 1:3,
        col.names = c("primaryid", "caseid", "outc_cod"),
        na.strings = "", blank.lines.skip = TRUE, strip.white = TRUE)
}))
death_ids <- unique(outc[outc_cod == "DE", primaryid])

# ---------- 宽泛筛查 (复现 1305) ----------
lung_terms_broad <- c(
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
  "Pneumonia haemophilus$", "Pneumonia legionella$",
  "Pneumonia chlamydial$", "Pneumonia moraxella$",
  "Eosinophilic pneumonia$", "Pneumonia anthrax$",
  "Radiation pneumonitis$", "Pulmonary mycosis$",
  "Pneumonia cytomegaloviral$", "Pneumonia respiratory syncytial viral$",
  "Pneumonia adenoviral$", "Pneumonia herpes viral$",
  "Hypersensitivity pneumonitis$", "Pulmonary nocardiosis$",
  "Idiopathic interstitial pneumonia$", "Acute interstitial pneumonitis$",
  "Interstitial pneumonitis$", "Pneumonia chemical$",
  "Pneumonia lipid$", "Metapneumovirus pneumonia$",
  "Pneumonia parainfluenzae viral$", "Upper respiratory tract infection$",
  "Lower respiratory tract infection viral$",
  "Respiratory tract infection viral$",
  "Respiratory tract infection bacterial$",
  "Viral upper respiratory tract infection$",
  "Respiratory moniliasis$", "Infective exacerbation of bronchiectasis$",
  "Pulmonary actinomycosis$", "Pneumonia bacterial NOS$",
  "Gallbladder empyema$")
lung_pattern <- paste(lung_terms_broad, collapse = "|")
reac_raw[, is_lung_broad := grepl(lung_pattern, pt, ignore.case = TRUE) &
              !grepl("Gallbladder empyema", pt, ignore.case = TRUE)]
lung_ids_broad <- unique(reac_raw[is_lung_broad == TRUE, primaryid])
cohort_broad <- length(intersect(cloz_ids, lung_ids_broad))
cat("GATE broad cohort =", cohort_broad, "(expected 1305)\n")
stopifnot(cohort_broad == 1305)

# ---------- 26 PT 筛查清单 (稿件如实清单) ----------
pt17 <- c("pneumonia","pneumonia aspiration","lower respiratory tract infection",
          "upper respiratory tract infection","covid-19 pneumonia",
          "respiratory tract infection","pneumonia bacterial","pneumonia viral",
          "pneumonitis","empyema","idiopathic interstitial pneumonia",
          "lung abscess","pulmonary tuberculosis","pneumonia klebsiella",
          "respiratory tract infection viral","pneumonia influenzal",
          "pneumonia staphylococcal")
pt_extra <- c("pneumonia necrotising","pneumonia anthrax","pulmonary sepsis",
              "acute interstitial pneumonitis","atypical pneumonia",
              "eosinophilic pneumonia","hypersensitivity pneumonitis",
              "organising pneumonia","pleural infection")
pt26 <- c(pt17, pt_extra)

reac_raw[, pt_l := tolower(trimws(pt))]
reac_raw[, is_26 := pt_l %in% pt26]
ids26 <- unique(reac_raw[is_26 == TRUE, primaryid])
cohort26 <- length(intersect(cloz_ids, ids26))
cat("GATE 26-PT cohort =", cohort26, "(expected 1305)\n")
stopifnot(cohort26 == 1305)

# ---------- 复合终点与亚组定义 ----------
pt_composite <- c("pneumonia","pneumonia aspiration","lower respiratory tract infection",
                  "pneumonia bacterial","pneumonia viral","empyema","lung abscess",
                  "pneumonia klebsiella","pneumonia influenzal","pneumonia staphylococcal")
pt_pneumonia_spectrum <- setdiff(pt_composite, "lower respiratory tract infection")
pt_lrti <- "lower respiratory tract infection"
pt_special_other <- setdiff(pt26, pt_composite)

# ---------- 信号计算 (与旧脚本同公式, 保证可比) ----------
calc_signals <- function(a, b, c, d, label = "") {
  ROR <- (a / b) / (c / d)
  se_ln_ROR <- sqrt(1/a + 1/b + 1/c + 1/d)
  ln_ROR <- log(ROR)
  ROR_ci_low <- exp(ln_ROR - 1.96 * se_ln_ROR)
  ROR_ci_high <- exp(ln_ROR + 1.96 * se_ln_ROR)
  ROR_p <- 2 * pnorm(-abs(ln_ROR / se_ln_ROR))
  PRR <- (a / (a + b)) / (c / (c + d))
  se_ln_PRR <- sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d))
  ln_PRR <- log(PRR)
  PRR_ci_low <- exp(ln_PRR - 1.96 * se_ln_PRR)
  PRR_ci_high <- exp(ln_PRR + 1.96 * se_ln_PRR)
  N <- a + b + c + d
  E_ic <- (a + b) * (a + c) / N
  IC <- log2((a + 0.5) / (E_ic + 0.5))
  IC_sd <- sqrt(1 / (log(2)^2) * (1/(a+0.5) + 1/(E_ic+0.5)))
  IC_025 <- IC - 1.96 * IC_sd
  IC_975 <- IC + 1.96 * IC_sd
  # 新规则: 对照组 n<5 不报 ROR/PRR; 阳性 = IC025>0 且对照组 n>=5
  report_ror <- c >= 5 && a >= 1
  if (is.infinite(ROR) || !report_ror) {
    ROR_s <- "NR"; ROR_low_s <- "NR"; ROR_high_s <- "NR"
    PRR_s <- "NR"; PRR_low_s <- "NR"; PRR_high_s <- "NR"
  } else {
    ROR_s <- round(ROR, 2); ROR_low_s <- round(ROR_ci_low, 2); ROR_high_s <- round(ROR_ci_high, 2)
    PRR_s <- round(PRR, 2); PRR_low_s <- round(PRR_ci_low, 2); PRR_high_s <- round(PRR_ci_high, 2)
  }
  positive <- (c >= 5) && (IC_025 > 0) && is.finite(IC_025)
  data.frame(
    Label = label, a = a, b = b, c = c, d = d, N_total = N,
    ROR = ROR_s, ROR_CI95_low = ROR_low_s, ROR_CI95_high = ROR_high_s,
    PRR = PRR_s, PRR_CI95_low = PRR_low_s, PRR_CI95_high = PRR_high_s,
    IC = round(IC, 3), IC_025 = round(IC_025, 3), IC_975 = round(IC_975, 3),
    P_value = if (report_ror) format.pval(ROR_p, digits = 3, eps = 1e-10) else "NR",
    Signal = ifelse(positive, "Positive", "Not positive"),
    ROR_reported = report_ror,
    ROR_raw = ifelse(is.finite(ROR), ROR, NA_real_),
    stringsAsFactors = FALSE)
}

two_by_two <- function(pt_set, cloz_ids, risp_ids, reac_raw) {
  ids <- unique(reac_raw[pt_l %in% pt_set, primaryid])
  a <- length(intersect(cloz_ids, ids))
  b <- length(setdiff(cloz_ids, ids))
  c <- length(intersect(risp_ids, ids))
  d <- length(setdiff(risp_ids, ids))
  list(a = a, b = b, c = c, d = d)
}

# ---------- 验证闸门: 复现旧逐 PT 数字 ----------
cat("\n=== 验证闸门: 旧 17 PT 数字复现 (vs clozapine_lung_signal_detection.csv) ===\n")
gate_file <- file.path(data_dir, "clozapine_lung_signal_detection.csv")
if (!file.exists(gate_file)) {
  cat("  闸门文件未找到; 先运行 01_signal_detection.R 生成该文件, 本次跳过复现验证\n")
} else {
  old_csv <- fread(gate_file)
  old_csv[, PT_l := tolower(trimws(PT_Term))]
  gate_ok <- TRUE
  for (i in seq_len(nrow(old_csv))) {
    p <- old_csv$PT_l[i]
    tt <- two_by_two(p, cloz_ids, risp_ids, reac_raw)
    exp_a <- old_csv$Clozapine_Event[i]; exp_c <- old_csv$Risperidone_Event[i]
    ok <- (tt$a == exp_a) && (tt$c == exp_c)
    if (!ok) gate_ok <- FALSE
    cat(sprintf("  %-40s a=%d(%d) c=%d(%d) %s\n", p, tt$a, exp_a, tt$c, exp_c,
                ifelse(ok, "OK", "MISMATCH")))
  }
  stopifnot(gate_ok)
  cat("  17 PT 复现: ALL OK\n")
}

# ---------- 新规则信号分析 ----------
cat("\n=== 新规则信号分析 (clozapine vs risperidone) ===\n")

rows <- list()
add_row <- function(pt_set, label) {
  tt <- two_by_two(pt_set, cloz_ids, risp_ids, reac_raw)
  rows[[length(rows) + 1]] <<- calc_signals(tt$a, tt$b, tt$c, tt$d, label)
}
add_row(pt_composite, "Composite pulmonary infection endpoint (primary)")
add_row(pt_pneumonia_spectrum, "Pneumonia spectrum (subgroup)")
add_row(pt_lrti, "Non-pneumonia lower respiratory tract infection (subgroup)")
add_row(pt_special_other, "Special/other pulmonary infection terms (subgroup)")
for (p in pt17) add_row(p, p)
add_row(pt_extra, "Other pulmonary infection terms (9 additional PTs)")

sig_df <- rbindlist(rows)
sig_df[, ROR_backsolve := ifelse(is.finite(ROR_raw),
                                 round((a * d) / (b * c), 2), NA_real_)]
sig_df[, backsolve_ok := is.na(ROR_raw) | (ROR == "NR") | (abs(ROR_backsolve - ROR_raw) < 0.011)]
print(sig_df[, .(Label, a, c, ROR, ROR_CI95_low, ROR_CI95_high, IC, IC_025, Signal, ROR_backsolve, backsolve_ok)])
stopifnot(all(sig_df$backsolve_ok))
fwrite(sig_df, file.path(out_dir, "composite_endpoint_signal_results.csv"))
cat("\nROR 回算: ALL OK; CSV written\n")

# ---------- 主要诊断分配 (26 PT 全链条) ----------
cat("\n=== 主要诊断分配 (n=1305) ===\n")

priority_chain <- list(
  "Pneumonia aspiration" = "pneumonia aspiration",
  "Pneumonia" = c("pneumonia","pneumonia bacterial","pneumonia viral",
                  "pneumonia klebsiella","pneumonia influenzal",
                  "pneumonia staphylococcal","pneumonia necrotising",
                  "pneumonia anthrax","atypical pneumonia",
                  "eosinophilic pneumonia","organising pneumonia"),
  "Lower respiratory tract infection" = "lower respiratory tract infection",
  "COVID-19 pneumonia" = "covid-19 pneumonia",
  "Upper respiratory tract infection" = "upper respiratory tract infection",
  "Respiratory tract infection" = c("respiratory tract infection",
                                    "respiratory tract infection viral"),
  "Pneumonitis" = c("pneumonitis","acute interstitial pneumonitis",
                    "hypersensitivity pneumonitis"),
  "Empyema" = c("empyema","pleural infection"),
  "Lung abscess" = "lung abscess",
  "Pulmonary tuberculosis" = "pulmonary tuberculosis",
  "Pulmonary sepsis" = "pulmonary sepsis",
  "Idiopathic interstitial pneumonia" = "idiopathic interstitial pneumonia")

cohort_ids <- intersect(cloz_ids, ids26)
pt_by_id <- reac_raw[primaryid %in% cohort_ids & is_26 == TRUE,
                     .(pt_l = unique(pt_l)), by = primaryid]

assign_principal <- function(pts, chain) {
  for (nm in names(chain)) {
    if (any(chain[[nm]] %in% pts)) return(nm)
  }
  return(NA_character_)
}
pt_by_id[, principal := assign_principal(pt_l, priority_chain), by = primaryid]
stopifnot(sum(is.na(pt_by_id$principal)) == 0)

pt_by_id[, is_fatal := primaryid %in% death_ids]
pt_uniq <- unique(pt_by_id, by = "primaryid")
cat("GATE cohort rows =", nrow(pt_uniq), "(expected 1305)\n")
cat("GATE fatal total =", sum(pt_uniq$is_fatal), "(expected 309)\n")
stopifnot(nrow(pt_uniq) == 1305, sum(pt_uniq$is_fatal) == 309)

tab1 <- pt_uniq[, .(n = .N, fatal = sum(is_fatal)), by = principal]
tab1[, nonfatal := n - fatal]
tab1[, pct_fatal := round(100 * fatal / n, 1)]
tab1 <- tab1[order(-fatal)]
cat("\n主要诊断分类 (Table 1 感染类型行):\n")
print(tab1)
cat(sprintf("\n分区合计 = %d (expected 1305)\n", sum(tab1$n)))
cat(sprintf("死亡合计 = %d (expected 309)\n", sum(tab1$fatal)))
stopifnot(sum(tab1$n) == 1305, sum(tab1$fatal) == 309)
fwrite(tab1, file.path(out_dir, "principal_diagnosis_table1.csv"))

# ---------- 复合终点队列关键数字 ----------
cat("\n=== 复合终点关键数字 ===\n")
comp_ids <- unique(reac_raw[pt_l %in% pt_composite, primaryid])
comp_cloz <- length(intersect(cloz_ids, comp_ids))
comp_risp <- length(intersect(risp_ids, comp_ids))
comp_in_cohort <- length(intersect(cohort_ids, comp_ids))
cat(sprintf("复合终点: cloz %d / risp %d; 1305 队列内 %d (%.1f%%)\n",
            comp_cloz, comp_risp, comp_in_cohort, 100 * comp_in_cohort / 1305))
comp_fatal <- sum(intersect(cohort_ids, comp_ids) %in% death_ids)
cat(sprintf("复合终点队列内 fatal = %d\n", comp_fatal))

# ---------- 补充数字 (稿件引用) ----------
cat("\n=== 补充数字 ===\n")
ids17 <- unique(reac_raw[pt_l %in% pt17, primaryid])
risp_cohort17 <- length(intersect(risp_ids, ids17))
risp_cohort26 <- length(intersect(risp_ids, ids26))
cat(sprintf("risp 17-PT cohort = %d; risp 26-PT cohort = %d\n",
            risp_cohort17, risp_cohort26))
cloz17 <- length(intersect(cloz_ids, ids17))
cloz26 <- length(intersect(cloz_ids, ids26))
cat(sprintf("cloz 17-PT cohort = %d; cloz 26-PT cohort = %d\n", cloz17, cloz26))
spec_ids <- unique(reac_raw[pt_l %in% pt_pneumonia_spectrum, primaryid])
cat(sprintf("pneumonia spectrum: cloz %d / risp %d\n",
            length(intersect(cloz_ids, spec_ids)), length(intersect(risp_ids, spec_ids))))
lrti_ids <- unique(reac_raw[pt_l %in% pt_lrti, primaryid])
cat(sprintf("LRTI subgroup: cloz %d / risp %d\n",
            length(intersect(cloz_ids, lrti_ids)), length(intersect(risp_ids, lrti_ids))))
so_ids <- unique(reac_raw[pt_l %in% pt_special_other, primaryid])
cat(sprintf("special/other: cloz %d / risp %d\n",
            length(intersect(cloz_ids, so_ids)), length(intersect(risp_ids, so_ids))))
# 复合终点各 PT 的 risp 计数 (供 S4 方案A)
for (p in sort(unique(c(pt17, pt_extra)))) {
  c_p <- length(intersect(risp_ids, unique(reac_raw[pt_l == p, primaryid])))
  cat(sprintf("risp %-40s = %d\n", p, c_p))
}

sink(); sink(type = "message")
cat("DONE\n")
