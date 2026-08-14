# ============================================================
# 氯氮平所致肺部感染 — 深度分析
# 包含: 信号检测(ROR/PRR/IC)、亚组分析、时间趋势、利培酮对比、可视化
# ============================================================

library(data.table)
library(dplyr)
library(tidyr)

data_dir <- "F:/faersdata"
setwd(data_dir)

# ============================================================
# Part 0: 数据导入 + 去重
# ============================================================
message("=== 导入数据 ===")

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

# ============================================================
# Part 1: 信号检测 — 氯氮平 vs 利培酮 (肺部感染)
# ============================================================
message("\n=========================================")
message("Part 1: 药物警戒信号检测 (Disproportionality Analysis)")
message("=========================================")

# 标记药物类型
drug_raw[, drug_type := ifelse(
  grepl("CLOZAPINE|CLOZARIL|FAZACLO|VERSACLOZ|LEPONEX|ZAPONEX|CLOPINE|DENZAPINE",
        drugname, ignore.case = TRUE) &
    !grepl("RISPERIDONE|RISPERDAL|PERSERIS", drugname, ignore.case = TRUE),
  "Clozapine", "Other")]

drug_raw[grepl("RISPERIDONE|RISPERDAL|PERSERIS|RISPERDAL CONSTA",
               drugname, ignore.case = TRUE) &
           drug_type != "Clozapine",
         drug_type := "Risperidone"]

message(sprintf("药物分布: Clozapine=%d, Risperidone=%d, Other=%d",
                sum(drug_raw$drug_type == "Clozapine"),
                sum(drug_raw$drug_type == "Risperidone"),
                sum(drug_raw$drug_type == "Other")))

# 肺部感染 PT 术语列表
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
  "Pneumonia haemophilus$", "Pneumonia legionella$",
  "Pneumonia chlamydial$", "Pneumonia moraxella$",
  "Eosinophilic pneumonia$", "Pneumonia anthrax$",
  "Radiation pneumonitis$", "Pulmonary mycosis$",
  "Pneumonia cytomegaloviral$",
  "Pneumonia respiratory syncytial viral$",
  "Pneumonia adenoviral$", "Pneumonia herpes viral$",
  "Hypersensitivity pneumonitis$", "Pulmonary nocardiosis$",
  "Idiopathic interstitial pneumonia$",
  "Acute interstitial pneumonitis$",
  "Interstitial pneumonitis$", "Pneumonia chemical$",
  "Pneumonia lipid$", "Metapneumovirus pneumonia$",
  "Pneumonia parainfluenzae viral$",
  "Upper respiratory tract infection$",  # 上呼吸道感染
  "Lower respiratory tract infection viral$",
  "Respiratory tract infection viral$",
  "Respiratory tract infection bacterial$",
  "Viral upper respiratory tract infection$",
  "Respiratory moniliasis$",
  "Infective exacerbation of bronchiectasis$",
  "Pulmonary actinomycosis$",
  "Pneumonia bacterial NOS$",
  "Gallbladder empyema$"  # 排除
)

lung_pattern <- paste(lung_terms, collapse = "|")

# 标记肺部感染
reac_raw[, is_lung := grepl(lung_pattern, pt, ignore.case = TRUE) &
              !grepl("Gallbladder empyema", pt, ignore.case = TRUE)]

message(sprintf("肺部感染记录: %d / %d", sum(reac_raw$is_lung), nrow(reac_raw)))

# --- 信号检测: 2x2 列联表 ---

# 策略: 在药物-反应水平上计算(每个primaryid可有多条drug和reac)
# 使用唯一 primaryid 去重

# 获取所有 PS 药物记录
drug_ps <- drug_raw[role_cod == "PS"]

# 每个 primaryid 是否涉及氯氮平/利培酮
cloz_ids <- unique(drug_ps[drug_type == "Clozapine", primaryid])
risp_ids <- unique(drug_ps[drug_type == "Risperidone", primaryid])

# 每个 primaryid 是否报告了肺部感染
lung_ids <- unique(reac_raw[is_lung == TRUE, primaryid])

# 所有出现过的 primaryid
all_ids <- unique(c(drug_ps$primaryid))

message(sprintf("\n整体 primaryid 数: %d", length(all_ids)))
message(sprintf("氯氮平相关 primaryid: %d", length(cloz_ids)))
message(sprintf("利培酮相关 primaryid: %d", length(risp_ids)))
message(sprintf("肺部感染相关 primaryid: %d", length(lung_ids)))

# --- ROR & PRR 计算函数 ---
calc_signals <- function(a, b, c, d, label = "") {
  # a: drug+event, b: drug+no_event, c: no_drug+event, d: no_drug+no_event

  # ROR
  ROR <- (a / b) / (c / d)
  se_ln_ROR <- sqrt(1/a + 1/b + 1/c + 1/d)
  ln_ROR <- log(ROR)
  ROR_ci_low <- exp(ln_ROR - 1.96 * se_ln_ROR)
  ROR_ci_high <- exp(ln_ROR + 1.96 * se_ln_ROR)
  ROR_p <- 2 * pnorm(-abs(ln_ROR / se_ln_ROR))

  # PRR
  PRR <- (a / (a + b)) / (c / (c + d))
  se_ln_PRR <- sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d))
  ln_PRR <- log(PRR)
  PRR_ci_low <- exp(ln_PRR - 1.96 * se_ln_PRR)
  PRR_ci_high <- exp(ln_PRR + 1.96 * se_ln_PRR)

  # Chi-square
  N <- a + b + c + d
  expected <- (a + b) * (a + c) / N
  chi_sq <- (a - expected)^2 / expected
  if (a < expected) chi_sq <- -chi_sq  # 负值表示低于预期

  # Information Component (IC) - Bayesian
  # IC = log2((a + 0.5) / E), E = (a+b)*(a+c)/(N)
  E_ic <- (a + b) * (a + c) / N
  IC <- log2((a + 0.5) / (E_ic + 0.5))
  # IC 标准差简化
  IC_sd <- sqrt(1 / (log(2)^2) * (1/(a+0.5) + 1/(E_ic+0.5)))
  IC_025 <- IC - 1.96 * IC_sd
  IC_975 <- IC + 1.96 * IC_sd

  data.frame(
    Label = label,
    a = a, b = b, c = c, d = d,
    N_total = N,
    ROR = round(ROR, 2),
    ROR_CI95_low = round(ROR_ci_low, 2),
    ROR_CI95_high = round(ROR_ci_high, 2),
    PRR = round(PRR, 2),
    PRR_CI95_low = round(PRR_ci_low, 2),
    PRR_CI95_high = round(PRR_ci_high, 2),
    Chi_sq = round(chi_sq, 2),
    IC = round(IC, 3),
    IC_025 = round(IC_025, 3),
    IC_975 = round(IC_975, 3),
    P_value = format.pval(ROR_p, digits = 3, eps = 1e-10),
    Signal = ifelse(ROR > 1 & ROR_ci_low > 1 & a >= 3 & IC_025 > 0,
                    "Positive Signal", "No Signal"),
    stringsAsFactors = FALSE
  )
}

# 氯氮平 vs 利培酮 肺部感染信号
cloz_lung_n  <- length(intersect(cloz_ids, lung_ids))
cloz_nolung  <- length(setdiff(cloz_ids, lung_ids))
risp_lung_n  <- length(intersect(risp_ids, lung_ids))
risp_nolung  <- length(setdiff(risp_ids, lung_ids))

message("\n--- 氯氮平 vs 利培酮 肺部感染 2x2表 ---")
message(sprintf("              肺部感染(+)  肺部感染(-)  合计"))
message(sprintf("氯氮平         %-10d  %-10d  %d", cloz_lung_n, cloz_nolung, cloz_lung_n+cloz_nolung))
message(sprintf("利培酮         %-10d  %-10d  %d", risp_lung_n, risp_nolung, risp_lung_n+risp_nolung))
message(sprintf("合计           %-10d  %-10d  %d",
                cloz_lung_n+risp_lung_n, cloz_nolung+risp_nolung,
                cloz_lung_n+cloz_nolung+risp_lung_n+risp_nolung))

signal_result <- calc_signals(cloz_lung_n, cloz_nolung, risp_lung_n, risp_nolung,
                               "Clozapine vs Risperidone — 肺部感染(总)")

# --- 按具体 PT 术语做信号检测 ---
message("\n--- 各肺部感染 PT 的 ROR 信号 ---")

top_lung_pts <- names(sort(table(reac_raw[is_lung == TRUE]$pt), decreasing = TRUE))
# 过滤掉样本太少的
top_lung_pts <- top_lung_pts[1:min(20, length(top_lung_pts))]

signal_by_pt <- list()
for (this_pt in top_lung_pts) {
  pt_lung_ids <- unique(reac_raw[pt == this_pt, primaryid])

  a <- length(intersect(cloz_ids, pt_lung_ids))
  b <- length(setdiff(cloz_ids, pt_lung_ids))
  c <- length(intersect(risp_ids, pt_lung_ids))
  d <- length(setdiff(risp_ids, pt_lung_ids))

  if (a >= 1) {
    signal_by_pt[[this_pt]] <- calc_signals(a, b, c, d, this_pt)
  }
}

signal_pt_df <- rbindlist(signal_by_pt)
signal_pt_df <- signal_pt_df[order(-ROR)]

cat("\n各 PT 的信号检测结果 (氯氮平 vs 利培酮):\n")
print(signal_pt_df[, .(Label, a, b, c, d, ROR, ROR_CI95_low, ROR_CI95_high,
                         IC, IC_025, IC_975, Signal)])

# ============================================================
# Part 2: 氯氮平肺部感染 — 人口学亚组分析
# ============================================================
message("\n=========================================")
message("Part 2: 氯氮平肺部感染人口学亚组分析")
message("=========================================")

# 构建氯氮平肺部感染完整数据
drug_cloz_ps <- drug_ps[drug_type == "Clozapine"]
cloz_lung_data <- reac_raw[is_lung == TRUE] %>%
  inner_join(drug_cloz_ps, by = c("primaryid", "caseid"),
             relationship = "many-to-many") %>%
  inner_join(demo_raw, by = c("primaryid", "caseid"),
             relationship = "many-to-many") %>%
  as.data.table()

message(sprintf("氯氮平肺部感染合并数据: %d 行, %d 独特 primaryid",
                nrow(cloz_lung_data), uniqueN(cloz_lung_data$primaryid)))

# 数值化年龄
cloz_lung_data[, age_num := as.numeric(age)]
cloz_lung_data[age_cod == "DEC", age_num := age_num / 12]  # 月→岁(极少)
cloz_lung_data[age_cod == "WK" | age_cod == "WEEK", age_num := age_num / 52]
cloz_lung_data[age_cod == "DY" | age_cod == "DAY", age_num := age_num / 365]

# 年龄分组 (更细)
cloz_lung_data[, age_group := cut(age_num,
  breaks = c(-Inf, 18, 30, 40, 50, 60, 70, 80, Inf),
  labels = c("<18", "18-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"),
  right = TRUE
)]

cat("\n--- 年龄分组 × 性别 交叉表 ---\n")
age_sex_tbl <- cloz_lung_data[, .N, by = .(age_group, sex)][order(age_group, sex)]
print(dcast(age_sex_tbl, age_group ~ sex, value.var = "N", fill = 0))

cat("\n--- 年龄描述统计 ---\n")
print(summary(cloz_lung_data$age_num))

cat("\n--- 体重描述统计 ---\n")
cloz_lung_data[, wt_num := as.numeric(wt)]
print(summary(cloz_lung_data$wt_num))

# --- 严重结局分析 ---
cat("\n--- 严重结局 (I/F code) ---\n")
# F = Fatal, I = 非严重(但可能严重未报告), 实际上i_f_code中F表示严重
severity_tbl <- cloz_lung_data[, .N, by = .(i_f_code, pt)][order(i_f_code, -N)]
print(severity_tbl[1:30])

cat("\n严重结局 × PT交叉表:\n")
sev_pt <- dcast(cloz_lung_data[, .N, by = .(i_f_code, pt)],
                pt ~ i_f_code, value.var = "N", fill = 0)
sev_pt[, total := rowSums(.SD), .SDcols = patterns("^[FI]$|^[A-Z]$")]
print(sev_pt[order(-total)][1:20])

# --- 年龄与感染类型关联 ---
cat("\n--- 各年龄组 Top 3 PT ---\n")
age_pt_top3 <- cloz_lung_data[!is.na(age_group), .N, by = .(age_group, pt)][
  order(age_group, -N)
]
age_pt_top3 <- age_pt_top3[, .SD[1:min(3, .N)], by = age_group]
print(age_pt_top3)

# ============================================================
# Part 3: 时间趋势分析
# ============================================================
message("\n=========================================")
message("Part 3: 时间趋势分析")
message("=========================================")

# 提取年份和季度
cloz_lung_data[, event_year := as.integer(substr(event_dt, 1, 4))]
cloz_lung_data[, event_month := as.integer(substr(event_dt, 5, 6))]
cloz_lung_data[is.na(event_month) | event_month == 0 | event_month > 12, event_month := 6L]
cloz_lung_data[is.na(event_year) | event_year < 1900 | event_year > 2100, c("event_year", "event_month") := .(2022L, 6L)]
cloz_lung_data[, event_ym := as.Date(sprintf("%04d-%02d-15", event_year, event_month))]
cloz_lung_data[, event_qtr := paste0(event_year, "Q", ceiling(event_month / 3))]

# 只保留 2015-2025 的数据(有效时间窗口)
trend_data <- cloz_lung_data[event_year >= 2015 & event_year <= 2025]

cat("\n--- 年度报告趋势 ---\n")
yearly <- trend_data[, .(
  Reports = uniqueN(primaryid),
  Records = .N
), by = event_year][order(event_year)]
print(yearly)

cat("\n--- 季度报告趋势 ---\n")
qtrly <- trend_data[, .(
  Reports = uniqueN(primaryid),
  Records = .N
), by = event_qtr][order(event_qtr)]
print(qtrly)

# --- 2020年前后对比 (COVID前/后) ---
pre_covid  <- trend_data[event_year %in% 2015:2019]
post_covid <- trend_data[event_year %in% 2020:2025]

cat("\n--- COVID前后对比 ---\n")
cat(sprintf("COVID前 (2015-2019): %d 报告, %d 记录\n",
            uniqueN(pre_covid$primaryid), nrow(pre_covid)))
cat(sprintf("COVID后 (2020-2025): %d 报告, %d 记录\n",
            uniqueN(post_covid$primaryid), nrow(post_covid)))

# 各PT年度趋势
cat("\n--- 主要 PT 年度趋势 ---\n")
top_pts <- names(sort(table(cloz_lung_data$pt), decreasing = TRUE))[1:6]
for (tp in top_pts) {
  cat(sprintf("\n  %s:\n", tp))
  pt_yearly <- cloz_lung_data[pt == tp & event_year >= 2015 & event_year <= 2025,
                               .(Reports = uniqueN(primaryid)),
                               by = event_year][order(event_year)]
  print(pt_yearly)
}

# ============================================================
# Part 4: 吸入性肺炎专项分析
# ============================================================
message("\n=========================================")
message("Part 4: 吸入性肺炎 (Aspiration Pneumonia) 专项")
message("=========================================")

# 吸入性肺炎是氯氮平唾液分泌过多的已知严重并发症
asp_data <- cloz_lung_data[grepl("aspiration|aspiratio", pt, ignore.case = TRUE) |
                             pt %in% c("Pneumonia aspiration", "Aspiration pneumonia",
                                        "Pneumonitis aspiration")]

message(sprintf("吸入性肺炎记录: %d 行, %d 独特报告",
                nrow(asp_data), uniqueN(asp_data$primaryid)))

if (nrow(asp_data) > 0) {
  cat("\n吸入性肺炎 — 性别分布:\n")
  print(table(asp_data$sex, useNA = "ifany"))

  cat("\n吸入性肺炎 — 年龄组分布:\n")
  print(table(asp_data$age_group, useNA = "ifany"))

  cat("\n吸入性肺炎 — 年度趋势:\n")
  asp_yearly <- asp_data[event_year >= 2015, .(N = uniqueN(primaryid)), by = event_year]
  print(asp_yearly[order(event_year)])

  cat("\n吸入性肺炎 — 报告国家 Top 10:\n")
  asp_country <- sort(table(asp_data$reporter_country), decreasing = TRUE)
  print(asp_country[1:min(10, length(asp_country))])

  cat("\n吸入性肺炎 — 严重结局:\n")
  print(table(asp_data$i_f_code, useNA = "ifany"))

  cat("\n吸入性肺炎 — Dechallenge/Rechallenge:\n")
  cat("Dechallenge:\n")
  print(table(asp_data$dechal, useNA = "ifany"))
  cat("Rechallenge:\n")
  print(table(asp_data$rechal, useNA = "ifany"))

  # 吸入性肺炎信号 (vs 利培酮)
  asp_lung_ids <- unique(asp_data$primaryid)
  asp_all_pt_ids <- unique(reac_raw[grepl("aspiration|aspiratio", pt, ignore.case = TRUE) |
                                     pt %in% c("Pneumonia aspiration", "Aspiration pneumonia",
                                               "Pneumonitis aspiration"), primaryid])

  a_asp <- length(intersect(cloz_ids, asp_all_pt_ids))
  b_asp <- length(setdiff(cloz_ids, asp_all_pt_ids))
  c_asp <- length(intersect(risp_ids, asp_all_pt_ids))
  d_asp <- length(setdiff(risp_ids, asp_all_pt_ids))

  cat("\n吸入性肺炎信号检测 (氯氮平 vs 利培酮):\n")
  asp_signal <- calc_signals(a_asp, b_asp, c_asp, d_asp, "Aspiration Pneumonia")
  print(asp_signal[, c("Label", "a", "b", "c", "d", "ROR", "ROR_CI95_low",
                         "ROR_CI95_high", "PRR", "IC", "IC_025", "Signal")])
}

# ============================================================
# Part 5: 合并用药分析
# ============================================================
message("\n=========================================")
message("Part 5: 合并用药分析")
message("=========================================")

# 氯氮平肺部感染病例中,有哪些伴随药物?
cloz_lung_caseids <- unique(cloz_lung_data$caseid)

# 这些case中所有的药物记录(不仅仅是氯氮平)
conmed_data <- drug_raw[caseid %in% cloz_lung_caseids &
                          !grepl("CLOZAPINE|CLOZARIL|FAZACLO|VERSACLOZ|LEPONEX|ZAPONEX|CLOPINE",
                                 drugname, ignore.case = TRUE)]

message(sprintf("氯氮平肺部感染病例中的伴随药物记录: %d 行", nrow(conmed_data)))

cat("\n--- 最常见伴随药物 Top 30 ---\n")
conmed_top <- sort(table(conmed_data$drugname), decreasing = TRUE)
print(conmed_top[1:min(30, length(conmed_top))])

cat("\n--- 最常见伴随活性成分 Top 20 ---\n")
conmed_ai <- sort(table(conmed_data$prod_ai), decreasing = TRUE)
print(conmed_ai[1:min(20, length(conmed_ai))])

# 抗精神病药合并使用
antipsychotics <- c("OLANZAPINE", "QUETIAPINE", "ARIPIPRAZOLE", "HALOPERIDOL",
                     "PALIPERIDONE", "LURASIDONE", "BREXPIPRAZOLE", "ZIPRASIDONE",
                     "AMISULPRIDE", "CHLORPROMAZINE", "LEVOMEPROMAZINE",
                     "ASENAPINE", "CARIPRAZINE", "ILOPERIDONE")
ap_pattern <- paste(antipsychotics, collapse = "|")
conmed_ap <- conmed_data[grepl(ap_pattern, prod_ai, ignore.case = TRUE)]
ap_tbl <- sort(table(conmed_ap$prod_ai), decreasing = TRUE)
cat("\n--- 合并抗精神病药 ---\n")
if (length(ap_tbl) > 0) print(ap_tbl)

# 抗胆碱能药(可能加重误吸)
anticholinergics <- c("BENZTROPINE", "BIPERIDEN", "TRIHEXYPHENIDYL", "PROCYCLIDINE",
                       "HYOSCINE", "SCOPOLAMINE", "ATROPINE", "OXYBUTYNIN",
                       "TOLTERODINE", "SOLIFENACIN", "DARIFENACIN")
ac_pattern <- paste(anticholinergics, collapse = "|")
conmed_ac <- conmed_data[grepl(ac_pattern, prod_ai, ignore.case = TRUE)]
ac_tbl <- sort(table(conmed_ac$prod_ai), decreasing = TRUE)
cat("\n--- 合并抗胆碱能药 (可能影响唾液/误吸) ---\n")
if (length(ac_tbl) > 0) print(ac_tbl)

# ============================================================
# Part 6: 可视化
# ============================================================
message("\n=========================================")
message("Part 6: 生成可视化")
message("=========================================")

# 使用基础R绘图,不依赖ggplot2
output_prefix <- "clozapine_lung_fig"

# --- 图1: 信号检测森林图 ---
png(file.path(data_dir, paste0(output_prefix, "_forest.png")),
    width = 1000, height = 700, res = 120)

par(mar = c(5, 10, 4, 2))
sig_data <- signal_pt_df[order(ROR)]

# 取ROR在合理范围的PT
plot_data <- sig_data[ROR < 100]

plot_xlim <- c(0.1, max(c(plot_data$ROR_CI95_high, 30), na.rm = TRUE))

plot(1, type = "n", xlim = plot_xlim, ylim = c(0.5, nrow(plot_data) + 0.8),
     xlab = "ROR (95% CI)", ylab = "", yaxt = "n", log = "x",
     main = "氯氮平 vs 利培酮 肺部感染信号 (ROR Forest Plot)")

abline(v = 1, lty = 2, col = "grey50", lwd = 2)

for (i in 1:nrow(plot_data)) {
  ror_val <- plot_data$ROR[i]
  ci_low <- plot_data$ROR_CI95_low[i]
  ci_high <- plot_data$ROR_CI95_high[i]
  y_pos <- nrow(plot_data) - i + 1

  color <- ifelse(ci_low > 1, "#E41A1C", "#377EB8")
  points(ror_val, y_pos, pch = 16, col = color, cex = 1.2)
  segments(ci_low, y_pos, ci_high, y_pos, col = color, lwd = 2.5)

  label_text <- sprintf("%-45s (n=%d)", plot_data$Label[i], plot_data$a[i])
  axis(2, at = y_pos, labels = label_text, las = 2, cex.axis = 0.7,
       col.axis = color)

  # 显示数值
  text(ci_high * 1.15, y_pos,
       labels = sprintf("%.1f (%.1f-%.1f)", ror_val, ci_low, ci_high),
       cex = 0.55, adj = 0)
}

legend("bottomright",
       legend = c("Positive Signal (ROR lower CI > 1)",
                  "No Significant Signal"),
       col = c("#E41A1C", "#377EB8"), pch = 16, cex = 0.8)
dev.off()
message("  保存: *_forest.png")

# --- 图2: 年度趋势 ---
png(file.path(data_dir, paste0(output_prefix, "_yearly.png")),
    width = 1200, height = 600, res = 120)

par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

# 总报告趋势
yr_all <- trend_data[, .(Reports = uniqueN(primaryid)), by = event_year][order(event_year)]
barplot(yr_all$Reports, names.arg = yr_all$event_year,
        col = "#4575B4", border = "white",
        main = "氯氮平肺部感染 年度报告趋势",
        xlab = "年份", ylab = "报告数", las = 2)

# Top 6 PT 年度堆叠面积
top6 <- names(sort(table(cloz_lung_data$pt), decreasing = TRUE))[1:6]
yr_pt <- cloz_lung_data[pt %in% top6 & event_year >= 2018 & event_year <= 2025,
                         .(N = uniqueN(primaryid)), by = .(event_year, pt)]

# 转换为宽表
yr_pt_wide <- dcast(yr_pt, event_year ~ pt, value.var = "N", fill = 0)
yr_mat <- as.matrix(yr_pt_wide[, -1])
rownames(yr_mat) <- yr_pt_wide$event_year

colors_pt <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")
barplot(t(yr_mat), beside = FALSE, col = colors_pt, border = "white",
        main = "主要肺部感染类型年度分布",
        xlab = "年份", ylab = "报告数",
        legend.text = colnames(yr_mat),
        args.legend = list(x = "topleft", cex = 0.6, bty = "n"))
dev.off()
message("  保存: *_yearly.png")

# --- 图3: 年龄×性别分析 ---
png(file.path(data_dir, paste0(output_prefix, "_age_sex.png")),
    width = 1000, height = 600, res = 120)

par(mfrow = c(1, 2), mar = c(6, 5, 4, 2))

# 年龄分布直方图
hist(cloz_lung_data$age_num[cloz_lung_data$age_num > 0 & cloz_lung_data$age_num < 100],
     breaks = 30, col = "#4575B4", border = "white",
     main = "氯氮平肺部感染 年龄分布",
     xlab = "年龄 (岁)", ylab = "报告数")

abline(v = median(cloz_lung_data$age_num, na.rm = TRUE), col = "red", lwd = 2, lty = 2)
text(median(cloz_lung_data$age_num, na.rm = TRUE) + 5, par("usr")[4] * 0.9,
     sprintf("Median = %.1f", median(cloz_lung_data$age_num, na.rm = TRUE)), col = "red")

# 性别饼图
sex_tbl <- table(cloz_lung_data$sex, useNA = "ifany")
sex_labels <- c("Female", "Male", "Unknown")
sex_pct <- round(100 * sex_tbl / sum(sex_tbl), 1)
pie(sex_tbl,
    labels = sprintf("%s\n%d (%.1f%%)", sex_labels, sex_tbl, sex_pct),
    col = c("#F781BF", "#377EB8", "#BEBEBE"),
    main = "性别分布")
dev.off()
message("  保存: *_age_sex.png")

# --- 图4: 信号检测热力图矩阵 ---
png(file.path(data_dir, paste0(output_prefix, "_heatmap.png")),
    width = 1000, height = 800, res = 120)

par(mar = c(5, 15, 4, 2))

# 取Top信号
heat_data <- signal_pt_df[order(-IC)][1:min(25, nrow(signal_pt_df))]

# 绘制条形图展示IC值
bar_colors <- ifelse(heat_data$IC_025 > 0, "#E41A1C", "#377EB8")
y_positions <- barplot(heat_data$IC, horiz = TRUE,
                       col = bar_colors, border = "white",
                       xlab = "Information Component (IC)",
                       main = "氯氮平 vs 利培酮 肺部感染信号强度 (IC值)",
                       names.arg = rep("", nrow(heat_data)),
                       xlim = c(min(c(heat_data$IC_025, 0)) - 0.5,
                                max(c(heat_data$IC_975, 2)) + 0.5))

text(y = y_positions, x = par("usr")[1],
     labels = heat_data$Label, adj = 1, cex = 0.7, xpd = TRUE)

# 添加IC 95%CI 线段
for (i in 1:nrow(heat_data)) {
  segments(heat_data$IC_025[i], y_positions[i],
           heat_data$IC_975[i], y_positions[i],
           col = bar_colors[i], lwd = 2)
}

abline(v = 0, lty = 2, col = "grey50")

# 图例
legend("bottomright",
       legend = c(sprintf("IC025 > 0 (n=%d)", sum(heat_data$IC_025 > 0)),
                  sprintf("IC025 <= 0 (n=%d)", sum(heat_data$IC_025 <= 0))),
       fill = c("#E41A1C", "#377EB8"), cex = 0.8, bty = "n")
dev.off()
message("  保存: *_heatmap.png")

# ============================================================
# Part 7: 导出汇总报告
# ============================================================

# 生成汇总数据集
cat("\n============================================================\n")
cat("              生成汇总报告\n")
cat("============================================================\n")

# 核心信号检测结果导出
signal_export <- signal_pt_df[, .(
  PT_Term = Label,
  Clozapine_Event = a,
  Clozapine_NoEvent = b,
  Risperidone_Event = c,
  Risperidone_NoEvent = d,
  ROR, ROR_CI95_low, ROR_CI95_high,
  PRR, PRR_CI95_low, PRR_CI95_high,
  Chi_sq, IC, IC_025, IC_975,
  P_value, Signal
)]

fwrite(signal_export, file.path(data_dir, "clozapine_lung_signal_detection.csv"), bom = TRUE)

# 年度趋势导出
fwrite(yearly, file.path(data_dir, "clozapine_lung_yearly_trend.csv"), bom = TRUE)

# Excel 完整报告
if (requireNamespace("openxlsx", quietly = TRUE)) {
  library(openxlsx)
  output_xlsx <- file.path(data_dir, "clozapine_lung_infection_advanced_report.xlsx")

  wb <- createWorkbook()

  addWorksheet(wb, "信号检测")
  writeData(wb, "信号检测", signal_export)

  addWorksheet(wb, "年度趋势")
  writeData(wb, "年度趋势", yearly)

  addWorksheet(wb, "季度趋势")
  writeData(wb, "季度趋势", qtrly)

  addWorksheet(wb, "伴随药物 Top50")
  conmed_df <- data.frame(
    药物名称 = names(conmed_top)[1:min(50, length(conmed_top))],
    报告数 = as.integer(conmed_top[1:min(50, length(conmed_top))])
  )
  writeData(wb, "伴随药物 Top50", conmed_df)

  addWorksheet(wb, "年龄性别交叉表")
  age_sex_wide <- dcast(age_sex_tbl, age_group ~ sex, value.var = "N", fill = 0)
  writeData(wb, "年龄性别交叉表", age_sex_wide)

  addWorksheet(wb, "吸入性肺炎")
  if (nrow(asp_data) > 0) {
    asp_summary <- data.frame(
      Metric = c("总记录", "独特报告", "独特病例", "男性%", "女性%"),
      Value = c(nrow(asp_data), uniqueN(asp_data$primaryid),
                uniqueN(asp_data$caseid),
                sprintf("%.1f%%", 100 * sum(asp_data$sex == "M", na.rm = TRUE) / nrow(asp_data)),
                sprintf("%.1f%%", 100 * sum(asp_data$sex == "F", na.rm = TRUE) / nrow(asp_data)))
    )
    writeData(wb, "吸入性肺炎", asp_summary)
  }

  saveWorkbook(wb, output_xlsx, overwrite = TRUE)
  message(sprintf("Excel报告: %s", output_xlsx))
}

# ============================================================
# 最终汇总
# ============================================================
cat("\n============================================================\n")
cat("           氯氮平肺部感染 深度分析完成\n")
cat("============================================================\n")
cat(sprintf("信号检测: %d 个 PT 术语参与分析\n", nrow(signal_pt_df)))
cat(sprintf("阳性信号: %d 个 (ROR>1 & ROR_CI_low>1 & n>=3 & IC025>0)\n",
            sum(signal_pt_df$Signal == "Positive Signal")))

cat("\n--- 核心发现 ---\n")

# 最强的阳性信号
pos_signals <- signal_pt_df[Signal == "Positive Signal"][order(-IC)]
if (nrow(pos_signals) > 0) {
  cat("阳性信号 (按IC降序):\n")
  for (i in 1:nrow(pos_signals)) {
    cat(sprintf("  %s: ROR=%.1f (%.1f-%.1f), IC=%.2f (%.2f-%.2f), n=%d\n",
                pos_signals$Label[i],
                pos_signals$ROR[i], pos_signals$ROR_CI95_low[i], pos_signals$ROR_CI95_high[i],
                pos_signals$IC[i], pos_signals$IC_025[i], pos_signals$IC_975[i],
                pos_signals$a[i]))
  }
}

cat(sprintf("\n年度趋势: %d → %d (2018 → 2025)\n",
            yr_all[event_year == 2018, Reports],
            yr_all[event_year == 2025, Reports]))
cat(sprintf("吸入性肺炎: %d 报告, 占肺部感染的 %.1f%%\n",
            uniqueN(asp_data$primaryid),
            100 * uniqueN(asp_data$primaryid) / uniqueN(cloz_lung_data$primaryid)))
cat(sprintf("报告主要来源: 英国(%d), 加拿大(%d), 美国(%d)\n",
            sum(cloz_lung_data$reporter_country == "GB"),
            sum(cloz_lung_data$reporter_country == "CA"),
            sum(cloz_lung_data$reporter_country == "US")))

cat("\n导出文件:\n")
cat("  clozapine_lung_signal_detection.csv  — 信号检测结果\n")
cat("  clozapine_lung_yearly_trend.csv      — 年度趋势\n")
cat("  clozapine_lung_infection_advanced_report.xlsx — 完整报告(多sheet)\n")
cat("  clozapine_lung_fig_*.png             — 可视化图表\n")
cat("============================================================\n")
