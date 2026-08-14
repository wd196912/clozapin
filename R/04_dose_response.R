# ============================================================
# 氯氮平肺部感染 — 剂量-反应关系深度分析
# 包含: RCS非线性建模、多变量Logistic回归、亚组分析、可视化
# ============================================================

library(data.table)
library(dplyr)

data_dir <- "F:/faersdata"
setwd(data_dir)

# ============================================================
# Part 0: 数据导入 (复用已有数据管道)
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

# Fatal outcome = outc_cod == "DE" in quarterly OUTC files (manuscript section 2.5)
outc_files <- list.files("faers_ascii_2025q4/ASCII", pattern = "^OUTC.*\\.txt$",
                         full.names = TRUE)
outc <- rbindlist(lapply(outc_files, function(f) {
  fread(f, sep = "$", header = FALSE, select = 1:3,
        col.names = c("primaryid", "caseid", "outc_cod"),
        na.strings = "", blank.lines.skip = TRUE, strip.white = TRUE)
}))
death_ids <- unique(outc[outc_cod == "DE", primaryid])

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

# 肺部感染术语
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

# 构建氯氮平肺部感染完整数据
drug_cloz_ps <- drug_raw[drug_type == "Clozapine" & role_cod == "PS"]

cloz_lung <- reac_raw[is_lung == TRUE] %>%
  inner_join(drug_cloz_ps, by = c("primaryid", "caseid"),
             relationship = "many-to-many") %>%
  inner_join(demo_raw, by = c("primaryid", "caseid"),
             relationship = "many-to-many") %>%
  as.data.table()

message(sprintf("氯氮平肺部感染: %d 行, %d 独特primaryid",
                nrow(cloz_lung), uniqueN(cloz_lung$primaryid)))

# ============================================================
# Part 1: 剂量数据深度清洗
# ============================================================
message("\n=========================================")
message("Part 1: 剂量数据深度清洗")
message("=========================================")

dose_raw <- copy(cloz_lung)

# 1.1 提取数值剂量
dose_raw[, dose_amt_num := as.numeric(dose_amt)]

# 1.2 标准化剂量单位 → mg
dose_raw[, dose_unit_clean := toupper(gsub("[^A-Z]", "", dose_unit))]
dose_raw[dose_unit_clean %in% c("MG","MGS","MILLIGRAMS","MILLIGRAM"), dose_mg := dose_amt_num]
dose_raw[dose_unit_clean %in% c("G","GM","GRAMS","GRAM"), dose_mg := dose_amt_num * 1000]
dose_raw[dose_unit_clean %in% c("UG","MCG","MICROGRAMS"), dose_mg := dose_amt_num / 1000]
dose_raw[dose_unit_clean %in% c("ML","MILLILITERS","MILLILITER"), dose_mg := dose_amt_num * 100]  # 口服液: 100mg/mL
dose_raw[dose_unit_clean == "NG", dose_mg := dose_amt_num / 1e6]
dose_raw[dose_unit_clean %in% c("IU","UNITS","UNIT"), dose_mg := NA]  # 无法转换

# 1.3 清洗剂量频率 → 日剂量
dose_raw[, dose_freq_clean := toupper(gsub("[^A-Z0-9]", "", dose_freq))]
dose_raw[, daily_mult := 1]  # 默认每日1次

# 识别频率模式
dose_raw[grepl("BID|BID$|TWICE.*DAILY|TWODAILY|2XD|2X.*DAILY|TWICEADAY|EVERY12H|Q12H|BIDLY",
               dose_freq_clean), daily_mult := 2]
dose_raw[grepl("TID|TID$|THREE.*DAILY|3XD|3X.*DAILY|THREETIMESDAILY|EVERY8H|Q8H|TIDLY",
               dose_freq_clean), daily_mult := 3]
dose_raw[grepl("QID|QID$|FOUR.*DAILY|4XD|4X.*DAILY|FOURTIMESDAILY|EVERY6H|Q6H|QIDLY",
               dose_freq_clean), daily_mult := 4]
dose_raw[grepl("QD|ONCE.*DAILY|1XD|1X.*DAILY|ONCEDAILY|ONCEADAY|DAILY|PERDAY|EVERYDAY|EVERY24H|Q24H|QDAM|QDPM|QDAY|HS|BEDTIME",
               dose_freq_clean), daily_mult := 1]
dose_raw[grepl("QOD|EVERYOTHER|ALTERNATE|EVERY2D|EVERY48H|Q48H",
               dose_freq_clean), daily_mult := 0.5]
dose_raw[grepl("QW|WEEKLY|ONCEAWEEK|EVERYWEEK|Q7D|EVERY7D",
               dose_freq_clean), daily_mult := 1/7]
dose_raw[grepl("BIW|TWICEWEEKLY|TWICEAWEEK|2XWEEKLY",
               dose_freq_clean), daily_mult := 2/7]
dose_raw[grepl("TIW|THREETIMESWEEKLY|3XWEEKLY",
               dose_freq_clean), daily_mult := 3/7]
dose_raw[grepl("QHS|ATBEDTIME|ATNIGHT|NOCTE|NOCTURNAL",
               dose_freq_clean), daily_mult := 1]

# 1.4 计算日剂量
dose_raw[, daily_dose_mg := dose_mg * daily_mult]

# 1.5 按 primaryid 取中位日剂量（去重）
dose_by_id <- dose_raw[!is.na(daily_dose_mg) & daily_dose_mg > 0 & daily_dose_mg <= 2000,
                        .(daily_dose = median(daily_dose_mg, na.rm = TRUE),
                          dose_amt = first(dose_amt),
                          dose_unit = first(dose_unit),
                          dose_freq = first(dose_freq),
                          route = first(route),
                          dechal = first(dechal),
                          rechal = first(rechal),
                          cum_dose_chr = first(cum_dose_chr),
                          cum_dose_unit = first(cum_dose_unit)),
                        by = primaryid]

message(sprintf("有效剂量 primaryid: %d / %d (%.1f%%)",
                nrow(dose_by_id), uniqueN(dose_raw$primaryid),
                100 * nrow(dose_by_id) / uniqueN(dose_raw$primaryid)))

# 1.6 合并回人口学数据
demo_unique <- cloz_lung[, .(
  age = first(age), age_cod = first(age_cod), sex = first(sex),
  i_f_code = first(i_f_code), reporter_country = first(reporter_country),
  event_dt = first(event_dt), wt = first(wt), wt_cod = first(wt_cod),
  e_sub = first(e_sub), rept_cod = first(rept_cod)
), by = primaryid]

# PT 信息（每个 primaryid 可能有多个 PT）
pt_by_id <- cloz_lung[, .(
  pts = paste(sort(unique(pt)), collapse = "; "),
  pt_primary = first(pt),
  is_aspiration = any(grepl("aspiration|aspiratio", pt, ignore.case = TRUE)),
  is_pneumonia = any(pt == "Pneumonia"),
  is_LRTI = any(grepl("lower respiratory", pt, ignore.case = TRUE))
), by = primaryid]

dose_full <- dose_by_id %>%
  inner_join(demo_unique, by = "primaryid") %>%
  inner_join(pt_by_id, by = "primaryid") %>%
  as.data.table()

# 1.7 计算年龄
dose_full[, age_num := as.numeric(age)]
dose_full[age_cod %in% c("DEC"), age_num := age_num / 12]
dose_full[age_cod %in% c("WK","WEEK"), age_num := age_num / 52]
dose_full[age_cod %in% c("DY","DAY"), age_num := age_num / 365]

dose_full[, age_group := cut(age_num,
  breaks = c(-Inf, 40, 50, 60, 70, Inf),
  labels = c("<40","40-49","50-59","60-69","70+"))]

dose_full[, sex_bin := ifelse(sex == "M", 1, ifelse(sex == "F", 0, NA))]

# 年份
dose_full[, event_year := as.integer(substr(event_dt, 1, 4))]

# 结局变量
dose_full[, fatal := ifelse(primaryid %in% death_ids, 1, 0)]

message(sprintf("最终分析数据集: %d 条记录", nrow(dose_full)))

cat("\n--- 日剂量分布 ---\n")
cat(sprintf("  Mean (SD): %.1f (%.1f) mg/d\n", mean(dose_full$daily_dose), sd(dose_full$daily_dose)))
cat(sprintf("  Median (IQR): %.1f (%.1f-%.1f) mg/d\n",
            median(dose_full$daily_dose),
            quantile(dose_full$daily_dose, 0.25),
            quantile(dose_full$daily_dose, 0.75)))
cat(sprintf("  Range: %.0f-%.0f mg/d\n", min(dose_full$daily_dose), max(dose_full$daily_dose)))

# 剂量分组
dose_full[, dose_group := cut(daily_dose,
  breaks = c(0, 100, 200, 300, 400, 500, 600, 900, 2001),
  labels = c("<100","100-199","200-299","300-399","400-499","500-599","600-899","900+"),
  include.lowest = TRUE)]

# ============================================================
# Part 2: Logistic 回归 — 剂量-反应关系
# ============================================================
message("\n=========================================")
message("Part 2: Logistic 回归分析")
message("=========================================")

# 2.1 连续剂量对吸入性肺炎的影响
cat("\n--- 模型1: 连续剂量 → 吸入性肺炎 (Crude) ---\n")
m1 <- glm(is_aspiration ~ daily_dose, data = dose_full, family = binomial)
print(summary(m1))
or_per_100mg <- exp(coef(m1)["daily_dose"] * 100)
cat(sprintf("OR per 100mg increase = %.3f (95%%CI: %.3f-%.3f)\n",
            or_per_100mg,
            exp(confint(m1)["daily_dose", 1] * 100),
            exp(confint(m1)["daily_dose", 2] * 100)))

# 2.2 调整年龄和性别
cat("\n--- 模型2: 连续剂量 → 吸入性肺炎 (Adjusted: age + sex) ---\n")
m2 <- glm(is_aspiration ~ daily_dose + age_num + sex_bin,
          data = dose_full, family = binomial)
print(summary(m2))
or_adj_100 <- exp(coef(m2)["daily_dose"] * 100)
cat(sprintf("Adjusted OR per 100mg = %.3f (95%%CI: %.3f-%.3f)\n",
            or_adj_100,
            exp(confint(m2)["daily_dose", 1] * 100),
            exp(confint(m2)["daily_dose", 2] * 100)))

# 2.3 连续剂量对致死结局
cat("\n--- 模型3: 连续剂量 → 致死结局 (Adjusted) ---\n")
m3 <- glm(fatal ~ daily_dose + age_num + sex_bin,
          data = dose_full, family = binomial)
print(summary(m3))
or_fatal_100 <- exp(coef(m3)["daily_dose"] * 100)
cat(sprintf("Adjusted OR (fatal) per 100mg = %.3f (95%%CI: %.3f-%.3f)\n",
            or_fatal_100,
            exp(confint(m3)["daily_dose", 1] * 100),
            exp(confint(m3)["daily_dose", 2] * 100)))

# 2.4 连续剂量对肺炎(总)
cat("\n--- 模型4: 连续剂量 → Pneumonia (Adjusted) ---\n")
m4 <- glm(is_pneumonia ~ daily_dose + age_num + sex_bin,
          data = dose_full, family = binomial)
print(summary(m4))
or_pneu_100 <- exp(coef(m4)["daily_dose"] * 100)
cat(sprintf("Adjusted OR (pneumonia) per 100mg = %.3f (95%%CI: %.3f-%.3f)\n",
            or_pneu_100,
            exp(confint(m4)["daily_dose", 1] * 100),
            exp(confint(m4)["daily_dose", 2] * 100)))

# 2.5 剂量分组 Logistic（趋势检验的替代）
cat("\n--- 模型5: 剂量分组 → 吸入性肺炎 (趋势检验) ---\n")
dose_full[, dose_group_ord := as.numeric(dose_group)]
m5 <- glm(is_aspiration ~ dose_group, data = dose_full, family = binomial)
print(summary(m5))

# 线性趋势检验（将剂量组作为连续变量）
m5_trend <- glm(is_aspiration ~ dose_group_ord, data = dose_full, family = binomial)
trend_p <- coef(summary(m5_trend))["dose_group_ord", "Pr(>|z|)"]
cat(sprintf("\n剂量分组线性趋势检验 (吸入性肺炎): p = %.6f\n", trend_p))

# 2.6 调整后的剂量分组模型
cat("\n--- 模型6: 剂量分组 → 吸入性肺炎 (Adjusted) ---\n")
m6 <- glm(is_aspiration ~ dose_group + age_num + sex_bin, data = dose_full, family = binomial)
print(summary(m6))

# 提取各剂量组的 adjusted OR (vs 最低剂量组)
cat("\n各剂量组 Adjusted OR (ref=<100mg):\n")
dose_grp_levels <- levels(dose_full$dose_group)
# 使用 confint.default 避免 profile CI 在分离数据时的收敛问题
m6_coef <- coef(m6)
m6_se <- sqrt(diag(vcov(m6)))
dose_coef_names <- grep("dose_group", names(m6_coef), value = TRUE)

or_list <- list()
or_list[[1]] <- data.frame(
  Dose_Group = "<100",
  OR = 1, CI_low = 1, CI_high = 1, P_value = NA, stringsAsFactors = FALSE)

for (i in seq_along(dose_grp_levels)[-1]) {
  coef_name <- paste0("dose_group", dose_grp_levels[i])
  if (coef_name %in% names(m6_coef)) {
    beta <- m6_coef[coef_name]
    se <- m6_se[coef_name]
    z <- beta / se
    p <- 2 * pnorm(-abs(z))
    or_list[[length(or_list) + 1]] <- data.frame(
      Dose_Group = dose_grp_levels[i],
      OR = exp(beta),
      CI_low = exp(beta - 1.96 * se),
      CI_high = exp(beta + 1.96 * se),
      P_value = p,
      stringsAsFactors = FALSE)
  }
}
or_m6 <- do.call(rbind, or_list)
print(or_m6, row.names = FALSE)

# 报告最高剂量组 vs 最低
cat(sprintf("\n  900+ mg vs <100mg: OR=%.1f, p=%.4f\n",
            or_m6$OR[nrow(or_m6)], or_m6$P_value[nrow(or_m6)]))

# ============================================================
# Part 3: 限制性立方样条 (RCS) — 非线性剂量-反应
# ============================================================
message("\n=========================================")
message("Part 3: 限制性立方样条 (RCS) 非线性建模")
message("=========================================")

# 3.1 手动实现 RCS 基础函数
rcs_basis <- function(x, knots) {
  # 限制性立方样条基函数: k knots → k-1 basis functions
  k <- length(knots)
  if (k < 2) stop("Need at least 2 knots")

  n <- length(x)
  ncol_out <- k - 1  # k-1 个基函数
  X <- matrix(0, nrow = n, ncol = ncol_out)

  # 第一个基: 线性项
  X[, 1] <- x

  # 非线性基 (k-2 个)
  if (k >= 3) {
    for (j in 1:(k-2)) {
      tj <- knots[j+1]
      X[, j+1] <- pmax(x - tj, 0)^3 -
        ((knots[k-1] - tj) / (knots[k] - knots[k-1])) * pmax(x - knots[k-1], 0)^3 +
        ((knots[k] - tj) / (knots[k] - knots[k-1])) * pmax(x - knots[k], 0)^3
    }
  }

  colnames(X) <- paste0("rcs", 1:ncol_out)
  X
}

# 3.2 选择 knot 位置（在剂量分布的 10th, 50th, 90th 百分位）
knots_3 <- quantile(dose_full$daily_dose, probs = c(0.10, 0.50, 0.90), na.rm = TRUE)
cat(sprintf("RCS Knots (10th, 50th, 90th percentile): %.0f, %.0f, %.0f mg/d\n",
            knots_3[1], knots_3[2], knots_3[3]))

# 创建 RCS 基函数
rcs_mat <- rcs_basis(dose_full$daily_dose, knots_3)
dose_full[, rcs1 := rcs_mat[, 1]]
dose_full[, rcs2 := rcs_mat[, 2]]

# 3.3 RCS Logistic 模型 — 吸入性肺炎
cat("\n--- RCS模型1: 吸入性肺炎 (Unadjusted) ---\n")
m_rcs1 <- glm(is_aspiration ~ rcs1 + rcs2, data = dose_full, family = binomial)

# 计算预测概率（在剂量范围内）
dose_range <- seq(25, 900, by = 5)
rcs_pred <- rcs_basis(dose_range, knots_3)
X_pred <- cbind(1, rcs_pred)
eta <- X_pred %*% coef(m_rcs1)
pred_prob <- 1 / (1 + exp(-eta))
se_eta <- sqrt(diag(X_pred %*% vcov(m_rcs1) %*% t(X_pred)))
eta_low <- eta - 1.96 * se_eta
eta_high <- eta + 1.96 * se_eta
pred_low <- 1 / (1 + exp(-eta_low))
pred_high <- 1 / (1 + exp(-eta_high))

rcs_asp_result <- data.frame(
  dose = dose_range,
  prob = pred_prob,
  lower = pred_low,
  upper = pred_high
)

cat("RCS模型 (吸入性肺炎) — 非线性检验:\n")
# 嵌套模型比较: RCS vs 纯线性
m_rcs_linear <- glm(is_aspiration ~ daily_dose, data = dose_full, family = binomial)
lrt_nonlinear <- anova(m_rcs_linear, m_rcs1, test = "Chisq")
print(lrt_nonlinear)
cat(sprintf("Nonlinearity test p-value: %s\n",
            format.pval(lrt_nonlinear$`Pr(>Chi)`[2], digits = 4)))

# 3.4 RCS Logistic 模型 — 致死结局 (Adjusted)
cat("\n--- RCS模型2: 致死结局 (Adjusted: age+sex) ---\n")
m_rcs2 <- glm(fatal ~ rcs1 + rcs2 + age_num + sex_bin,
              data = dose_full, family = binomial)

# 预测 (固定协变量为均值/参考)
age_med <- median(dose_full$age_num, na.rm = TRUE)
X_pred2 <- cbind(1, rcs_pred, age_med, 0)  # 男性参考
eta2 <- X_pred2 %*% coef(m_rcs2)
pred_prob2 <- 1 / (1 + exp(-eta2))
se_eta2 <- sqrt(diag(X_pred2 %*% vcov(m_rcs2) %*% t(X_pred2)))
eta2_low <- eta2 - 1.96 * se_eta2
eta2_high <- eta2 + 1.96 * se_eta2

rcs_fatal_result <- data.frame(
  dose = dose_range,
  prob = pred_prob2,
  lower = 1 / (1 + exp(-eta2_low)),
  upper = 1 / (1 + exp(-eta2_high))
)

# 3.5 RCS Logistic — 肺炎
cat("\n--- RCS模型3: Pneumonia (Adjusted) ---\n")
m_rcs3 <- glm(is_pneumonia ~ rcs1 + rcs2 + age_num + sex_bin,
              data = dose_full, family = binomial)

eta3 <- cbind(1, rcs_pred, age_med, 0) %*% coef(m_rcs3)
pred_prob3 <- 1 / (1 + exp(-eta3))
se3 <- sqrt(diag(cbind(1, rcs_pred, age_med, 0) %*% vcov(m_rcs3) %*%
                 t(cbind(1, rcs_pred, age_med, 0))))
rcs_pneu_result <- data.frame(
  dose = dose_range,
  prob = pred_prob3,
  lower = 1 / (1 + exp(-eta3 - 1.96 * se3)),
  upper = 1 / (1 + exp(-eta3 + 1.96 * se3))
)

m_rcs_linear3 <- glm(is_pneumonia ~ daily_dose + age_num + sex_bin,
                     data = dose_full, family = binomial)
lrt_nonlinear3 <- anova(m_rcs_linear3, m_rcs3, test = "Chisq")
cat(sprintf("Pneumonia nonlinearity p: %s\n",
            format.pval(lrt_nonlinear3$`Pr(>Chi)`[2], digits = 4)))

# ============================================================
# Part 4: 亚组分析 — 剂量-反应异质性
# ============================================================
message("\n=========================================")
message("Part 4: 亚组剂量-反应分析")
message("=========================================")

# 4.1 按年龄组分层: 剂量→吸入性肺炎
cat("\n--- 年龄亚组: 剂量→吸入性肺炎 ---\n")
age_groups <- c("<40","40-49","50-59","60-69","70+")
subgroup_or <- data.frame()

for (ag in age_groups) {
  sub <- dose_full[age_group == ag]
  if (nrow(sub) >= 20) {
    m_sub <- glm(is_aspiration ~ daily_dose, data = sub, family = binomial)
    ci <- confint(m_sub)["daily_dose", ]
    subgroup_or <- rbind(subgroup_or, data.frame(
      Subgroup = ag,
      N = nrow(sub),
      OR_per_100mg = exp(coef(m_sub)["daily_dose"] * 100),
      CI_low = exp(ci[1] * 100),
      CI_high = exp(ci[2] * 100),
      P_value = coef(summary(m_sub))["daily_dose", "Pr(>|z|)"]
    ))
  }
}
print(subgroup_or)

# 4.2 按性别分层
cat("\n--- 性别亚组: 剂量→吸入性肺炎 ---\n")
for (sx in c("F", "M")) {
  sub <- dose_full[sex == sx]
  m_sub <- glm(is_aspiration ~ daily_dose + age_num, data = sub, family = binomial)
  cat(sprintf("\n%s (n=%d):\n", ifelse(sx=="M","男性","女性"), nrow(sub)))
  or_val <- exp(coef(m_sub)["daily_dose"] * 100)
  ci <- exp(confint(m_sub)["daily_dose", ] * 100)
  cat(sprintf("  Adjusted OR per 100mg = %.3f (%.3f-%.3f), p=%.4f\n",
              or_val, ci[1], ci[2],
              coef(summary(m_sub))["daily_dose", "Pr(>|z|)"]))
}

# 4.3 交互作用检验: 剂量×年龄
cat("\n--- 交互作用: 剂量 × 年龄 (吸入性肺炎) ---\n")
m_interact <- glm(is_aspiration ~ daily_dose * age_num + sex_bin,
                  data = dose_full, family = binomial)
print(summary(m_interact))
cat(sprintf("Dose × Age interaction p = %.6f\n",
            coef(summary(m_interact))["daily_dose:age_num", "Pr(>|z|)"]))

# 4.4 交互作用: 剂量×性别
cat("\n--- 交互作用: 剂量 × 性别 (吸入性肺炎) ---\n")
dose_full[, dose_c100 := daily_dose / 100]
m_sex_interact <- glm(is_aspiration ~ dose_c100 * sex_bin + age_num,
                      data = dose_full, family = binomial)
print(summary(m_sex_interact))
cat(sprintf("Dose × Sex interaction p = %.6f\n",
            coef(summary(m_sex_interact))["dose_c100:sex_bin", "Pr(>|z|)"]))

# ============================================================
# Part 5: 剂量-频率联合分析
# ============================================================
message("\n=========================================")
message("Part 5: 剂量-频率联合分析")
message("=========================================")

# 5.1 给药频率分布
cat("\n--- 给药频率分布 ---\n")
freq_tbl <- sort(table(dose_full$dose_freq), decreasing = TRUE)
freq_pct <- round(100 * freq_tbl / sum(freq_tbl), 1)
freq_df <- data.frame(Frequency = names(freq_tbl), N = as.integer(freq_tbl), Pct = freq_pct)
print(head(freq_df, 15))

# 5.2 单次剂量 vs 分次给药
dose_full[, dosing_schedule := ifelse(
  grepl("BID|TID|QID|TWO|THREE|FOUR|TWICE|TID|QID", dose_freq, ignore.case = TRUE),
  "Multiple daily", "Once daily")]
dose_full[is.na(dose_freq) | dose_freq == "", dosing_schedule := "Unknown"]

cat("\n--- 给药方案 vs 吸入性肺炎 ---\n")
sched_tbl <- dose_full[, .(
  N = .N,
  Aspiration = sum(is_aspiration),
  Asp_pct = round(100 * sum(is_aspiration) / .N, 1),
  Mean_Dose = mean(daily_dose, na.rm = TRUE),
  Fatal_pct = round(100 * sum(fatal) / .N, 1)
), by = dosing_schedule]
print(sched_tbl)

# 5.3 单次大剂量 vs 多次分服的效应差异
cat("\n--- 分层: 给药方案 × 剂量水平 ---\n")
dose_full[, dose_level := ifelse(daily_dose >= 400, "High (>=400mg)", "Low (<400mg)")]
sched_dose <- dose_full[, .(
  N = .N,
  Asp_pct = round(100 * sum(is_aspiration) / .N, 1),
  Fatal_pct = round(100 * sum(fatal) / .N, 1)
), by = .(dose_level, dosing_schedule)][order(dose_level, dosing_schedule)]
print(sched_dose)

# ============================================================
# Part 6: 剂量-反应的绝对风险估计
# ============================================================
message("\n=========================================")
message("Part 6: 绝对风险与剂量阈值分析")
message("=========================================")

# 6.1 寻找吸入性肺炎风险的剂量阈值 (分段回归)
cat("\n--- 剂量阈值探索 (吸入性肺炎) ---\n")
thresholds <- c(100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600)
threshold_results <- data.frame()

for (th in thresholds) {
  dose_full[, above_th := daily_dose >= th]
  tbl <- table(dose_full$above_th, dose_full$is_aspiration)
  if (nrow(tbl) == 2 && ncol(tbl) == 2) {
    ft <- fisher.test(tbl)
    threshold_results <- rbind(threshold_results, data.frame(
      Threshold = th,
      OR = ft$estimate,
      CI_low = ft$conf.int[1],
      CI_high = ft$conf.int[2],
      P_value = ft$p.value,
      Above_N = sum(dose_full$above_th),
      Asp_Above = sum(dose_full$above_th & dose_full$is_aspiration),
      Asp_Below = sum(!dose_full$above_th & dose_full$is_aspiration)
    ))
  }
}
print(threshold_results[order(threshold_results$OR, decreasing = TRUE), ])

# 6.2 最佳阈值 (OR最大且有显著性)
best_th <- threshold_results[which.max(threshold_results$OR), ]
cat(sprintf("\n最佳阈值: %.0f mg/d, OR=%.2f (%.2f-%.2f), p=%.6f\n",
            best_th$Threshold, best_th$OR, best_th$CI_low, best_th$CI_high, best_th$P_value))

# 6.3 各剂量组的绝对风险 (NNH近似)
cat("\n--- 各剂量组吸入性肺炎绝对风险 ---\n")
abs_risk <- dose_full[, .(
  Total = .N,
  Aspiration_N = sum(is_aspiration),
  Asp_Risk = round(100 * sum(is_aspiration) / .N, 2),
  Fatal_N = sum(fatal),
  Fatal_Risk = round(100 * sum(fatal) / .N, 2)
), by = dose_group][order(dose_group)]
abs_risk[, NNH_asp := round(100 / (Asp_Risk - abs_risk$Asp_Risk[1]), 1)]
print(abs_risk)

# ============================================================
# Part 7: 死亡病例的剂量特征
# ============================================================
message("\n=========================================")
message("Part 7: 死亡病例剂量特征")
message("=========================================")

cat("\n--- 死亡 vs 存活 剂量对比 ---\n")
fatal_dose <- dose_full[fatal == 1]
nonfatal_dose <- dose_full[fatal == 0]

cat(sprintf("死亡组: n=%d, Median dose=%.0f (IQR: %.0f-%.0f)\n",
            nrow(fatal_dose),
            median(fatal_dose$daily_dose),
            quantile(fatal_dose$daily_dose, 0.25),
            quantile(fatal_dose$daily_dose, 0.75)))
cat(sprintf("存活组: n=%d, Median dose=%.0f (IQR: %.0f-%.0f)\n",
            nrow(nonfatal_dose),
            median(nonfatal_dose$daily_dose),
            quantile(nonfatal_dose$daily_dose, 0.25),
            quantile(nonfatal_dose$daily_dose, 0.75)))

# Mann-Whitney
mw_test <- wilcox.test(fatal_dose$daily_dose, nonfatal_dose$daily_dose)
cat(sprintf("Mann-Whitney p = %.6f\n", mw_test$p.value))

# 死亡病例的剂量分布
cat("\n--- 死亡病例剂量分组 ---\n")
fatal_by_dose <- fatal_dose[, .(
  Deaths = .N,
  Asp_Deaths = sum(is_aspiration),
  Asp_pct = round(100 * sum(is_aspiration) / .N, 1)
), by = dose_group][order(dose_group)]
print(fatal_by_dose)

# ============================================================
# Part 8: 综合可视化
# ============================================================
message("\n=========================================")
message("Part 8: 剂量-反应可视化")
message("=========================================")

output_dir <- "F:/faersdata"
pfx <- "clozapine_dose_response"

# --- 图1: RCS 剂量-反应曲线 (三合一) ---
png(file.path(output_dir, paste0(pfx, "_rcs_curves.png")),
    width = 1400, height = 500, res = 130)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

# 1A: 吸入性肺炎
plot(rcs_asp_result$dose, rcs_asp_result$prob * 100,
     type = "l", lwd = 3, col = "#E41A1C",
     xlim = c(0, 900), ylim = c(0, max(rcs_asp_result$upper) * 100 * 1.2),
     main = "吸入性肺炎剂量-反应曲线",
     xlab = "氯氮平日剂量 (mg)", ylab = "预测概率 (%)")
polygon(c(rcs_asp_result$dose, rev(rcs_asp_result$dose)),
        c(rcs_asp_result$lower * 100, rev(rcs_asp_result$upper * 100)),
        col = rgb(0.89, 0.10, 0.11, 0.2), border = NA)
lines(rcs_asp_result$dose, rcs_asp_result$prob * 100, lwd = 3, col = "#E41A1C")
rug(dose_full$daily_dose[dose_full$is_aspiration == 1], col = rgb(0,0,0,0.3), ticksize = 0.02)
abline(v = c(300, 600), lty = c(2,2), col = c("orange","darkred"), lwd = 1.5)
text(300, par("usr")[4] * 0.95, "300mg", col = "orange", cex = 0.7)
text(600, par("usr")[4] * 0.95, "600mg", col = "darkred", cex = 0.7)

# 1B: 致死结局
plot(rcs_fatal_result$dose, rcs_fatal_result$prob * 100,
     type = "l", lwd = 3, col = "#D73027",
     xlim = c(0, 900), ylim = c(0, max(rcs_fatal_result$upper) * 100 * 1.2),
     main = "致死结局剂量-反应曲线\n(Adjusted: age + sex)",
     xlab = "氯氮平日剂量 (mg)", ylab = "预测致死概率 (%)")
polygon(c(rcs_fatal_result$dose, rev(rcs_fatal_result$dose)),
        c(rcs_fatal_result$lower * 100, rev(rcs_fatal_result$upper * 100)),
        col = rgb(0.84, 0.19, 0.15, 0.2), border = NA)
lines(rcs_fatal_result$dose, rcs_fatal_result$prob * 100, lwd = 3, col = "#D73027")
abline(v = c(300, 600), lty = c(2,2), col = c("orange","darkred"), lwd = 1.5)

# 1C: 肺炎(总)
plot(rcs_pneu_result$dose, rcs_pneu_result$prob * 100,
     type = "l", lwd = 3, col = "#4575B4",
     xlim = c(0, 900), ylim = c(0, max(rcs_pneu_result$upper) * 100 * 1.2),
     main = "肺炎(总)剂量-反应曲线\n(Adjusted: age + sex)",
     xlab = "氯氮平日剂量 (mg)", ylab = "预测概率 (%)")
polygon(c(rcs_pneu_result$dose, rev(rcs_pneu_result$dose)),
        c(rcs_pneu_result$lower * 100, rev(rcs_pneu_result$upper * 100)),
        col = rgb(0.27, 0.46, 0.71, 0.2), border = NA)
lines(rcs_pneu_result$dose, rcs_pneu_result$prob * 100, lwd = 3, col = "#4575B4")
abline(v = c(300, 600), lty = c(2,2), col = c("orange","darkred"), lwd = 1.5)

dev.off()
message("  -> *_rcs_curves.png")

# --- 图2: 亚组剂量-反应森林图 ---
png(file.path(output_dir, paste0(pfx, "_subgroup_forest.png")),
    width = 1000, height = 600, res = 130)

par(mar = c(5, 12, 4, 2))
# 构建总表
main_or <- exp(coef(m1)["daily_dose"] * 100)
main_ci <- exp(confint(m1)["daily_dose", ] * 100)
main_adj_or <- exp(coef(m2)["daily_dose"] * 100)
main_adj_ci <- exp(confint(m2)["daily_dose", ] * 100)

forest_data <- rbind(
  data.frame(Label = "Overall (Crude)", OR = main_or, Low = main_ci[1], High = main_ci[2],
             Color = "#E41A1C", stringsAsFactors = FALSE),
  data.frame(Label = "Overall (Adjusted*)", OR = main_adj_or, Low = main_adj_ci[1], High = main_adj_ci[2],
             Color = "#377EB8", stringsAsFactors = FALSE),
  data.frame(Label = sprintf("Age <40 (n=%d)", subgroup_or$N[1]),
             OR = subgroup_or$OR_per_100mg[1], Low = subgroup_or$CI_low[1], High = subgroup_or$CI_high[1],
             Color = ifelse(subgroup_or$P_value[1] < 0.05, "#E41A1C", "#BEBEBE"),
             stringsAsFactors = FALSE),
  data.frame(Label = sprintf("Age 40-49 (n=%d)", subgroup_or$N[2]),
             OR = subgroup_or$OR_per_100mg[2], Low = subgroup_or$CI_low[2], High = subgroup_or$CI_high[2],
             Color = ifelse(subgroup_or$P_value[2] < 0.05, "#E41A1C", "#BEBEBE"),
             stringsAsFactors = FALSE),
  data.frame(Label = sprintf("Age 50-59 (n=%d)", subgroup_or$N[3]),
             OR = subgroup_or$OR_per_100mg[3], Low = subgroup_or$CI_low[3], High = subgroup_or$CI_high[3],
             Color = ifelse(subgroup_or$P_value[3] < 0.05, "#E41A1C", "#BEBEBE"),
             stringsAsFactors = FALSE),
  data.frame(Label = sprintf("Age 60-69 (n=%d)", subgroup_or$N[4]),
             OR = subgroup_or$OR_per_100mg[4], Low = subgroup_or$CI_low[4], High = subgroup_or$CI_high[4],
             Color = ifelse(subgroup_or$P_value[4] < 0.05, "#E41A1C", "#BEBEBE"),
             stringsAsFactors = FALSE),
  data.frame(Label = sprintf("Age 70+ (n=%d)", subgroup_or$N[5]),
             OR = subgroup_or$OR_per_100mg[5], Low = subgroup_or$CI_low[5], High = subgroup_or$CI_high[5],
             Color = ifelse(subgroup_or$P_value[5] < 0.05, "#E41A1C", "#BEBEBE"),
             stringsAsFactors = FALSE)
)

n_items <- nrow(forest_data)
x_lim <- c(0.5, max(c(forest_data$High, 3), na.rm = TRUE))

plot(1, type = "n", xlim = x_lim, ylim = c(0.5, n_items + 0.5),
     xlab = "OR per 100mg increase (95% CI)", ylab = "", yaxt = "n", log = "x",
     main = "氯氮平剂量-吸入性肺炎风险: 亚组分析")

abline(v = 1, lty = 2, col = "grey50", lwd = 2)

for (i in 1:n_items) {
  y_pos <- n_items - i + 1
  points(forest_data$OR[i], y_pos, pch = 16, col = forest_data$Color[i], cex = 1.5)
  segments(forest_data$Low[i], y_pos, forest_data$High[i], y_pos,
           col = forest_data$Color[i], lwd = 3)
  axis(2, at = y_pos, labels = forest_data$Label[i], las = 2, cex.axis = 0.8,
       col.axis = forest_data$Color[i])
  text(forest_data$High[i] * 1.05, y_pos,
       labels = sprintf("%.2f (%.2f-%.2f)", forest_data$OR[i], forest_data$Low[i], forest_data$High[i]),
       cex = 0.65, adj = 0)
}

legend("bottomright", legend = c("Significant (p<0.05)", "Not significant"),
       col = c("#E41A1C", "#BEBEBE"), pch = 16, cex = 0.8, bty = "n")

dev.off()
message("  -> *_subgroup_forest.png")

# --- 图3: 剂量分组堆叠条形图 + 风险热图 ---
png(file.path(output_dir, paste0(pfx, "_dose_group_analysis.png")),
    width = 1400, height = 600, res = 130)

par(mfrow = c(1, 3), mar = c(7, 5, 4, 2))

# 3A: 各剂量组吸入性肺炎比例
bp <- barplot(abs_risk$Asp_Risk, names.arg = abs_risk$dose_group,
              col = ifelse(abs_risk$dose_group %in% c("500-599","600-899","900+"),
                           "#E41A1C", "#4575B4"),
              border = "white",
              main = "各剂量组吸入性肺炎风险 (%)",
              xlab = "日剂量 (mg)", ylab = "吸入性肺炎比例 (%)", las = 2)
text(bp, abs_risk$Asp_Risk + 1.2,
     labels = sprintf("%.1f%%", abs_risk$Asp_Risk), cex = 0.75)
text(bp, abs_risk$Asp_Risk - 2,
     labels = sprintf("n=%d", abs_risk$Total), cex = 0.6, col = "white")

# 3B: 各剂量组致死率趋势
bp2 <- barplot(abs_risk$Fatal_Risk, names.arg = abs_risk$dose_group,
               col = "#D73027", border = "white",
               main = "各剂量组致死率 (%)",
               xlab = "日剂量 (mg)", ylab = "致死率 (%)", las = 2)
text(bp2, abs_risk$Fatal_Risk + 1.5,
     labels = sprintf("%.1f%%", abs_risk$Fatal_Risk), cex = 0.75)

# 3C: 剂量阈值 OR (吸入性肺炎)
plot(threshold_results$Threshold, threshold_results$OR,
     type = "b", pch = 16, col = ifelse(threshold_results$P_value < 0.05, "#E41A1C", "#BEBEBE"),
     lwd = 2, main = "剂量阈值 OR (吸入性肺炎)",
     xlab = "剂量阈值 (mg/d)", ylab = "OR (Above vs Below)")
segments(threshold_results$Threshold,
         threshold_results$CI_low,
         threshold_results$Threshold,
         threshold_results$CI_high,
         col = ifelse(threshold_results$P_value < 0.05, "#E41A1C", "#BEBEBE"),
         lwd = 1.5)
abline(h = 1, lty = 2, col = "grey50")
# 标注显著点
sig_th <- threshold_results[threshold_results$P_value < 0.05, ]
if (nrow(sig_th) > 0) {
  text(sig_th$Threshold, sig_th$OR + sig_th$OR * 0.08,
       labels = sprintf("%.1f*", sig_th$OR), cex = 0.7, col = "#E41A1C")
} else {
  # 标最接近显著性的
  best <- threshold_results[which.min(threshold_results$P_value), ]
  text(best$Threshold, best$OR + best$OR * 0.08,
       labels = sprintf("%.1f (p=%.2f)", best$OR, best$P_value), cex = 0.7, col = "darkorange")
}

dev.off()
message("  -> *_dose_group_analysis.png")

# --- 图4: 剂量分布 + 给药方案 ---
png(file.path(output_dir, paste0(pfx, "_distribution.png")),
    width = 1200, height = 500, res = 130)

par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

# 4A: 剂量分布直方图 + 密度(分组: 吸入性肺炎/非)
hist(dose_full$daily_dose[dose_full$is_aspiration == 0],
     breaks = 30, col = rgb(0.27, 0.46, 0.71, 0.5), border = "white",
     main = "剂量分布: 吸入性肺炎 vs 非吸入性肺炎",
     xlab = "日剂量 (mg)", ylab = "频数", xlim = c(0, 1000))
hist(dose_full$daily_dose[dose_full$is_aspiration == 1],
     breaks = 20, col = rgb(0.89, 0.10, 0.11, 0.5), border = "white", add = TRUE)
legend("topright",
       legend = c("非吸入性肺炎", "吸入性肺炎"),
       fill = c(rgb(0.27, 0.46, 0.71, 0.5), rgb(0.89, 0.10, 0.11, 0.5)),
       cex = 0.8, bty = "n")
# 添加密度线
dens_asp <- density(dose_full$daily_dose[dose_full$is_aspiration == 1], n = 512)
dens_non <- density(dose_full$daily_dose[dose_full$is_aspiration == 0], n = 512)
# 缩放密度到频数
scale_factor <- length(dose_full$daily_dose[dose_full$is_aspiration == 0]) * 5
lines(dens_non$x, dens_non$y * scale_factor / max(dens_non$y) * 0.3, col = "#4575B4", lwd = 2)
if (sum(dose_full$is_aspiration == 1) > 10) {
  lines(dens_asp$x, dens_asp$y * scale_factor / max(dens_asp$y) * 0.3, col = "#E41A1C", lwd = 2)
}

# 4B: 给药方案饼图
par(mar = c(1, 1, 4, 1))
schedule_counts <- table(dose_full$dosing_schedule)
pie(schedule_counts,
    labels = sprintf("%s\nn=%d (%.1f%%)",
                     names(schedule_counts), schedule_counts,
                     100 * schedule_counts / sum(schedule_counts)),
    col = c("#4575B4", "#E41A1C", "#BEBEBE"),
    main = "给药方案分布",
    cex = 0.8)

dev.off()
message("  -> *_distribution.png")

# ============================================================
# Part 9: 导出结果
# ============================================================
message("\n=========================================")
message("Part 9: 导出分析结果")
message("=========================================")

# CSV导出
fwrite(dose_full[, .(primaryid, daily_dose, dose_group, dose_freq, dosing_schedule,
                      age_num, age_group, sex, reporter_country, i_f_code, fatal,
                      is_aspiration, is_pneumonia, is_LRTI,
                      pt_primary, dechal, rechal)],
       file.path(output_dir, "dose_response_analysis_data.csv"), bom = TRUE)

fwrite(abs_risk, file.path(output_dir, "dose_response_absolute_risk.csv"), bom = TRUE)
fwrite(threshold_results, file.path(output_dir, "dose_threshold_analysis.csv"), bom = TRUE)
fwrite(subgroup_or, file.path(output_dir, "dose_subgroup_or.csv"), bom = TRUE)
fwrite(rcs_asp_result, file.path(output_dir, "rcs_aspiration_curve_data.csv"), bom = TRUE)
fwrite(rcs_fatal_result, file.path(output_dir, "rcs_fatal_curve_data.csv"), bom = TRUE)
fwrite(rcs_pneu_result, file.path(output_dir, "rcs_pneumonia_curve_data.csv"), bom = TRUE)

# Excel 报告
if (requireNamespace("openxlsx", quietly = TRUE)) {
  library(openxlsx)
  output_xlsx <- file.path(data_dir, "clozapine_dose_response_report.xlsx")

  wb <- createWorkbook()

  # Sheet 1: 总览
  addWorksheet(wb, "分析总览")
  overview <- data.frame(
    指标 = c("有效剂量记录数", "中位日剂量(mg)", "IQR(mg)",
            "吸入性肺炎OR_per_100mg", "吸入性肺炎OR_95%CI",
            "吸入性肺炎Adjusted_OR", "吸入性肺炎Adj_95%CI",
            "致死OR_per_100mg", "致死OR_95%CI",
            "非线性检验_p(吸入)", "非线性检验_p(致死)", "非线性检验_p(肺炎)",
            "剂量×年龄交互p", "剂量×性别交互p"),
    数值 = c(
      nrow(dose_full),
      sprintf("%.0f", median(dose_full$daily_dose)),
      sprintf("%.0f-%.0f", quantile(dose_full$daily_dose, 0.25), quantile(dose_full$daily_dose, 0.75)),
      sprintf("%.3f", main_or),
      sprintf("%.3f-%.3f", main_ci[1], main_ci[2]),
      sprintf("%.3f", main_adj_or),
      sprintf("%.3f-%.3f", main_adj_ci[1], main_adj_ci[2]),
      sprintf("%.3f", or_fatal_100),
      sprintf("%.3f-%.3f", exp(confint(m3)["daily_dose", 1] * 100), exp(confint(m3)["daily_dose", 2] * 100)),
      format.pval(lrt_nonlinear$`Pr(>Chi)`[2], digits = 4),
      format.pval(NA, digits = 4),
      format.pval(lrt_nonlinear3$`Pr(>Chi)`[2], digits = 4),
      format.pval(coef(summary(m_interact))["daily_dose:age_num", "Pr(>|z|)"], digits = 4),
      format.pval(coef(summary(m_sex_interact))["dose_c100:sex_bin", "Pr(>|z|)"], digits = 4)
    )
  )
  writeData(wb, "分析总览", overview)

  # Sheet 2: 剂量分组风险
  addWorksheet(wb, "剂量分组绝对风险")
  writeData(wb, "剂量分组绝对风险", abs_risk)

  # Sheet 3: 剂量阈值分析
  addWorksheet(wb, "剂量阈值分析")
  writeData(wb, "剂量阈值分析", threshold_results)

  # Sheet 4: 亚组分析
  addWorksheet(wb, "亚组分析")
  writeData(wb, "亚组分析", subgroup_or)

  # Sheet 5: RCS预测曲线
  addWorksheet(wb, "RCS吸入性肺炎")
  writeData(wb, "RCS吸入性肺炎", rcs_asp_result)

  addWorksheet(wb, "RCS致死")
  writeData(wb, "RCS致死", rcs_fatal_result)

  addWorksheet(wb, "RCS肺炎")
  writeData(wb, "RCS肺炎", rcs_pneu_result)

  # Sheet 6: 给药方案
  addWorksheet(wb, "给药方案分析")
  writeData(wb, "给药方案分析", sched_dose)

  # Sheet 7: 原始剂量数据
  addWorksheet(wb, "分析数据集")
  writeData(wb, "分析数据集", dose_full[, .(primaryid, daily_dose, dose_group, dose_freq,
                                        dosing_schedule, age_num, age_group, sex,
                                        reporter_country, fatal, is_aspiration,
                                        is_pneumonia, is_LRTI, pt_primary)])

  saveWorkbook(wb, output_xlsx, overwrite = TRUE)
  message(sprintf("Excel报告: %s", output_xlsx))
}

# ============================================================
# 最终汇总
# ============================================================
cat("\n============================================================\n")
cat("       氯氮平肺部感染 — 剂量-反应关系深度分析完成\n")
cat("============================================================\n\n")

cat("--- 核心发现 ---\n\n")

cat("1. 连续剂量效应:\n")
cat(sprintf("   每增加100mg/d, 吸入性肺炎风险:\n"))
cat(sprintf("     Crude OR = %.3f (95%%CI: %.3f-%.3f)\n", main_or, main_ci[1], main_ci[2]))
cat(sprintf("     Adjusted OR = %.3f (95%%CI: %.3f-%.3f)\n", main_adj_or, main_adj_ci[1], main_adj_ci[2]))
cat(sprintf("   每增加100mg/d, 致死风险:\n"))
m3_ci <- confint.default(m3)["daily_dose", ]
cat(sprintf("     Adjusted OR = %.3f (95%%CI: %.3f-%.3f)\n",
            or_fatal_100,
            exp(m3_ci[1] * 100),
            exp(m3_ci[2] * 100)))

cat("\n2. 非线性检验:\n")
cat(sprintf("   吸入性肺炎剂量-反应非线性: p=%s\n",
            format.pval(lrt_nonlinear$`Pr(>Chi)`[2], digits = 4)))
cat(sprintf("   肺炎(总)剂量-反应非线性: p=%s\n",
            format.pval(lrt_nonlinear3$`Pr(>Chi)`[2], digits = 4)))

cat("\n3. 剂量阈值:\n")
cat(sprintf("   吸入性肺炎最佳阈值: %.0f mg/d (OR=%.2f, 95%%CI: %.2f-%.2f)\n",
            best_th$Threshold, best_th$OR, best_th$CI_low, best_th$CI_high))

cat("\n4. 亚组异质性:\n")
cat(sprintf("   剂量×年龄交互: p=%s\n",
            format.pval(coef(summary(m_interact))["daily_dose:age_num", "Pr(>|z|)"], digits = 4)))
cat(sprintf("   剂量×性别交互: p=%s\n",
            format.pval(coef(summary(m_sex_interact))["dose_c100:sex_bin", "Pr(>|z|)"], digits = 4)))

cat("\n5. 高剂量组 (>=600mg/d):\n")
high_dose_risk <- abs_risk[dose_group %in% c("600-899","900+")]
cat(sprintf("   吸入性肺炎风险: %.1f%%-%.1f%%\n",
            min(high_dose_risk$Asp_Risk), max(high_dose_risk$Asp_Risk)))
cat(sprintf("   致死率: %.1f%%-%.1f%%\n",
            min(high_dose_risk$Fatal_Risk), max(high_dose_risk$Fatal_Risk)))

cat("\n导出文件:\n")
cat("  clozapine_dose_response_report.xlsx    — 完整Excel报告(7 sheets)\n")
cat("  dose_response_analysis_data.csv        — 分析数据集\n")
cat("  dose_response_absolute_risk.csv        — 剂量组绝对风险\n")
cat("  dose_threshold_analysis.csv            — 阈值分析\n")
cat("  dose_subgroup_or.csv                   — 亚组OR\n")
cat("  rcs_*_curve_data.csv                   — RCS预测曲线数据\n")
cat("  clozapine_dose_response_*.png          — 可视化\n")
cat("============================================================\n")
