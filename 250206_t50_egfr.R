# Association of T50 and eGFR at AKS cohort
# Creation: 2025-02-06
# Author: Shun Ishibashi, Claude 3.5 Sonnet
# Purpose: Analyse the association between T50 and eGFR

# sessionInfo()
# R version 4.3.2 (2023-10-31)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS 15.2

################################################################################
# Install and Load----
################################################################################
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  rio,        # データを読み込むためのパッケージ
  here,       # 相対ファイルパスを設定するためのパッケージ
  tidyverse,  # データ管理と可視化のためのパッケージ
  janitor,    # データ前処理と表の作成のためのパッケージ
  skimr,    　# データの概要を把握
  lme4,       # 混合効果モデル
  lmerTest,   # 混合効果モデルにおけるp value
  performance,# check_model
  naniar,     # 欠測の可視化
  tableone,   # Table1作成
  gtsummary   # 表作成
)

################################################################################
# Data import and cleaning----
################################################################################
# data setをimport
linelist_raw <- import(
  here("Data", "0. AKS_dataset.xlsx"), sheet = "Dataset")

# R用にdata cleaning
linelist <- linelist_raw %>% 
  # 詳細情報やmemo情報は重なっているので今回は削除
  select(!starts_with(c("memo_", "Memo_", "詳細_"))) %>% 
  clean_names() %>% # 列名標準化
  # KDIGO分類
  mutate(kdigo_cr_18 = as.factor(case_when(
    is.na(date_18)          ~ NA_character_,
    e_gf_rcr_18 >= 90       ~ "1",
    e_gf_rcr_18 >= 60       ~ "2",
    e_gf_rcr_18 >= 45       ~ "3a",
    e_gf_rcr_18 >= 30       ~ "3b",
    e_gf_rcr_18 >= 15       ~ "4",
    TRUE                    ~ "5"
  ))) %>%
  mutate(kdigo_cys_18 = as.factor(case_when(
    is.na(date_18)          ~ NA_character_,
    e_gf_rcys_c_18 >= 90    ~ "1",
    e_gf_rcys_c_18 >= 60    ~ "2",
    e_gf_rcys_c_18 >= 45    ~ "3a",
    e_gf_rcys_c_18 >= 30    ~ "3b",
    e_gf_rcys_c_18 >= 15    ~ "4",
    TRUE                    ~ "5"
  ))) %>%
  mutate(kdigo_crcys_18 = as.factor(case_when(
    is.na(date_18)          ~ NA_character_,
    e_gf_rcr_cys_c_18 >= 90 ~ "1",
    e_gf_rcr_cys_c_18 >= 60 ~ "2",
    e_gf_rcr_cys_c_18 >= 45 ~ "3a",
    e_gf_rcr_cys_c_18 >= 30 ~ "3b",
    e_gf_rcr_cys_c_18 >= 15 ~ "4",
    TRUE                    ~ "5"
  ))) %>%
  mutate(kdigo_uacr_18 = as.factor(case_when(
    is.na(date_18) ~ NA_character_,
    uacr_18 < 30   ~ "1",
    uacr_18 < 300  ~ "2",
    TRUE           ~ "3"
  ))) %>% 
  mutate(ckd_cr_18 = factor(case_when(
    is.na(date_18) ~ NA_character_,
    (kdigo_cr_18 == "1" | kdigo_cr_18 == "2") & kdigo_uacr_18 == "1" ~ "0",
    TRUE                                                             ~ "1"
  ), levels = c("0", "1"), labels = c("non-CKD", "CKD"))) %>% 
  mutate(ckd_cys_18 = factor(case_when(
    is.na(date_18) ~ NA_character_,
    (kdigo_cys_18 == "1" | kdigo_cys_18 == "2") & kdigo_uacr_18 == "1" ~ "0",
    TRUE                                                               ~ "1"
  ), levels = c("0", "1"), labels = c("non-CKD", "CKD"))) %>% 
  mutate(ckd_crcys_18 = factor(case_when(
    is.na(date_18) ~ NA_character_,
    (kdigo_crcys_18 == "1" | kdigo_crcys_18 == "2") & kdigo_uacr_18 == "1" ~ "0",
    TRUE                                                                   ~ "1"
  ), levels = c("0", "1"), labels = c("non-CKD", "CKD"))) %>%
  # 喫煙を0 Never、1 Former, 2 Currentに変更
  mutate(across(
    .cols = c(smoking_18, smoking_19, smoking_20, smoking_21, smoking_22, smoking_23),
    .names = "{.col}",
    ~ factor(case_when(
      . == 0 ~ "Never",  # never smoker
      . == 1 & get(paste0("current_", cur_column())) == 1 ~ "Current",  # current smoker
      . == 1 & get(paste0("current_", cur_column())) == 0 ~ "Former",  # ex-smoker
      TRUE ~ NA_character_  # その他の場合（矛盾するデータなど）
    ), levels = c("Never", "Former", "Current")))) %>% 
  mutate(across(
    .cols = starts_with("current_smoking_"),
    ~ factor(.)
  )) %>% 
  # select(-starts_with("current_")) %>% 
  # DM(Diabetes)有無を判定
  mutate(dm_18 = factor(case_when(
    is.na(hb_a1c_ngsp_18) & is.na(dm_med_18) ~ NA_character_,
    hb_a1c_ngsp_18 >= 6.5 | dm_med_18 == 1   ~ "Yes",
    TRUE                                     ~ "No"
  ))) %>%
  # HTN(Hypertension)有無を判定
  mutate(htn_18 = factor(case_when(
    is.na(b_sbp_18) & is.na(anti_ht_med_18) ~ NA_character_,
    b_sbp_18 >= 130 | anti_ht_med_18 == 1   ~ "Yes",
    TRUE                                    ~ "No"
  ))) %>%
  # DL(Dyslipidemia)有無を判定
  mutate(dl_18 = factor(case_when(
    is.na(ldl_c_18) & is.na(ld_med_18)  ~ NA_character_,
    hdl_c_18 < 40 | ldl_c_18 >= 140 |
      tg_18 >= 150 | ld_med_18 == 1     ~ "Yes",
    TRUE                                ~ "No"
  ))) %>% 
  # T50 fitは2018年と2019年データがあるため、averageも作成
  mutate(t50_fit_ave_18 = case_when(
    is.na(t50_fit_18) & !is.na(t50_fit_19) ~ t50_fit_19, # 18NA、19あり ~ 19
    is.na(t50_fit_19) ~ t50_fit_18, # 19NA ~ 18
    # 18あり、19あり ~ 18と19の平均値
    !is.na(t50_fit_18) & !is.na(t50_fit_19) ~ round((t50_fit_18 + t50_fit_19) / 2, digits = 2),
    TRUE ~ NA_integer_
  )) %>% 
  # T50の変化量
  mutate(t50_change_absolute_18 = case_when(
    is.na(t50_fit_18) | is.na(t50_fit_19) ~ NA_real_,
    TRUE                                  ~ t50_fit_19 - t50_fit_18
  )) %>% 
  mutate(t50_change_percentage_18 = case_when(
    is.na(t50_fit_18) | is.na(t50_fit_19) ~ NA_real_,
    TRUE                                  ~ (t50_change_absolute_18 / t50_fit_18) * 100
  )) %>% 
  # T50のカテゴリ作成
  # 2分位
  mutate(t50_2tile_18 = case_when(
    t50_fit_18 <= median(t50_fit_18, na.rm = TRUE) ~ "Low",
    TRUE ~ "High"
  ) %>% 
    factor(levels = c("High", "Low"))
  ) %>% 
  # 3分位
  mutate(t50_3tile_18 = case_when(
    t50_fit_18 <= quantile(t50_fit_18, 1/3, na.rm = TRUE) ~ "Low",
    t50_fit_18 <= quantile(t50_fit_18, 2/3, na.rm = TRUE) ~ "Middle",
    TRUE ~ "High"
  ) %>% 
    factor(levels = c("High", "Middle", "Low"))
  ) %>% 
  mutate(t50_4tile_18 = case_when(
    t50_fit_18 <= quantile(t50_fit_18, 0.25, na.rm = TRUE) ~ "Q1",
    t50_fit_18 <= quantile(t50_fit_18, 0.50, na.rm = TRUE) ~ "Q2",
    t50_fit_18 <= quantile(t50_fit_18, 0.75, na.rm = TRUE) ~ "Q3",
    TRUE ~ "Q4"
  ) %>% 
    factor(levels = c("Q4", "Q3", "Q2", "Q1"))
  ) %>% 
  filter(!is.na(date_18)) %>% # 初回検査なし(date_18がNA)の参加者を除外
  # AKSP-030Cは_23に腎移植ありのため、同患者の2023年の測定データは全てNAに
  mutate(across(ends_with("_23"), ~ if_else(grepl("^AKSP-030C", id), NA, .))) %>%
  mutate(row_id = row_number()) %>% # 通し番号を作成
  select(row_id, everything()) %>% # 通し番号を先頭
  { print(names(.)[1:1000]); invisible(.) } # パイプ後に表示可能

linelist %>% 
  select(contains("t50")) %>% 
  glimpse()
################################################################################
# T50用の変数選択----
################################################################################
# 参考文献：https://doi.org/10.1093/ckj/sfae343
# T50と正の相関が予想
## Alb, Mg, HCO3-, TG, HDL-C, Cre, 1,25(OH)2D, 25(OH)D, fetuin-A, pyrophosphate
# T50と負の相関が予想
## P, iCa, PTH, Urea, CRP
# CKDと関連が予想
## UACR

# 共変量のリスト作成
# IDリスト
id_vars <- c("row_id", "id")
# ベースライン連続変数
baseline_continuous_vars <- c("age_18")
# ベースラインカテゴリ変数
baseline_categorical_vars <- c(
  "sex_18", "smoking_18", "current_smoking_18", "dm_18", "htn_18", "dl_18",
  "kdigo_crcys_18", "kdigo_cr_18", "kdigo_cys_18",
  "kdigo_uacr_18",
  "ckd_crcys_18", "ckd_cr_18", "ckd_cys_18",
  "t50_2tile_18", "t50_3tile_18", "t50_4tile_18"
)
# 時間依存性(全て連続変数)
time_dependent_vars <- c(
  "egfr_cr", "egfr_cys", "egfr_crcys",
  "u_b2_mg_cr", "u_nag_cr", "ulfcr",
  "s_alb", "s_mg", "tg", "hdl_c", "ldl_c", "t_c", "x1_25_oh_2d",
  "s_pi", "s_ca", "pth_int", "bun", "hs_crp", "uacr",
  "b_sbp", "b_dbp", "bmi"
)
# 時間依存性連続変数のベースライン
time_dependent_vars_18 <- paste0(time_dependent_vars, "_18")
# ベースライン連続変数(P関連)
baseline_p_vars <- c(
  "t50_fit_18", "t50_fit_ave_18",
  "s_fgf23_18", "fe_pi_18", "fe_pi_fgf23_18", 
  "p_cpp_18", "s_cpp_18", "trp_18", "tmp_gfr_18",
  "t50_change_absolute_18", "t50_change_percentage_18"
)
# ベースライン連続変数全て
baseline_vars <- c(time_dependent_vars_18, baseline_p_vars)
# slope変数
slope_vars <- c("n", "duration",
                "slope_cr", "slope_cys", "slope_crcys")

# 新たな列名
list_cont <- baseline_continuous_vars %>% 
  c(baseline_p_vars) %>% 
  {sub("_18$", "", .)} %>%
  c(time_dependent_vars) %>% 
  {paste0("base_", .)}
list_cat <- baseline_categorical_vars %>%
  {sub("_18$", "", .)} %>%
  {paste0("base_cat_", .)}

# ラベルの定義
var_labels <- c(
  # ベースライン連続変数
  base_age = "Age (y)",
  # ベースラインP関連変数
  base_t50_fit = "T50 (min)",
  base_t50_fit_ave = "T50 2y Average (min)",
  base_s_fgf23 = "sFGF23 (pg/mL)",
  base_fe_pi = "FEPi (%)",
  base_fe_pi_fgf23 = "FEPi/FGF23 (mL/pg%)",
  base_p_cpp = "Plasma CPP (a.u.)",
  base_s_cpp = "Serum CPP (a.u.)",
  base_trp = "TRP",
  base_tmp_gfr = "TmP/GFR (mg/dL)",
  base_s_alb = "Albumin (g/dL)",
  base_s_mg = "Magnesium (mg/dL)",
  base_tg = "Triglycerides (mg/dL)",
  base_hdl_c = "HDL Cholesterol (mg/dL)",
  base_ldl_c = "LDL Cholesterol (mg/dL)",
  base_t_c = "Total Cholesterol (mg/dL)",
  base_x1_25_oh_2d = "1,25(OH)2D (pg/mL)",
  base_s_pi = "Serum Phosphate (mg/dL)",
  base_s_ca = "Serum Calcium (mg_dL)",
  base_pth_int = "intact PTH (pg/mL)",
  base_bun = "BUN (mg/dL)",
  base_hs_crp = "hs-CRP (mg/dL)",
  base_egfr_cr = "eGFRcr (mL/min/1.73m²)",
  base_egfr_cys = "eGFRcys (mL/min/1.73m²)",
  base_egfr_crcys = "eGFRcrcys (mL/min/1.73m²)",
  base_uacr = "UACR (mg/gCr)",
  base_b_sbp = "Systolic Blood Pressure (mmHg)",
  base_b_dbp = "Diastolic Bood Pressure (mmHg)",
  base_bmi = "Body Mass Index (kg/m²)",
  base_u_b2_mg_cr = "Urinary β2-microglobulin (μg/gCr)",
  base_u_nag_cr = "Urinary NAG (U/gCr)",
  base_ulfcr = "Urinary L-FABP (μg/gCr)",
  base_t50_change_absolute = "Absolute Change of T50 (min/year)",
  base_t50_change_percentage = "Change of T50 (%)",
  # 時間依存性変数
  s_alb = "Albumin (g/dL)",
  s_mg = "Magnesium (mg/dL)",
  tg = "Triglycerides (mg/dL)",
  hdl_c = "HDL Cholesterol (mg/dL)",
  ldl_c = "LDL Cholesterol (mg/dL)",
  t_c = "Total Cholesterol (mg/dL)",
  x1_25_oh_2d = "1,25(OH)2D (pg/mL)",
  s_pi = "Serum Phosphate (mg/dL)",
  s_ca = "Serum Calcium (mg_dL)",
  pth_int = "intact PTH (pg/mL)",
  bun = "BUN (mg/dL)",
  hs_crp = "hs-CRP (mg/dL)",
  egfr_cr = "eGFRcr (mL/min/1.73m²)",
  egfr_cys = "eGFRcys (mL/min/1.73m²)",
  egfr_crcys = "eGFRcrcys (mL/min/1.73m²)",
  uacr = "UACR (mg/gCr)",
  b_sbp = "Systolic Blood Pressure (mmHg)", 
  b_dbp = "Diastolic Bood Pressure (mmHg)",
  bmi = "Body Mass Index (kg/m²)",
  u_b2_mg_cr = "Urinary β2-microglobulin (μg/gCr)",
  u_nag_cr = "Urinary NAG (U/gCr)",
  ulfcr = "Urinary L-FABP (μg/gCr)",
  # カテゴリカル変数
  base_cat_sex = "Sex",
  base_cat_smoking = "Smoking",
  base_cat_current_smoking = "Smoking(current)",
  base_cat_dm = "Diabetes",
  base_cat_htn = "Hypertension",
  base_cat_dl = "Dyslipidemia",
  base_cat_kdigo_crcys = "GFR Category(eGFRcrcys)",
  base_cat_kdigo_cr = "GFR Category(eGFRcr)",
  base_cat_kdigo_cys = "GFR Category(eGFRcys)",
  base_cat_kdigo_uacr = "UACR Category",
  base_cat_ckd_crcys = "CKD Status(eGFRcrcys)",
  base_cat_ckd_cr = "CKD Status(eGFRcr)",
  base_cat_ckd_cys = "CKD Status(eGFRcys)",
  base_cat_t50_2tile = "T50 Category by Median",
  base_cat_t50_3tile = "T50 Category by Tertiles",
  base_cat_t50_4tile = "T50 Category by Quartiles",
  # slope変数
  n = "Number of measurements",
  duration = "Measurement period (y)",
  slope_cr = "eGFRcr slope (mL/min/1.73m²/year)",
  slope_cys = "eGFRcys slope (mL/min/1.73m²/year)",
  slope_crcys = "eGFRcrcys slope (mL/min/1.73m²/year)"
)
# 
var_labels_tilde <- gsub(" = ", " ~ ", var_labels)
# 上記変数を抽出したデータセット作成
t50_df <- linelist %>% 
  # 必要な共変量の抽出
  select(
    all_of(id_vars),  # ID情報
    all_of(baseline_continuous_vars),  # Base line 共変量(連続変数)
    all_of(baseline_categorical_vars),  # Base line 共変量(カテゴリー変数)
    starts_with("date_"),
    contains("e_gf_r"),
    starts_with(paste0(time_dependent_vars, "_")),  # Time dependent 共変量
    all_of(baseline_p_vars)  # baseline 共変量(P関連)
  ) %>% 
  select(!starts_with("bun_s_cr_")) %>% 
  # Sexをfactor型へ
  mutate(sex_18 = factor(sex_18)) %>% 
  # eGFRの列名変更
  rename_with(~ gsub("e_gf_r", "egfr_", .)) %>% 
  rename_with(~ gsub("_cr_cys_c_", "_crcys_", .)) %>% 
  rename_with(~ gsub("_cys_c_", "_cys_", .)) %>% 
  # { glimpse(.); invisible(.) } # パイプ後に表示可能
  { print(names(.)[1:170]); invisible(.) } # パイプ後に表示可能

# id, dateをlong形式にして、ベースラインからの測定年数を計算
date_df <- t50_df %>% 
  select(all_of(id_vars), starts_with("date_")) %>% 
  pivot_longer(
    cols = -all_of(id_vars),
    names_pattern = "date_(.+)",
    names_to = "time_point",
    values_to = "date"
  ) %>% 
  mutate(time_point = as.numeric(time_point)) %>% 
  group_by(row_id) %>% 
  mutate(
    date_base = date[time_point == 18]) %>% 
  mutate(
    years = time_length(
      interval(date_base, date),
      "years"
    ) %>% round(2)
  ) %>%
  ungroup() %>% 
  select(all_of(id_vars), time_point, years) %>% 
  { glimpse(.); invisible(.) } # パイプ後に表示可能

# 時間依存性の変数をwide dataに
time_dependent_df <- t50_df %>%
  select(
    all_of(id_vars),
    starts_with(c(
      paste0(time_dependent_vars, "_")
    )),
    contains("egfr_")
  ) %>%
  pivot_longer(
    cols = -all_of(id_vars),
    names_pattern = "(.+)_([0-9]{2})",
    names_to = c("variable", "time_point"),
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = variable,
    values_from = value
  ) %>% 
  mutate(time_point = as.numeric(time_point)) %>% 
  { glimpse(.); invisible(.) } # パイプ後に表示可能

# ベースライン連続変数をwide dataに
baseline_continuous_df <- t50_df %>%
  select(
    all_of(id_vars),
    all_of(baseline_continuous_vars),
    all_of(baseline_vars)
  ) %>%
  pivot_longer(
    cols = -all_of(id_vars),
    names_pattern = "(.+)_([0-9]{2})",
    names_to = c("variable", "time_point"),
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = variable,
    values_from = value
  ) %>% 
  rename_with(~paste0("base_", .), -c(all_of(id_vars), "time_point")) %>% 
  mutate(time_point = as.numeric(time_point)) %>% 
  { glimpse(.); invisible(.) } # パイプ後に表示可能

# ベースラインカテゴリ変数をwide dataに
baseline_categorical_df <- t50_df %>%
  select(
    all_of(id_vars),
    all_of(baseline_categorical_vars)
  ) %>%
  pivot_longer(
    cols = -all_of(id_vars),
    names_pattern = "(.+)_([0-9]{2})",
    names_to = c("variable", "time_point"),
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = variable,
    values_from = value
  ) %>% 
  rename_with(~paste0("base_cat_", .), -c(all_of(id_vars), "time_point")) %>% 
  mutate(time_point = as.numeric(time_point)) %>% 
  { glimpse(.); invisible(.) } # パイプ後に表示可能

# それぞれの変数のデータリストを結合
long_df <- time_dependent_df %>%
  left_join(baseline_continuous_df, by = c(id_vars, "time_point")) %>%
  left_join(baseline_categorical_df, by = c(id_vars, "time_point")) %>% 
  # linelist_dateを併合
  left_join(
    date_df,
    by = c(id_vars, "time_point")
  ) %>% 
  # eGFRcrがNAの行を削除
  filter(!is.na(egfr_cr)) %>% 
  # 以下の処理でID毎にbaseline共変量を他のtimepointにも補完
  group_by(id) %>%
  mutate(
    across(
      starts_with("base_"),
      ~ case_when(
        time_point == "18" ~ .x,
        TRUE ~ first(.x[time_point == "18"])
      )
    )
  ) %>%
  ungroup() %>% 
  # 列の並び替え
  select(row_id, id, time_point, years,
         starts_with("base_"),
         everything()) %>% 
  { glimpse(.); invisible(.) } # パイプ後に表示可能

# eGFRの測定回数と期間を計算(eGFRcrとeGFRcysの測定は同時なのでcrのみで算定)
measurement_info <- long_df %>%
  group_by(row_id) %>%
  summarise(
    # 測定回数
    n = sum(!is.na(egfr_cr)),
    # 測定期間（年数）
    duration = case_when(
      n > 1 ~ diff(range(years[!is.na(egfr_cr)])),
      TRUE ~ 0
    )) %>% 
  { glimpse(.); invisible(.) } # パイプ後に表示可能

# 全ての情報を取り入れたlong data set
linelist_long <- long_df %>% 
  left_join(measurement_info,
            by = "row_id") %>% 
  select(row_id, id, n, duration, time_point, years, everything()) %>% 
  # eGFRを2回以上測定した患者に限定
  group_by(row_id) %>% 
  filter(n >= 2) %>% 
  ungroup() %>% 
  { glimpse(.); invisible(.) } # パイプ後に表示可能

# yearの分布確認
linelist_long %>% 
  filter(!is.na(egfr_cr)) %>% 
  group_by(time_point) %>% 
  summarise(
    n = n(),
    mean = mean(years, na.rm = TRUE),
    sd = sd(years, na.rm = TRUE),
    min = min(years, na.rm = TRUE),
    q1 = quantile(years, 0.25, na.rm = TRUE),
    median = median(years, na.rm = TRUE),
    q3 = quantile(years, 0.75, na.rm = TRUE),
    max = max(years, na.rm = TRUE)
  ) %>% 
  filter(time_point >= 18)

################################################################################
# eGFR slopeの計算----
################################################################################
glimpse(linelist_long)
# LMM fit at egfr_cr
model_cr <- lmer(egfr_cr ~ 1 + years + (1 + years | row_id),
                 data = linelist_long
                 , control = lmerControl(optimizer = "bobyqa")
)
# summary(model_cr)
## model check
# check_model(
#   model_cr,
#   check = c("reqq", "qq", "linearity", "homogeneity", "outliers", "pp_check")
# )

# LMM fit at egfr_cys
model_cys <- lmer(egfr_cys ~ 1 + years + (1 + years | row_id),
                  data = linelist_long
                  , control = lmerControl(optimizer = "bobyqa")
)
# summary(model_cys)
## model check
# check_model(
#   model_cys,
#   check = c("reqq", "qq", "linearity", "homogeneity", "outliers", "pp_check")
# )

# LMM fit at egfr_crcys
model_crcys <- lmer(egfr_crcys ~ 1 + years + (1 + years | row_id),
                    data = linelist_long
                    , control = lmerControl(optimizer = "bobyqa")
)
# summary(model_crcys)
## model check
# check_model(
#   model_crcys,
#   check = c("reqq", "qq", "linearity", "homogeneity", "outliers", "pp_check")
# )

# 条件付き係数（conditional coefficients）の出力
# eGFRcr
individual_coef_cr <- coef(model_cr)$row_id %>%
  as.data.frame() %>%
  rownames_to_column("row_id") %>%
  rename(
    intercept_cr = `(Intercept)`,
    slope_cr = years
  ) %>%
  mutate(row_id = as.numeric(row_id))
# eGFRcys
individual_coef_cys <- coef(model_cys)$row_id %>%
  as.data.frame() %>%
  rownames_to_column("row_id") %>%
  rename(
    intercept_cys = `(Intercept)`,
    slope_cys = years
  ) %>%
  mutate(row_id = as.numeric(row_id))
# eGFRcrcys
individual_coef_crcys <- coef(model_crcys)$row_id %>%
  as.data.frame() %>%
  rownames_to_column("row_id") %>%
  rename(
    intercept_crcys = `(Intercept)`,
    slope_crcys = years
  ) %>%
  mutate(row_id = as.numeric(row_id))
# 3つのデータフレームを結合
individual_coef <- individual_coef_cr %>%
  left_join(individual_coef_cys, by = "row_id") %>%
  left_join(individual_coef_crcys, by = "row_id") %>% 
  select(row_id, slope_cr, slope_cys, slope_crcys)
individual_egfr <- individual_coef %>% 
  left_join(measurement_info,
            by = "row_id") %>% 
  select(row_id, n, duration, slope_cr, slope_cys, slope_crcys)
glimpse(individual_coef)
# eGFR slopeデータの統合
linelist_long_slope <- linelist_long %>% 
  left_join(individual_coef,
            by = "row_id")
linelist_long_slope %>% 
  glimpse()

# Factorのlevels確認
levels(linelist_long_slope$base_cat_sex)

# FactorのLeverls再定義
linelist_long_slope <- linelist_long_slope %>%
  mutate(
    # 性別
    base_cat_sex = factor(base_cat_sex, levels = c("M", "F")),
    
    # 喫煙
    base_cat_smoking = factor(base_cat_smoking, 
                              levels = c("Never", "Former", "Current")),
    
    # 現在の喫煙
    base_cat_current_smoking = factor(base_cat_current_smoking,
                                      levels = c("0", "1")),
    
    # 疾患
    base_cat_dm = factor(base_cat_dm, levels = c("No", "Yes")),
    base_cat_htn = factor(base_cat_htn, levels = c("No", "Yes")),
    base_cat_dl = factor(base_cat_dl, levels = c("No", "Yes")),
    
    # CKD関連
    base_cat_kdigo_crcys = factor(base_cat_kdigo_crcys, 
                                  levels = c("1", "2", "3a", "3b", "4")),
    base_cat_kdigo_cr = factor(base_cat_kdigo_cr, 
                               levels = c("1", "2", "3a", "3b", "4")),
    base_cat_kdigo_cys = factor(base_cat_kdigo_cys, 
                                levels = c("1", "2", "3a", "3b", "4")),
    base_cat_kdigo_uacr = factor(base_cat_kdigo_uacr,
                                 levels = c("1", "2", "3")),
    base_cat_ckd_crcys = factor(base_cat_ckd_crcys,
                                levels = c("non-CKD", "CKD")),
    base_cat_ckd_cr = factor(base_cat_ckd_cys,
                             levels = c("non-CKD", "CKD")),
    base_cat_ckd_cys = factor(base_cat_ckd_cys,
                              levels = c("non-CKD", "CKD")),
    # T50分位
    base_cat_t50_2tile = factor(base_cat_t50_2tile,
                                levels = c("High", "Low")),
    base_cat_t50_3tile = factor(base_cat_t50_3tile,
                                levels = c("High", "Middle", "Low")),
    base_cat_t50_4tile = factor(base_cat_t50_4tile,
                                levels = c("Q4", "Q3", "Q2", "Q1"))
  )
################################################################################
# Visualize----
################################################################################
# ベースラインデータの抽出と数値変数の可視化
baseline_num_plot <- linelist_long_slope %>%
  filter(time_point == 18) %>%
  select(base_cat_ckd_cr, any_of(list_cont), any_of(slope_vars)) %>%
  pivot_longer(cols = -base_cat_ckd_cr, 
               names_to = "variable", 
               values_to = "value") %>%
  ggplot(aes(x = base_cat_ckd_cr, y = value)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = base_cat_ckd_cr), 
              alpha = 0.5, width = 0.2) +
  facet_wrap(~ variable, scales = "free_y", ncol = 4,
             labeller = labeller(variable = var_labels)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.text = element_text(size = 8)) +
  labs(x = "CKD Status", y = "Value",
       title = "Distribution of Numeric Variables by CKD Status at Baseline",
       caption = "Box plots show median, IQR, and range. Points represent individual observations.")

# カテゴリカル変数の可視化
baseline_cat_plot <- linelist_long %>%
  filter(time_point == 18) %>%
  select(base_cat_ckd_cr, any_of(list_cat)) %>%
  select(-contains("t50"),
         -contains("_crcys"),
         -contains("_cys")) %>% 
  pivot_longer(cols = -base_cat_ckd_cr, 
               names_to = "variable", 
               values_to = "value") %>%
  ggplot(aes(x = base_cat_ckd_cr, fill = value)) +
  geom_bar(position = "fill") +
  facet_wrap(~ variable, scales = "free_y", ncol = 4,
             labeller = labeller(variable = var_labels)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10)) +
  labs(x = "CKD Status", y = "Proportion",
       title = "Distribution of Categorical Variables by CKD Status at Baseline",
       fill = "Category")

glimpse(linelist_long)
# プロットの表示
print(baseline_num_plot)
print(baseline_cat_plot)


################################################################################
# 分布確認と変換----
################################################################################
glimpse(linelist_long_slope)

# ベースラインデータのヒストグラム
baseline_hist_plot <- linelist_long_slope %>%
  filter(time_point == 18) %>%
  select(any_of(list_cont), any_of(slope_vars)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable", 
    values_to = "value"
  ) %>%
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 30,
                 fill = "steelblue",
                 color = "white",
                 alpha = 0.7) +
  geom_density(color = "red", linewidth = 0.8) +
  facet_wrap(
    ~ variable, 
    scales = "free", 
    ncol = 4
    # ,labeller = labeller(variable = var_labels) # ここを入れるとラベリングあり
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)
  ) +
  labs(
    x = "Value",
    y = "Density",
    title = "Distribution of Numeric Variables at Baseline",
    caption = "Blue histogram shows distribution, red line shows kernel density estimate."
  )
print(baseline_hist_plot)

# Variables that need log transformation
log_vars <- c(
  "base_hs_crp", "base_uacr", "base_p_cpp", "base_s_cpp", "base_pth_int", "base_tg",
  "base_s_fgf23", "base_fe_pi", "base_fe_pi_fgf23",
  "base_u_b2_mg_cr", "base_u_nag_cr", "base_ulfcr"
)
normal_vars <- setdiff(c(list_cont, slope_vars), log_vars)
glimpse(linelist_long_slope)
# Log ransformation and Standardization
linelist_transformed <- linelist_long_slope %>% 
  mutate(
    across(
      all_of(log_vars),
      log,
      .names = "log_{.col}"
    )
  ) %>% 
  mutate(
    across(
      c(all_of(normal_vars), starts_with("log_")),
      \(x) as.numeric(scale(x)),
      .names = "scaled_{.col}"
    )
  )

# 標準化前の要約統計量を保存
scaling_info <- linelist_long_slope %>%
  summarise(
    across(
      c(all_of(normal_vars), starts_with("log_")),
      list(
        mean = \(x) mean(x, na.rm = TRUE),
        sd = \(x) sd(x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  )
glimpse(linelist_transformed)
glimpse(scaling_info)

# log変換前後のヒストグラム比較
log_comparison_plot <- linelist_transformed %>%
  filter(time_point == 18) %>%
  select(
    all_of(log_vars),
    starts_with("log")
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    transform_type = if_else(
      str_detect(variable, "^log_"),
      "After",
      "Before"
    ),
    variable = str_remove(variable, "^log_"),
    # 表示用の結合変数を作成
    plot_group = paste(variable, transform_type)
  ) %>%
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 30,
                 fill = "steelblue",
                 color = "white",
                 alpha = 0.7) +
  geom_density(color = "red", linewidth = 0.8) +
  facet_wrap(
    ~ plot_group,
    scales = "free",
    ncol = 2
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    x = "Value",
    y = "Density",
    title = "Log Transformation Comparison",
    caption = "Blue histogram shows distribution, red line shows kernel density estimate."
  )
print(log_comparison_plot)

################################################################################
# Table1作成----
################################################################################
# visualizeプロットの表示
print(baseline_num_plot)
print(baseline_cat_plot)

glimpse(linelist_long_slope)
levels(linelist_long_slope$base_cat_t50_3tile)

# CKD有無で層別化
Table1_ckd <- CreateTableOne(vars = c(list_cont, list_cat, slope_vars),
                             factorVars = list_cat,
                             includeNA = TRUE,
                             addOverall = TRUE,
                             smd = TRUE,
                             data = subset(linelist_transformed,
                                           time_point == 18),
                             testExact = fisher.test,
                             testNonNormal = kruskal.test,
                             strata = "base_cat_ckd_cr")

print(Table1_ckd,
      showAllLevels = TRUE,
      varLabels = TRUE,
      # smd = TRUE,
      explain = FALSE,
      catDigits = 1,
      contDigits = 1,
      pDigits = 3,
      test = TRUE,
      format = "pf",
      nonnormal = log_vars,
      printToggle	= TRUE
) 

# TableOne to DataFrame to Tibble
Table1_ckd_df <- print(Table1_ckd,
                       showAllLevels = TRUE,
                       explain = FALSE,
                       catDigits = 1,
                       contDigits = 1,
                       pDigits = 3,
                       test = TRUE,
                       # missing = TRUE,
                       format = "pf",
                       nonnormal = log_vars,
                       printToggle = FALSE) %>%
  as.data.frame() %>%  # まずdata.frameに変換（行名を保持するため）
  rownames_to_column(var = "Variable") %>%  # 行名を列に変換
  as_tibble() %>%  # tibbleに変換
  mutate(Variable = case_when(
    Variable == "n" ~ "Number",
    Variable == "n.1" ~ "n",
    str_starts(Variable, "X") ~ "",  # Xで始まる変数を空文字に
    TRUE ~ Variable  # その他の変数はそのまま
  ))
  

# 変数名を置換
for (name in names(var_labels)) {
  Table1_ckd_df$Variable <- gsub(paste0("^", name, "$"), var_labels[name], Table1_ckd_df$Variable)
}

# 出力
write.csv(Table1_ckd_df, 
          file = here("Archive", 
                      paste0(format(Sys.Date(), "%y%m%d"), "_Table1_ckd.csv")))

################################################################################
# Table1 T50で層別化----
################################################################################
# T50で層別化
# 二分位
Table1_t50_2 <- CreateTableOne(vars = c(list_cont, list_cat, slope_vars),
                               factorVars = list_cat,
                               includeNA = TRUE,
                               addOverall = TRUE,
                               smd = TRUE,
                               data = subset(linelist_transformed,
                                             time_point == 18),
                               testExact = fisher.test,
                               testNonNormal = kruskal.test,
                               strata = "base_cat_t50_2tile")
print(Table1_t50_2,
      showAllLevels = TRUE,
      smd = TRUE,
      explain = TRUE,
      catDigits = 1,
      contDigits = 1,
      pDigits = 3,
      test = TRUE,
      format = "pf"
)
# 三分位
Table1_t50_3 <- CreateTableOne(vars = c(list_cont, list_cat, slope_vars),
                               factorVars = list_cat,
                               includeNA = TRUE,
                               addOverall = TRUE,
                               smd = TRUE,
                               data = subset(linelist_transformed,
                                             time_point == 18),
                               testExact = fisher.test,
                               testNonNormal = kruskal.test,
                               strata = "base_cat_t50_3tile")
print(Table1_t50_3,
      showAllLevels = TRUE,
      smd = TRUE,
      explain = TRUE,
      catDigits = 1,
      contDigits = 1,
      pDigits = 3,
      test = TRUE,
      format = "pf"
)

# 四分位
Table1_t50_4 <- CreateTableOne(vars = c(list_cont, list_cat, slope_vars),
                               factorVars = list_cat,
                               includeNA = TRUE,
                               addOverall = TRUE,
                               smd = TRUE,
                               data = subset(linelist_transformed,
                                             time_point == 18),
                               testExact = fisher.test,
                               testNonNormal = kruskal.test,
                               strata = "base_cat_t50_4tile")
print(Table1_t50_4,
      showAllLevels = TRUE,
      smd = TRUE,
      explain = TRUE,
      catDigits = 1,
      contDigits = 1,
      pDigits = 3,
      test = TRUE,
      format = "pf"
)


Table1_df <- as_tibble(print(Table1_ckd,
                             showAllLevels = TRUE,
                             smd = TRUE,
                             explain = FALSE,
                             catDigits = 1,
                             contDigits = 1,
                             pDigits = 3,
                             test = TRUE,
                             format = "pf"))

################################################################################
# T50と経時的なeGFRとの関連性----
################################################################################
glimpse(linelist_transformed)
glimpse(scaling_info)
# Level 1 : eGFR_ij = β_0j + β_1j×Time_ij + ε_ij
# Level 2 : β_0j = γ_00 + γ_01×T50_j + μ_0j
#           β_1j = γ_10 + γ_11×T50_j + μ_1j
# γ_00 = intercepts, γ_01 = main effect of T50,
# γ_10 = main effect of time, γ_11 = interaction

# T50:正規化
model_1 <- lmer(egfr_cr ~ 1 +
                  scaled_base_t50_fit * years +
                  (1 + years | row_id),
                data = linelist_transformed
)
summary(model_1)
confint(model_1)
tbl_regression(model_1,
               label = var_labels_tilde) 
?tbl_regression
# CKDなしの群
model_1_non_ckd <- lmer(egfr_crcys ~ 1 +
                          scaled_base_t50_fit * years +
                          (1 + years | row_id),
                        data = subset(linelist_transformed,
                                      base_cat_ckd_crcys == "non-CKD"))
summary(model_1_non_ckd)
# CKDありの群
model_1_ckd <- lmer(egfr_crcys ~ 1 +
                      scaled_base_t50_fit * years +
                      (1 + years | row_id),
                    data = subset(linelist_transformed,
                                  base_cat_ckd_crcys == "CKD"))
summary(model_2_ckd)
# Level 1 : eGFR_ij = β_0j + β_1j×Time_ij + ε_ij
# Level 2 : β_0j = γ_00 + γ_01×T50_j + γ_02×CKD_j + γ_03×(T50_j × CKD_j) + μ_0j
# 　　　　: β_1j = γ_10 + γ_11×T50_j + γ_12×CKD_j + γ_13×(T50_j × CKD_j) + μ_0j
# γ_00 = intercepts, γ_01 = main effect of T50, γ_02 = main effect of CKD、γ_03 = interaction between T50 and CKD
# γ_10 = main effect of Time, γ_11 = interaction between T50 and Time, γ_12 = interaction between CKD and Time, γ_13 = 3-way interaction
model_3 <- lmer(egfr_cr ~ 1 +
                  scaled_base_t50_fit * years * base_cat_ckd_cr +
                  (1 + years | row_id),
                data = linelist_transformed)
summary(model_3)

# model_3にベースライン年齢、性別、log_UACRを追加
model_4 <- lmer(egfr_cr ~ 1 +
                  (scaled_base_t50_fit * base_cat_ckd_cr +
                  scaled_base_age +
                  base_cat_sex +
                  scaled_log_base_uacr) * years +
                  (1 + years | row_id),
                data = linelist_transformed)
summary(model_4)

#
model_5 <- lmer(egfr_crcys ~ 1 +
                  (scaled_base_t50_fit * base_cat_ckd_crcys +
                     scaled_base_age +
                     base_cat_sex +
                     scaled_log_base_uacr +
                     base_cat_current_smoking +
                     base_cat_dm +
                     base_cat_htn +
                     scaled_base_hdl_c +
                     scaled_base_s_pi +
                     scaled_base_s_alb
                  ) * years +
                  (1 + years | row_id),
                data = linelist_transformed)
summary(model_5)
print(model_5, correlation=TRUE)
check_model(
  model_5,
  check = c("reqq", "qq", "linearity", "homogeneity", "outliers", "pp_check")
)
check_model(
  model_5,
  check = "all"
)

################################################################################
# T50とeGFR slopeとの相関
################################################################################
glimpse(linelist_long_slope)

# データの準備（linelist_transformedがデータフレームとして既に読み込まれていると仮定）
# 必要な変数の抽出と整理
plot_data <- linelist_transformed %>%
  select(base_t50_fit, slope_cr, base_cat_ckd_cr) %>%
  rename(T50 = base_t50_fit,
         eGFR_slope = slope_cr,
         CKD_status = base_cat_ckd_cr)

# 散布図の作成
p <- ggplot(plot_data, aes(x = T50, y = eGFR_slope, color = CKD_status)) +
  # 散布図
  geom_point(alpha = 0.5) +
  # 予測線と95%信頼区間
  geom_smooth(method = "lm", formula = y ~ x) +
  # テーマと色の設定
  theme_minimal() +
  scale_color_manual(values = c("non-CKD" = "#0072B2", "CKD" = "#D55E00")) +
  # ラベルの設定
  labs(
    title = "Association between T50 and eGFR Slope by CKD Status",
    x = "T50 (minutes)",
    y = "eGFR Slope (mL/min/1.73m²/year)",
    color = "CKD Status"
  ) +
  # テーマの詳細設定
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# グラフの表示
print(p)

# 層別化した回帰分析の結果を表示
fit_non_ckd <- lm(eGFR_slope ~ T50, 
                  data = subset(plot_data, CKD_status == "non-CKD"))
fit_ckd <- lm(eGFR_slope ~ T50, 
              data = subset(plot_data, CKD_status == "CKD"))

# 結果のサマリー表示
cat("\nNon-CKD group regression summary:\n")
print(summary(fit_non_ckd))
cat("\nCKD group regression summary:\n")
print(summary(fit_ckd))

# 相関係数の計算
cor_non_ckd <- with(subset(plot_data, CKD_status == "non-CKD"),
                    cor.test(T50, eGFR_slope, method = "pearson"))
cor_ckd <- with(subset(plot_data, CKD_status == "CKD"),
                cor.test(T50, eGFR_slope, method = "pearson"))

cat("\nCorrelation coefficients:\n")
cat("Non-CKD group: r =", round(cor_non_ckd$estimate, 3),
    ", p =", round(cor_non_ckd$p.value, 4), "\n")
cat("CKD group: r =", round(cor_ckd$estimate, 3),
    ", p =", round(cor_ckd$p.value, 4), "\n")