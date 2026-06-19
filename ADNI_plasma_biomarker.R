# ADNI_plasma_biomarker.R
# Plasma p-tau217 and Aβ42/40 classification of amyloid-PET status in ADNI.
# Includes participant-grouped cross-validation, GEE/GLMM repeated-measures
# sensitivity analyses, clinical threshold sensitivity, temporal-window
# sensitivity, collinearity diagnostics, and interaction-model checks.

suppressPackageStartupMessages({
  library(tidyverse)
  library(tidymodels)
  library(janitor)
  library(lubridate)
  library(pROC)
  library(broom)
})

set.seed(2026)

# -------------------------------------------------------------
# 0. Tunable parameters
# -------------------------------------------------------------
GAP_DAYS_MAX <- 365    # max |plasma draw date - PET scan date|, in days
V_FOLDS      <- 5      # number of cross-validation folds
BOOT_N       <- 2000   # bootstrap resamples for AUC CIs

# Optional sensitivity windows for plasma-PET temporal matching.
GAP_SENSITIVITY_WINDOWS <- c(30, 90, 180, 365)

# Sensitivity-targeted thresholds for the combined model.
SENSITIVITY_TARGETS <- c(0.90, 0.95)

# --- Figure display options (COSMETIC ONLY -- do not affect the analysis) ---
LEGACY_AXES <- TRUE    # TRUE  => ROC x-axis as "1 - Specificity" (0->1)
                       # FALSE => pROC default (Specificity, 1.0 -> 0.0)

# NOTE: confirm the p-tau217 unit against the UPENN plasma data dictionary
# before publication-level use. Fujirebio Lumipulse plasma p-tau217 is typically
# reported in pg/mL, but verify it for this exact ADNI file.
PTAU_UNIT <- "pg/mL"

# Readable axis labels + consistent colors used by all figures below.
lab_ptau   <- bquote("Plasma p-tau217 (" * .(PTAU_UNIT) * ")")
lab_ratio  <- expression(paste("Plasma A", beta, "42/40 ratio"))
lab_status <- "Amyloid-PET status"
fill_cols  <- c(Negative = "#F8766D", Positive = "#00BFC4")

# -------------------------------------------------------------
# 1. Load raw ADNI files
# -------------------------------------------------------------
# This helper allows either the clean repository file names or local duplicate
# names such as "..._31Mar2026(2).csv".
find_csv <- function(primary_name, pattern) {
  if (file.exists(primary_name)) return(primary_name)
  hits <- list.files(pattern = pattern)
  if (length(hits) == 0) {
    stop("Could not find required CSV. Expected ", primary_name,
         " or a file matching pattern: ", pattern)
  }
  hits[[1]]
}

plasma_path <- find_csv(
  "UPENN_PLASMA_FUJIREBIO_QUANTERIX_31Mar2026.csv",
  "^UPENN_PLASMA_FUJIREBIO_QUANTERIX_31Mar2026.*\\.csv$"
)
pet_path <- find_csv(
  "UCBERKELEY_AMY_6MM_31Mar2026.csv",
  "^UCBERKELEY_AMY_6MM_31Mar2026.*\\.csv$"
)
dx_path <- find_csv(
  "DXSUM_31Mar2026.csv",
  "^DXSUM_31Mar2026.*\\.csv$"
)

cat("Using files:\n")
cat("  Plasma:", plasma_path, "\n")
cat("  Amyloid-PET:", pet_path, "\n")
cat("  Diagnosis:", dx_path, "\n\n")

plasma <- read_csv(plasma_path, show_col_types = FALSE) %>% clean_names()
pet    <- read_csv(pet_path,    show_col_types = FALSE) %>% clean_names()
dx     <- read_csv(dx_path,     show_col_types = FALSE) %>% clean_names()

# -------------------------------------------------------------
# 2. Clean each table; one row per (rid, viscode2)
#    clean_names() renders the plasma columns as
#    p_t217_f, ab42_ab40_f, nf_l_q, gfap_q
# -------------------------------------------------------------
plasma_clean <- plasma %>%
  transmute(
    rid, viscode2,
    examdate = ymd(examdate),
    ptau217  = p_t217_f,
    ab_ratio = ab42_ab40_f,
    nfl      = nf_l_q,
    gfap     = gfap_q
  ) %>%
  mutate(across(c(ptau217, ab_ratio, nfl, gfap),
                ~ ifelse(.x < 0, NA, .x))) %>%   # below-detection codes -> NA
  drop_na(ptau217, ab_ratio) %>%
  arrange(rid, viscode2, examdate) %>%
  distinct(rid, viscode2, .keep_all = TRUE)

pet_clean <- pet %>%
  transmute(
    rid, viscode2,
    scandate = ymd(scandate),
    amyloid_status, centiloids, summary_suvr
  ) %>%
  filter(!is.na(amyloid_status)) %>%
  arrange(rid, viscode2, scandate) %>%
  distinct(rid, viscode2, .keep_all = TRUE)

dx_clean <- dx %>%
  transmute(rid, viscode2, diagnosis) %>%
  arrange(rid, viscode2) %>%
  distinct(rid, viscode2, .keep_all = TRUE)

# -------------------------------------------------------------
# 3. Merge on participant + visit; keep an unfiltered merged dataset
#    so temporal-window sensitivity analyses can be run later.
# -------------------------------------------------------------
dat_merged <- plasma_clean %>%
  inner_join(pet_clean, by = c("rid", "viscode2")) %>%
  left_join(dx_clean,  by = c("rid", "viscode2")) %>%
  mutate(gap_days = abs(as.numeric(examdate - scandate)))

n_before <- nrow(dat_merged)

dat <- dat_merged %>%
  filter(gap_days <= GAP_DAYS_MAX)

cat("Visits before temporal filter:", n_before, "\n")
cat("Visits after  temporal filter (<=", GAP_DAYS_MAX, "days):", nrow(dat), "\n")
cat("Median |plasma-PET| gap (days):", median(dat$gap_days), "\n\n")

# -------------------------------------------------------------
# 4. Recode labels
# -------------------------------------------------------------
recode_analysis_labels <- function(x) {
  x %>%
    mutate(
      amyloid_status = factor(if_else(amyloid_status == 1, "Positive", "Negative"),
                              levels = c("Negative", "Positive")),
      diagnosis = factor(case_when(
        diagnosis == 1 ~ "CN",
        diagnosis == 2 ~ "MCI",
        diagnosis == 3 ~ "AD",
        TRUE ~ NA_character_
      ), levels = c("CN", "MCI", "AD"))
    )
}

dat <- recode_analysis_labels(dat)

# Non-independence check: how many participants contribute >1 visit?
mv <- dat %>% count(rid) %>%
  summarise(n_participants = n(), n_rows = sum(n),
            pct_multivisit = round(mean(n > 1) * 100, 1))
cat("Unique participants:", mv$n_participants,
    "| rows:", mv$n_rows,
    "| % with >1 visit:", mv$pct_multivisit, "%\n")
print(table(dat$amyloid_status)); cat("\n")

# -------------------------------------------------------------
# 5. Participant-grouped 5-fold CV
#    This prevents within-subject leakage between training and testing.
#    The SAME folds are reused for all three models so out-of-fold
#    predictions are paired on identical test cases.
# -------------------------------------------------------------
set.seed(2026)
folds <- group_vfold_cv(dat, group = rid, v = V_FOLDS)

glm_spec <- logistic_reg() %>% set_engine("glm")

make_recipes <- function(data_input) {
  list(
    both = recipe(amyloid_status ~ ptau217 + ab_ratio, data = data_input) %>%
      step_log(ptau217, offset = 0.001) %>%
      step_normalize(all_numeric_predictors()),
    ptau = recipe(amyloid_status ~ ptau217, data = data_input) %>%
      step_log(ptau217, offset = 0.001) %>%
      step_normalize(all_numeric_predictors()),
    ab = recipe(amyloid_status ~ ab_ratio, data = data_input) %>%
      step_normalize(all_numeric_predictors())
  )
}

recipes_main <- make_recipes(dat)
rec_both <- recipes_main$both
rec_ptau <- recipes_main$ptau
rec_ab   <- recipes_main$ab

ctrl <- control_resamples(save_pred = TRUE)

run_cv <- function(rec, folds_input = folds) {
  workflow() %>% add_recipe(rec) %>% add_model(glm_spec) %>%
    fit_resamples(folds_input, metrics = metric_set(roc_auc), control = ctrl)
}

res_both <- run_cv(rec_both)
res_ptau <- run_cv(rec_ptau)
res_ab   <- run_cv(rec_ab)

# Pooled out-of-fold predictions.
# Each visit-pair is predicted once by a model that did not train on that participant.
# All visits from the same participant are kept in the same fold.
pred_both <- collect_predictions(res_both) %>% arrange(.row)
pred_ptau <- collect_predictions(res_ptau) %>% arrange(.row)
pred_ab   <- collect_predictions(res_ab)   %>% arrange(.row)

# -------------------------------------------------------------
# 6. Pooled CV-AUC with bootstrap 95% CIs
# -------------------------------------------------------------
roc_from <- function(p) roc(p$amyloid_status, p$.pred_Positive,
                            levels = c("Negative", "Positive"),
                            direction = "<", quiet = TRUE)
r_ptau <- roc_from(pred_ptau)
r_ab   <- roc_from(pred_ab)
r_both <- roc_from(pred_both)

auc_row <- function(r, label) {
  set.seed(2026)
  ci <- ci.auc(r, method = "bootstrap", boot.n = BOOT_N)
  tibble(model = label, cv_auc = as.numeric(auc(r)),
         lower95 = ci[1], upper95 = ci[3])
}
auc_table <- bind_rows(
  auc_row(r_ptau, "log(p-tau217) only"),
  auc_row(r_ab,   "Abeta42/40 only"),
  auc_row(r_both, "Combined")
) %>% mutate(across(c(cv_auc, lower95, upper95), ~ round(.x, 3)))

cat("===== Pooled CV-AUC with bootstrap 95% CI =====\n")
print(auc_table)

# -------------------------------------------------------------
# 7. DeLong: does Abeta42/40 add discrimination beyond p-tau217?
#    Paired on identical pooled out-of-fold visit-pair predictions.
#    Because some participants contribute repeated visits, treat this p-value
#    as supportive rather than fully participant-clustered inference.
# -------------------------------------------------------------
dl <- roc.test(r_ptau, r_both, method = "delong", paired = TRUE)
cat(sprintf("\n===== DeLong: p-tau217 alone vs combined =====\n"))
cat(sprintf("AUC %.3f vs %.3f | Z = %.3f | p = %.4g\n",
            as.numeric(auc(r_ptau)), as.numeric(auc(r_both)),
            dl$statistic, dl$p.value))
cat(if (dl$p.value < 0.05)
      "-> Abeta42/40 adds statistically significant discrimination beyond p-tau217.\n"
    else
      "-> Abeta42/40 adds no statistically significant discrimination beyond p-tau217.\n")

# -------------------------------------------------------------
# 8. Pooled operating points for combined model
#    0.5 is a default ML threshold only. Youden and sensitivity-targeted
#    thresholds are exploratory and require external validation.
# -------------------------------------------------------------
threshold_metrics <- function(pred, threshold, threshold_type) {
  tmp <- pred %>%
    mutate(
      .pred_class = factor(
        if_else(.pred_Positive >= threshold, "Positive", "Negative"),
        levels = c("Negative", "Positive")
      )
    )

  bind_rows(
    accuracy(tmp, amyloid_status, .pred_class),
    sens(tmp, amyloid_status, .pred_class, event_level = "second"),
    spec(tmp, amyloid_status, .pred_class, event_level = "second"),
    ppv(tmp, amyloid_status, .pred_class, event_level = "second"),
    npv(tmp, amyloid_status, .pred_class, event_level = "second")
  ) %>%
    mutate(threshold_type = threshold_type,
           threshold = threshold,
           .before = .metric) %>%
    select(threshold_type, threshold, .metric, .estimate)
}

coords_tbl <- function(...) {
  as_tibble(coords(..., transpose = FALSE))
}

youden_tbl <- coords_tbl(
  r_both,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity")
)
youden_threshold <- as.numeric(youden_tbl$threshold[1])

sens_thresholds <- map_dfr(SENSITIVITY_TARGETS, function(target) {
  out <- coords_tbl(
    r_both,
    x = target,
    input = "sensitivity",
    ret = c("threshold", "sensitivity", "specificity")
  ) %>% slice(1)
  out %>% mutate(target_sensitivity = target, .before = 1)
})

threshold_table <- bind_rows(
  threshold_metrics(pred_both, 0.5, "0.5 default"),
  threshold_metrics(pred_both, youden_threshold, "Youden"),
  map2_dfr(
    sens_thresholds$threshold,
    sens_thresholds$target_sensitivity,
    ~ threshold_metrics(pred_both, as.numeric(.x), paste0("Sensitivity ", .y))
  )
) %>%
  mutate(across(c(threshold, .estimate), ~ round(.x, 3)))

cat("\n===== Combined model threshold sensitivity =====\n")
cat("Youden threshold:\n")
print(youden_tbl)
cat("\nSensitivity-targeted thresholds:\n")
print(sens_thresholds)
cat("\nThreshold comparison table:\n")
print(threshold_table)

# Retain the original 0.5 operating point and confusion matrix for continuity.
op <- pred_both %>%
  mutate(.pred_class = factor(if_else(.pred_Positive >= 0.5, "Positive", "Negative"),
                              levels = c("Negative", "Positive")))
cat("\nConfusion matrix at 0.5 threshold (pooled):\n")
print(op %>% conf_mat(amyloid_status, .pred_class))

# -------------------------------------------------------------
# 9. Final model coefficients and repeated-measures sensitivity
# -------------------------------------------------------------
final_fit <- workflow() %>% add_recipe(rec_both) %>% add_model(glm_spec) %>% fit(dat)
cat("\nFinal combined-model coefficients (log + standardized scale):\n")
print(broom::tidy(extract_fit_parsnip(final_fit)$fit))

# Cluster-aware coefficient inference. These sensitivity models are intended for
# inference about coefficients, not as the main out-of-fold prediction model.
dat_model <- dat %>%
  mutate(
    y = as.integer(amyloid_status == "Positive"),
    log_ptau217 = log(ptau217 + 0.001),
    log_ptau217_z = as.numeric(scale(log_ptau217)),
    ab_ratio_z = as.numeric(scale(ab_ratio))
  )

# -------------------------------------------------------------
# 9a. Multicollinearity check: p-tau217 and Abeta42/40
#     With two predictors, a one-line VIF is enough to show whether the
#     combined model is unstable due to collinearity.
# -------------------------------------------------------------
cat("\n===== Multicollinearity check: log(p-tau217) and Abeta42/40 =====\n")

cor_pearson <- cor.test(dat_model$log_ptau217_z, dat_model$ab_ratio_z, method = "pearson")
cor_spearman <- cor.test(dat_model$log_ptau217_z, dat_model$ab_ratio_z, method = "spearman")

cat(sprintf("Pearson r = %.3f | p = %.4g\n", cor_pearson$estimate, cor_pearson$p.value))
cat(sprintf("Spearman rho = %.3f | p = %.4g\n", cor_spearman$estimate, cor_spearman$p.value))

# Fit an additive glm on the same standardized scale used for coefficient checks.
additive_glm <- glm(
  y ~ log_ptau217_z + ab_ratio_z,
  data = dat_model,
  family = binomial
)

# Prefer performance::check_collinearity() if available; otherwise compute the
# equivalent VIF manually. With two predictors, both variables have the same VIF.
if (requireNamespace("performance", quietly = TRUE)) {
  collinearity_check <- performance::check_collinearity(additive_glm)
  print(collinearity_check)

  collinearity_table <- as_tibble(collinearity_check) %>%
    select(any_of(c("Parameter", "VIF", "VIF_CI_low", "VIF_CI_high",
                    "Tolerance", "Tolerance_CI_low", "Tolerance_CI_high")))
} else {
  vif_r2 <- summary(lm(log_ptau217_z ~ ab_ratio_z, data = dat_model))$r.squared
  vif_manual <- 1 / (1 - vif_r2)
  cat(sprintf("Manual VIF for both predictors = %.3f\n", vif_manual))
  cat("Package 'performance' is not installed; install with install.packages('performance') to run check_collinearity().\n")

  collinearity_table <- tibble(
    Parameter = c("log_ptau217_z", "ab_ratio_z"),
    VIF = vif_manual,
    VIF_CI_low = NA_real_,
    VIF_CI_high = NA_real_,
    Tolerance = 1 / vif_manual,
    Tolerance_CI_low = NA_real_,
    Tolerance_CI_high = NA_real_
  )
}

collinearity_summary <- tibble(
  statistic = c("pearson_r", "pearson_p", "spearman_rho", "spearman_p"),
  value = c(
    unname(cor_pearson$estimate),
    cor_pearson$p.value,
    unname(cor_spearman$estimate),
    cor_spearman$p.value
  )
)

# -------------------------------------------------------------
# 9b. Interaction check: log(p-tau217) * Abeta42/40
#     Tests whether the interaction adds meaningful information beyond the
#     additive combined model, using both LRT and paired out-of-fold ΔAUC.
# -------------------------------------------------------------
cat("\n===== Interaction check: log(p-tau217) * Abeta42/40 =====\n")

interaction_glm <- glm(
  y ~ log_ptau217_z * ab_ratio_z,
  data = dat_model,
  family = binomial
)

lrt_tbl <- anova(additive_glm, interaction_glm, test = "LRT")
cat("\nLikelihood-ratio test: additive vs interaction glm\n")
print(lrt_tbl)

rec_int <- recipe(amyloid_status ~ ptau217 + ab_ratio, data = dat) %>%
  step_log(ptau217, offset = 0.001) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_interact(terms = ~ ptau217:ab_ratio)

res_int <- run_cv(rec_int)
pred_int <- collect_predictions(res_int) %>% arrange(.row)
r_int <- roc_from(pred_int)
dl_int <- roc.test(r_both, r_int, method = "delong", paired = TRUE)

interaction_comparison_table <- tibble(
  comparison = "Additive combined vs interaction model",
  auc_additive = as.numeric(auc(r_both)),
  auc_interaction = as.numeric(auc(r_int)),
  delta_auc = auc_interaction - auc_additive,
  delong_z = as.numeric(dl_int$statistic),
  delong_p = dl_int$p.value,
  lrt_p = lrt_tbl$`Pr(>Chi)`[2]
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 6)))

cat("\nInteraction model CV-AUC comparison:\n")
print(interaction_comparison_table)

cat(if (interaction_comparison_table$delong_p < 0.05)
      "-> Interaction model improves AUC versus the additive combined model by DeLong test.\n"
    else
      "-> Interaction model does not provide a statistically significant AUC improvement versus the additive combined model.\n")

cat("\n===== Repeated-measures sensitivity: GEE / GLMM =====\n")

if (requireNamespace("geepack", quietly = TRUE)) {
  gee_fit <- geepack::geeglm(
    y ~ log_ptau217_z + ab_ratio_z,
    id = rid,
    data = dat_model,
    family = binomial,
    corstr = "exchangeable"
  )
  cat("\nGEE model clustered by participant ID:\n")
  print(summary(gee_fit))
} else {
  cat("\nPackage 'geepack' is not installed; skipping GEE sensitivity model.\n")
  cat("Install with install.packages('geepack') to run this section.\n")
}

if (requireNamespace("lme4", quietly = TRUE)) {
  glmm_fit <- lme4::glmer(
    y ~ log_ptau217_z + ab_ratio_z + (1 | rid),
    data = dat_model,
    family = binomial,
    control = lme4::glmerControl(optimizer = "bobyqa")
  )
  cat("\nGLMM model with participant random intercept:\n")
  print(summary(glmm_fit))
} else {
  cat("\nPackage 'lme4' is not installed; skipping GLMM sensitivity model.\n")
  cat("Install with install.packages('lme4') to run this section.\n")
}

# -------------------------------------------------------------
# 10. Temporal-window sensitivity analysis
#     Re-runs the full grouped-CV workflow using stricter plasma-PET windows.
# -------------------------------------------------------------
run_auc_by_gap <- function(max_gap) {
  dat_sub <- dat_merged %>%
    filter(gap_days <= max_gap) %>%
    recode_analysis_labels()

  n_rows <- nrow(dat_sub)
  n_participants <- n_distinct(dat_sub$rid)

  if (n_rows < 50 || n_distinct(dat_sub$amyloid_status) < 2 || n_participants < V_FOLDS) {
    return(tibble(
      max_gap_days = max_gap,
      n_rows = n_rows,
      n_participants = n_participants,
      median_gap_days = ifelse(n_rows > 0, median(dat_sub$gap_days), NA_real_),
      auc_ptau = NA_real_,
      auc_ab = NA_real_,
      auc_both = NA_real_
    ))
  }

  set.seed(2026)
  folds_sub <- group_vfold_cv(dat_sub, group = rid, v = V_FOLDS)
  recs_sub <- make_recipes(dat_sub)

  pred_ptau_sub <- collect_predictions(run_cv(recs_sub$ptau, folds_sub)) %>% arrange(.row)
  pred_ab_sub   <- collect_predictions(run_cv(recs_sub$ab,   folds_sub)) %>% arrange(.row)
  pred_both_sub <- collect_predictions(run_cv(recs_sub$both, folds_sub)) %>% arrange(.row)

  r_ptau_sub <- roc_from(pred_ptau_sub)
  r_ab_sub   <- roc_from(pred_ab_sub)
  r_both_sub <- roc_from(pred_both_sub)

  tibble(
    max_gap_days = max_gap,
    n_rows = n_rows,
    n_participants = n_participants,
    median_gap_days = median(dat_sub$gap_days),
    auc_ptau = as.numeric(auc(r_ptau_sub)),
    auc_ab = as.numeric(auc(r_ab_sub)),
    auc_both = as.numeric(auc(r_both_sub))
  )
}

gap_sensitivity <- map_dfr(GAP_SENSITIVITY_WINDOWS, run_auc_by_gap) %>%
  mutate(across(starts_with("auc"), ~ round(.x, 3)))

cat("\n===== Temporal-window sensitivity analysis =====\n")
print(gap_sensitivity)

# -------------------------------------------------------------
# 11. Model-comparison ROC figure (pooled out-of-fold curves)
# -------------------------------------------------------------
fmt_ci <- function(i) sprintf(": AUC %.3f (%.3f\u2013%.3f)",
                              auc_table$cv_auc[i], auc_table$lower95[i],
                              auc_table$upper95[i])

leg_list <- vector("list", nrow(auc_table) + 1L)
for (i in seq_len(nrow(auc_table))) {
  m <- auc_table$model[i]
  leg_list[[i]] <- if (grepl("Abeta", m))
    bquote(A * beta * "42/40 only" * .(fmt_ci(i)))
  else
    bquote(.(m) * .(fmt_ci(i)))
}
leg_list[[nrow(auc_table) + 1L]] <-
  bquote("Paired DeLong: p-tau217 vs combined, p = " *
           .(formatC(dl$p.value, format = "e", digits = 1)))
leg_expr <- do.call(expression, leg_list)

png("model_comparison_roc.png", width = 7, height = 6, units = "in", res = 300)
plot(r_ptau, col = "blue", lwd = 2, legacy.axes = LEGACY_AXES,
     main = "Model Comparison ROC (5-fold grouped CV, pooled)",
     xlab = if (LEGACY_AXES) "1 - Specificity (false-positive rate)" else "Specificity",
     ylab = "Sensitivity (true-positive rate)")
plot(r_ab,   col = "red",       add = TRUE, lwd = 2)
plot(r_both, col = "darkgreen", add = TRUE, lwd = 2)
legend("bottomright", legend = leg_expr,
       col = c("blue", "red", "darkgreen", NA),
       lwd = c(2, 2, 2, NA), seg.len = 1.4, cex = 0.85, bty = "n")
dev.off()
cat("\nSaved: model_comparison_roc.png\n")

# -------------------------------------------------------------
# 12. Descriptive figures
# -------------------------------------------------------------
p_box <- ggplot(dat, aes(amyloid_status, ptau217, fill = amyloid_status)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0.4) +
  scale_fill_manual(values = fill_cols) +
  labs(x = lab_status, y = lab_ptau) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
# Optional: de-compress the outliers with a log y-axis
# p_box <- p_box + scale_y_log10()
ggsave("ptau217_boxplot.png", p_box, width = 6, height = 6, dpi = 300)
cat("Saved: ptau217_boxplot.png\n")

p_sc <- ggplot(dat, aes(ptau217, ab_ratio, color = amyloid_status)) +
  geom_point(alpha = 0.5, size = 1.6) +
  scale_color_manual(values = fill_cols, name = lab_status) +
  labs(x = lab_ptau, y = lab_ratio) +
  theme_minimal(base_size = 14)
# NOTE: a single Abeta42/40 point near ~0.35 is a known QC-flagged extreme value.
# It is retained for transparency. To zoom the display without deleting data:
# p_sc <- p_sc + coord_cartesian(ylim = c(0, 0.16))
ggsave("biomarker_scatter.png", p_sc, width = 7.5, height = 6, dpi = 300)
cat("Saved: biomarker_scatter.png\n")

# Optional: write sensitivity tables to CSV for easier README updating.
write_csv(auc_table, "auc_table.csv")
write_csv(threshold_table, "threshold_sensitivity_table.csv")
write_csv(gap_sensitivity, "temporal_window_sensitivity_table.csv")
write_csv(collinearity_summary, "biomarker_correlation_table.csv")
write_csv(collinearity_table, "collinearity_vif_table.csv")
write_csv(interaction_comparison_table, "interaction_comparison_table.csv")
cat("Saved: auc_table.csv, threshold_sensitivity_table.csv, temporal_window_sensitivity_table.csv, biomarker_correlation_table.csv, collinearity_vif_table.csv, interaction_comparison_table.csv\n")
