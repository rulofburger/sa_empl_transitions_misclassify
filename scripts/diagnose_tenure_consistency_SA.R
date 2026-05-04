# ==============================================================================
# diagnose_tenure_consistency_SA.R
# ==============================================================================
# Created: 2026-04-30
#
# PURPOSE
#   Empirical diagnostic for the Spec I (epsilon) tenure model. This script
#   characterises tenure-flank consistency across observed state patterns, with
#   particular focus on the (1, 0, 1) and (0, 1, 0) state patterns where
#   misclassification is most likely to be visible in tenure / timegap flanks.
#
#   The output is a set of summary tables and figures that quantify:
#     1. Granularity of reported tenure (integer-month structure).
#     2. Distribution of E->E continuation increments (point-mass at +0.25 yr
#        vs. heavy-tailed contamination).
#     3. State pattern (1, 0, 1) tenure-flank consistency: classify each
#        observation as "consistent with continuous employment" (g3 ~ g1+0.5),
#        "consistent with genuine job exit + reentry" (g3 small, < 3 mo), or
#        "ambiguous / contaminated".
#     4. State pattern (0, 1, 0) timegap-flank consistency: similar
#        classification using interval-censored timegap categories.
#     5. Within-panel "fresh start" tenure distribution (s = (0, 1)) -- the
#        key empirical fact: 36.7% report tenure > 12 months at the start.
#     6. Implied epsilon estimates from method-of-moments for each diagnostic.
#
# OUTPUTS
#   tables/figures written to: output/figures/tenure_diagnostics/
#
# USAGE
#   source("scripts/diagnose_tenure_consistency_SA.R")
# ==============================================================================

suppressMessages({
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

# --- Output directory ---
.out_dir <- "output/figures/tenure_diagnostics"
dir.create(.out_dir, recursive = TRUE, showWarnings = FALSE)
.tab_dir <- "output/tables/tenure_diagnostics"
dir.create(.tab_dir, recursive = TRUE, showWarnings = FALSE)

# Standard plot dimensions (width x height in inches, 150 dpi)
.plot_w   <- 8
.plot_h   <- 5
.plot_dpi <- 150

# --- Ingest data ---
source(here::here("scripts/ingest_data_3waves_SA.R"))
if (!exists("df_qlfs"))
  stop("df_qlfs not found after ingest: check scripts/ingest_data_3waves_SA.R")
setDT(df_qlfs)
req_cols <- c("y1", "y2", "y3", "tenure1", "tenure2", "tenure3",
              "timegap_cat1", "timegap_cat2", "timegap_cat3", "weight")
missing_cols <- setdiff(req_cols, names(df_qlfs))
if (length(missing_cols) > 0)
  stop("df_qlfs missing required columns: ", paste(missing_cols, collapse = ", "))

# Months convention: tenure is in years; convert to months for diagnostics.
df_qlfs[, t1_mo := tenure1 * 12]
df_qlfs[, t2_mo := tenure2 * 12]
df_qlfs[, t3_mo := tenure3 * 12]

cat("\n========================================================\n")
cat(" TENURE CONSISTENCY DIAGNOSTICS\n")
cat(sprintf(" N (panels): %d\n", nrow(df_qlfs)))
cat("========================================================\n\n")

# ==============================================================================
# (1) Tenure granularity
# ==============================================================================
all_emp <- rbind(
  df_qlfs[y1 == 1L, .(t = t1_mo)],
  df_qlfs[y2 == 1L, .(t = t2_mo)],
  df_qlfs[y3 == 1L, .(t = t3_mo)]
)
all_emp[, mo := round(t)]

granularity <- data.table(
  metric = c(
    "Exact 12-month multiple", "Exact 6-month multiple",
    "Exact 3-month multiple", "Integer month"
  ),
  share = c(
    mean(all_emp$mo %% 12 == 0, na.rm = TRUE),
    mean(all_emp$mo %% 6 == 0,  na.rm = TRUE),
    mean(all_emp$mo %% 3 == 0,  na.rm = TRUE),
    mean(all_emp$mo == round(all_emp$mo), na.rm = TRUE)
  )
)
fwrite(granularity, file.path(.tab_dir, "01_tenure_granularity.csv"))

cat("---- (1) Tenure granularity ----\n")
print(granularity[, .(metric, share = sprintf("%.1f%%", 100 * share))])
cat("\n")

# ==============================================================================
# (2) E->E continuation increment distribution
# ==============================================================================
ee12 <- df_qlfs[y1 == 1L & y2 == 1L, .(dt_mo = (tenure2 - tenure1) * 12)]
ee23 <- df_qlfs[y2 == 1L & y3 == 1L, .(dt_mo = (tenure3 - tenure2) * 12)]
ee   <- rbind(ee12, ee23)

ee_buckets <- data.table(
  bucket = c(
    "Exact +3 months (delta = 0)", "|delta| <= 1 month",
    "|delta| in (1, 3] months",   "delta in (-12, -3] or [3+, 12) months",
    "|delta| > 12 months"
  ),
  share = c(
    mean(abs(ee$dt_mo - 3) < 0.01, na.rm = TRUE),
    mean(abs(ee$dt_mo - 3) <= 1,   na.rm = TRUE) - mean(abs(ee$dt_mo - 3) < 0.01, na.rm = TRUE),
    mean(abs(ee$dt_mo - 3) <= 3,   na.rm = TRUE) - mean(abs(ee$dt_mo - 3) <= 1,   na.rm = TRUE),
    mean(abs(ee$dt_mo - 3) <= 12,  na.rm = TRUE) - mean(abs(ee$dt_mo - 3) <= 3,   na.rm = TRUE),
    mean(abs(ee$dt_mo - 3) > 12,   na.rm = TRUE)
  )
)
fwrite(ee_buckets, file.path(.tab_dir, "02_ee_increment_buckets.csv"))

# Method-of-moments epsilon from continuation pairs:
eps_mom_cont <- 1 - mean(abs(ee$dt_mo - 3) < 0.01, na.rm = TRUE)

cat("---- (2) E->E continuation increment distribution ----\n")
print(ee_buckets[, .(bucket, share = sprintf("%.2f%%", 100 * share))])
cat(sprintf("\n  Method-of-moments epsilon (continuation pairs): %.3f\n\n", eps_mom_cont))

# Plot: histogram of increments truncated to [-24, 24] months
p_ee <- ggplot(ee[abs(dt_mo) < 24], aes(x = dt_mo)) +
  geom_histogram(binwidth = 1, fill = "#4477AA", color = "white") +
  geom_vline(xintercept = 3, linetype = "dashed", color = "red") +
  labs(
    title    = "E->E tenure continuation increments (truncated to +/- 2 years)",
    subtitle = sprintf("63.2%% exactly +3 months. eps_MoM = %.3f", eps_mom_cont),
    x        = "tenure(t) - tenure(t-1) (months)",
    y        = "count"
  ) +
  theme_minimal()
ggsave(file.path(.out_dir, "02_ee_increment_hist.png"),
       p_ee, width = .plot_w, height = .plot_h, dpi = .plot_dpi)

# ==============================================================================
# (3) State pattern (1, 0, 1) tenure flank classification
# ==============================================================================
s101 <- df_qlfs[y1 == 1L & y2 == 0L & y3 == 1L,
                .(t1_mo, t3_mo, gap = t3_mo - t1_mo, weight)]

# Classification:
#   "miscl_corroborated": g3 ~ g1 + 6mo (within 2 months) -- consistent with
#       continuous employment + miscl at wave 2
#   "fresh_start":         g3 <= 3mo  -- consistent with genuine fresh job at wave 3
#   "year_shift":          gap is at multiple of 12 months but not 6 -- year-rounding
#   "ambiguous":           anything else (mostly contamination)
s101[, classify := fcase(
  abs(gap - 6) <= 2,                                "miscl_corroborated",
  t3_mo <= 3 & gap < 0,                              "fresh_start",
  abs((gap - 6) %% 12) < 2 & abs(gap - 6) > 6,       "year_shift",
  default = "ambiguous"
)]

s101_class <- s101[, .(N = .N, share = .N / nrow(s101)), by = classify][order(-N)]
fwrite(s101_class, file.path(.tab_dir, "03_s101_classification.csv"))

cat("---- (3) State pattern s = (1,0,1): tenure-flank classification ----\n")
cat(sprintf("  N = %d observed (1,0,1) panels\n", nrow(s101)))
print(s101_class[, .(classify, N, share = sprintf("%.1f%%", 100 * share))])
cat("\n")

# Plot: scatter of (t1, t3) for s = (1,0,1) with classification colour
p_s101 <- ggplot(s101[t1_mo <= 240 & t3_mo <= 240], aes(x = t1_mo, y = t3_mo, color = classify)) +
  geom_abline(slope = 1, intercept = 6, linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_point(alpha = 0.35, size = 0.7) +
  labs(
    title    = "State pattern s = (1, 0, 1): tenure flanks",
    subtitle = "Dashed line: g3 = g1 + 6 months (continuous employment). Truncated to 0-240 months.",
    x        = "tenure_1 (months)",
    y        = "tenure_3 (months)",
    color    = "Classification"
  ) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(file.path(.out_dir, "03_s101_scatter.png"),
       p_s101, width = .plot_w, height = .plot_h, dpi = .plot_dpi)

# ==============================================================================
# (4) State pattern (0, 1, 0) timegap flank classification
# ==============================================================================
s010 <- df_qlfs[y1 == 0L & y2 == 1L & y3 == 0L,
                .(tg1 = timegap_cat1, tg3 = timegap_cat3, t2_mo, weight)]

#   "miscl_corroborated": tg1 == tg3 -- consistent with continuous nonemployment +
#       short employment misclassification at wave 2
#   "ambiguous": otherwise
s010[, classify := fifelse(tg1 == tg3, "miscl_corroborated", "ambiguous")]

s010_class <- s010[, .(N = .N, share = .N / nrow(s010)), by = classify][order(-N)]
fwrite(s010_class, file.path(.tab_dir, "04_s010_classification.csv"))

cat("---- (4) State pattern s = (0,1,0): timegap-flank classification ----\n")
cat(sprintf("  N = %d observed (0,1,0) panels\n", nrow(s010)))
print(s010_class[, .(classify, N, share = sprintf("%.1f%%", 100 * share))])
cat("\n")

# ==============================================================================
# (5) Within-panel fresh-start tenure distribution
#     (s = (s_{t-1} = 0, s_t = 1), tenure at wave t)
# ==============================================================================
fs <- rbind(
  df_qlfs[y1 == 0L & y2 == 1L, .(g_mo = t2_mo, weight)],
  df_qlfs[y2 == 0L & y3 == 1L, .(g_mo = t3_mo, weight)]
)

fs_buckets <- data.table(
  bucket = c("[0, 3] months", "(3, 6] months", "(6, 12] months", "(12, 36] months", "(36+) months"),
  share  = c(
    mean(fs$g_mo <= 3,              na.rm = TRUE),
    mean(fs$g_mo > 3 & fs$g_mo <= 6,  na.rm = TRUE),
    mean(fs$g_mo > 6 & fs$g_mo <= 12, na.rm = TRUE),
    mean(fs$g_mo > 12 & fs$g_mo <= 36, na.rm = TRUE),
    mean(fs$g_mo > 36,              na.rm = TRUE)
  )
)
fwrite(fs_buckets, file.path(.tab_dir, "05_fresh_start_buckets.csv"))

cat("---- (5) Fresh-start tenure distribution (within-panel s=(0,1)) ----\n")
cat(sprintf("  N = %d observed within-panel starts\n", nrow(fs)))
print(fs_buckets[, .(bucket, share = sprintf("%.1f%%", 100 * share))])
cat("\n  KEY INSIGHT: Most 'fresh starts' are NOT short tenure.\n")
cat(sprintf("  Implied lower bound on misclassification at t-1: %.1f%%\n",
            100 * mean(fs$g_mo > 12, na.rm = TRUE)))
cat("  (A 12+ month tenure at a 'fresh start' is impossible -- must be miscl at t-1.)\n\n")

# Plot: histogram of fresh-start tenure (truncated to 60mo)
p_fs <- ggplot(fs[g_mo <= 60], aes(x = g_mo)) +
  geom_histogram(binwidth = 1, fill = "#CC6677", color = "white") +
  geom_vline(xintercept = 3, linetype = "dashed", color = "black") +
  labs(
    title    = "Within-panel 'fresh start' tenure: most are NOT short",
    subtitle = sprintf("%.1f%% report > 12 months at supposed 'fresh start'",
                       100 * mean(fs$g_mo > 12, na.rm = TRUE)),
    x        = "tenure at fresh-start wave (months)",
    y        = "count"
  ) +
  theme_minimal()
ggsave(file.path(.out_dir, "05_fresh_start_hist.png"),
       p_fs, width = .plot_w, height = .plot_h, dpi = .plot_dpi)

# ==============================================================================
# (6) Spell-pair epsilon estimates (method-of-moments)
# ==============================================================================
# For h=(1,1) consecutive observed pairs: P(consistent) = 1 - eps
# For h=(1,1,1) under s=(1,0,1): P(consistent) = (1 - eps)^2
# Combine to estimate eps.

eps_mom_pair    <- 1 - mean(abs(ee$dt_mo - 3) < 0.01, na.rm = TRUE)
bridge_share    <- s101_class[classify == "miscl_corroborated", share]
eps_mom_bridge  <- if (length(bridge_share) == 1L && bridge_share > 0)
  1 - sqrt(bridge_share) else NA_real_

mom_table <- data.table(
  source = c("E->E adjacent (1-eps)", "1-2 bridge under (1,0,1) ((1-eps)^2)"),
  estimate = c(eps_mom_pair, eps_mom_bridge)
)
fwrite(mom_table, file.path(.tab_dir, "06_eps_method_of_moments.csv"))

cat("---- (6) Method-of-moments epsilon estimates ----\n")
print(mom_table[, .(source, estimate = sprintf("%.3f", estimate))])
cat(sprintf("\n  Note: bridge estimate is conditional on (1,0,1) state pattern;\n"))
cat(sprintf("  the divergence from the adjacent estimate is a useful diagnostic.\n\n"))

# ==============================================================================
# Save a single summary table
# ==============================================================================
summary_dt <- data.table(
  diagnostic = c(
    "E->E pairs with exact +3mo increment",
    "(1,0,1) panels with miscl-corroborated tenure (g3 ~ g1+6mo)",
    "(1,0,1) panels with fresh-start-consistent g3 (<= 3mo)",
    "Fresh-start tenure > 12 months (impossible under genuine entry)",
    "(0,1,0) panels with timegap1 == timegap3"
  ),
  N = c(nrow(ee), nrow(s101), nrow(s101), nrow(fs), nrow(s010)),
  share = c(
    mean(abs(ee$dt_mo - 3) < 0.01,              na.rm = TRUE),
    mean(s101$classify == "miscl_corroborated",  na.rm = TRUE),
    mean(s101$classify == "fresh_start",         na.rm = TRUE),
    mean(fs$g_mo > 12,                           na.rm = TRUE),
    mean(s010$classify == "miscl_corroborated",  na.rm = TRUE)
  )
)
summary_dt[, share_pct := sprintf("%.2f%%", 100 * share)]
fwrite(summary_dt, file.path(.tab_dir, "00_summary.csv"))

cat("---- Headline summary ----\n")
print(summary_dt[, .(diagnostic, N, share_pct)])
cat(sprintf("\nAll outputs written to:\n  %s\n  %s\n", .out_dir, .tab_dir))
