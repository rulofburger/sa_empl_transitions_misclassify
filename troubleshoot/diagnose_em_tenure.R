# ==============================================================================
# EM-tenure: Diagnostic script
# ==============================================================================
# Companion to troubleshoot/DIAGNOSIS.md. Run from the project root:
#   source("troubleshoot/diagnose_em_tenure.R")
#
# Structure:
#   §1  Synthetic validation    — confirms the EM code is correct
#   §2  Real data diagnostics   — quantifies data-model mismatches
#   §3  Progressive stress tests — isolates which pathology breaks estimation
#   §4  Targeted re-estimation  — applies remedies to real QLFS data
#
# NOTE: §2 and §4 require the QLFS data. §1 and §3 run on synthetic data only
# and can be executed without data access.
# ==============================================================================

library(data.table)
library(ggplot2)

# Load EM-tenure module (all functions become available in the current env)
source("EM-tenure/R/source_all.R")


# ==============================================================================
# HELPERS
# ==============================================================================

# Print a compact parameter comparison table to the console.
# true_params: named list — must contain the parameters present in est_params
# est_params:  fit$params from em_fit_tenure()
# label:       string identifying this run
# NOTE: sigma2_d is not estimated in discrete_timegap mode (the default).
#   The nonemployment duration is modelled as interval-censored Exp(lambda_d),
#   so only lambda_d (derived deterministically from theta0) appears in params.
#   Including sigma2_d in the comparison would cause vector recycling warnings.
.print_recovery <- function(true_params, est_params, label) {
  candidates <- c("alpha", "theta1", "theta0", "pi", "sigma2_g", "lambda_d", "sigma2_d")
  # Only compare parameters present in BOTH the truth list and the fit output
  params_of_interest <- intersect(candidates,
                                  intersect(names(true_params), names(est_params)))
  truth <- unlist(true_params[params_of_interest])
  est   <- unlist(est_params[params_of_interest])
  abs_err <- abs(est - truth)
  rel_err <- abs_err / (abs(truth) + 1e-12)

  dt <- data.table(
    param    = params_of_interest,
    true     = round(truth,   4),
    estimate = round(est,     4),
    abs_err  = round(abs_err, 4),
    rel_err  = round(rel_err, 3)
  )
  cat("\n", strrep("=", 60), "\n", label, "\n", strrep("=", 60), "\n", sep = "")
  print(dt)
  invisible(dt)
}

# Bucket continuous duration values into QLFS-style categories, then assign
# midpoints (in years). Mirrors the ingest script's timegap mapping.
.bucket_timegap <- function(x_years) {
  # Boundaries and midpoints in years (converted from months)
  breaks    <- c(0, 3, 6, 9, 12, 36, 60, Inf) / 12
  midpoints <- c(0, 1.5, 4.5, 7.5, 10.5, 24, 90) / 12  # 7 intervals

  idx <- findInterval(x_years, breaks, rightmost.closed = FALSE)
  idx <- pmin(pmax(idx, 1L), length(midpoints))
  midpoints[idx]
}


# ==============================================================================
# §1  SYNTHETIC VALIDATION
# ==============================================================================
# Purpose: confirm that em_fit_tenure() correctly recovers known true parameters
# on well-behaved synthetic data. If this fails, there is a code bug.
# If this succeeds, estimation problems on real data are data/misspecification
# issues, not code errors.
# ==============================================================================

cat("\n", strrep("#", 70), "\n",
    "§1  SYNTHETIC VALIDATION\n",
    strrep("#", 70), "\n", sep = "")

# True parameters (representative of South Africa; θ small, π small)
# sigma2_d is omitted: default discrete_timegap = TRUE uses interval-censored
# Exp(lambda_d) for nonemployment durations. lambda_d is derived from theta0
# via the CTMC link: lambda_d = -log(1 - theta0) / delta (delta = 3 quarters).
true_params <- list(
  alpha    = 0.55,
  theta1   = 0.97,
  theta0   = 0.03,
  pi       = 0.03,
  sigma2_g = 0.005,
  lambda_d = ctmc_lambda_from_theta(0.03)  # ≈ 0.01014; no sigma2_d in discrete mode
)

cat("\nTrue parameters:\n")
print(unlist(true_params))

# --- 1a. Small sample (n=500) ------------------------------------------------
cat("\n--- 1a. n = 500, with misclassification ---\n")
set.seed(42)
df_synth_500 <- simulate_panel(
  n        = 500,
  alpha    = true_params$alpha,
  theta1   = true_params$theta1,
  theta0   = true_params$theta0,
  pi       = true_params$pi,
  sigma2_g = true_params$sigma2_g,
  # sigma2_d omitted: discrete_timegap = TRUE (default) sets sigma_d = NA
  seed     = 42
)

fit_500 <- em_fit_tenure(df_synth_500, misclassification = TRUE, verbose = 1L)
.print_recovery(true_params, fit_500$params, "n=500, with misclassification")
cat("Converged:", fit_500$converged, "| Iterations:", fit_500$iterations, "\n")

# --- 1b. Larger sample (n=5000) ----------------------------------------------
cat("\n--- 1b. n = 5000, with misclassification ---\n")
set.seed(123)
df_synth_5000 <- simulate_panel(
  n        = 5000,
  alpha    = true_params$alpha,
  theta1   = true_params$theta1,
  theta0   = true_params$theta0,
  pi       = true_params$pi,
  sigma2_g = true_params$sigma2_g,
  # sigma2_d omitted: discrete_timegap = TRUE (default) sets sigma_d = NA
  seed     = 123
)

fit_5000 <- em_fit_tenure(df_synth_5000, misclassification = TRUE, verbose = 1L)
.print_recovery(true_params, fit_5000$params, "n=5000, with misclassification")
.print_recovery(true_params, fit_500$params, "n=500, with misclassification")
cat("Converged:", fit_5000$converged, "| Iterations:", fit_5000$iterations, "\n")

# --- 1c. No misclassification ------------------------------------------------
cat("\n--- 1c. n = 5000, without misclassification ---\n")
true_params_no_pi <- modifyList(true_params, list(pi = 0))
set.seed(456)
df_synth_nopi <- simulate_panel(
  n        = 5000,
  alpha    = true_params_no_pi$alpha,
  theta1   = true_params_no_pi$theta1,
  theta0   = true_params_no_pi$theta0,
  pi       = 0,
  sigma2_g = true_params_no_pi$sigma2_g,
  # sigma2_d omitted: discrete_timegap = TRUE (default) sets sigma_d = NA
  seed     = 456
)

fit_nopi <- em_fit_tenure(df_synth_nopi, misclassification = FALSE, verbose = 1L)
.print_recovery(true_params_no_pi, fit_nopi$params,
                "n=5000, without misclassification")
cat("Converged:", fit_nopi$converged, "| Iterations:", fit_nopi$iterations, "\n")

# --- 1d. Log-likelihood monotonicity check -----------------------------------
cat("\n--- 1d. Log-likelihood monotonicity (n=500) ---\n")
ll_history <- fit_500$history$loglik
ll_diffs   <- diff(ll_history[!is.na(ll_history)])
n_decreases <- sum(ll_diffs < -1e-6 * max(abs(ll_history), na.rm = TRUE))
cat("Number of LL decreases (beyond tolerance):", n_decreases, "\n")
cat("Min LL change:", round(min(ll_diffs), 6), "\n")
cat("All decreases within tolerance:",
    all(ll_diffs >= -1e-6 * max(abs(ll_history), na.rm = TRUE)), "\n")

cat("\n[§1 SUMMARY] If abs_err < 0.05 for n=5000 (for alpha, theta1, pi, sigma2_g),\n")
cat("  the code is correct. theta0 has higher relative error due to the rare-event\n")
cat("  rate (3%) and coarse discrete timegap categories — this is expected.\n")
cat("  sigma2_d does NOT appear: discrete_timegap = TRUE uses interval-censored\n")
cat("  Exp(lambda_d) for nonemployment durations; lambda_d is shown instead.\n")


# ==============================================================================
# §2  REAL DATA DIAGNOSTICS
# ==============================================================================
# Purpose: quantify the data-model mismatches in the post-ingest QLFS data.
# Requires: QLFS data ingested via scripts/ingest_data_3waves_SA.R
#
# Key distinction:
#   - Post-ingest data: from ingest script (tenure/timegap in years, zeros)
#   - Post-filter data: after estimate_pipeline.R ad-hoc filter
# ==============================================================================
library(tidyverse)
cat("\n\n", strrep("#", 70), "\n",
    "§2  REAL DATA DIAGNOSTICS\n",
    strrep("#", 70), "\n", sep = "")

if (!exists("df_qlfs")) {
  message("df_qlfs not found — sourcing ingest script ...")
  source("scripts/ingest_data_3waves_SA.R")
}

df <- as.data.table(df_qlfs)

cat("\n--- 2a. Sample size and employment rates ---\n")
cat("Total observations (post-ingest):", nrow(df), "\n")
cat("Employment rate wave 1:", round(mean(df$y1, na.rm = TRUE), 4), "\n")
cat("Employment rate wave 2:", round(mean(df$y2, na.rm = TRUE), 4), "\n")
cat("Employment rate wave 3:", round(mean(df$y3, na.rm = TRUE), 4), "\n")

# --- 2b. Raw observed transition matrix -------------------------------------
cat("\n--- 2b. Raw observed quarterly transition rates ---\n")
cat("(These are the 'naive' rates ignoring misclassification)\n\n")

trans_12 <- df[, .(
  n           = .N,
  pct         = round(.N / nrow(df) * 100, 2)
), by = .(y1, y2)][order(y1, y2)]

trans_23 <- df[, .(
  n           = .N,
  pct         = round(.N / nrow(df) * 100, 2)
), by = .(y2, y3)][order(y2, y3)]

cat("Wave 1 → 2 transitions:\n")
print(trans_12)
cat("\nWave 2 → 3 transitions:\n")
print(trans_23)

# Implied transition rates
df_emp1  <- df[y1 == 1]
df_nonemp1 <- df[y1 == 0]
obs_theta1 <- mean(df_emp1$y2,    na.rm = TRUE)  # P(y2=1 | y1=1)
obs_theta0 <- mean(df_nonemp1$y2, na.rm = TRUE)  # P(y2=1 | y1=0)
cat(sprintf(
  "\nObserved P(employed t+1 | employed t):    %.4f (%.1f%%)\n",
  obs_theta1, obs_theta1 * 100
))
cat(sprintf(
  "Observed P(employed t+1 | nonemployed t): %.4f (%.1f%%)\n",
  obs_theta0, obs_theta0 * 100
))
cat(sprintf(
  "Observed exit rate (1 - θ₁):              %.4f (%.1f%%)\n",
  1 - obs_theta1, (1 - obs_theta1) * 100
))
cat(sprintf(
  "Observed entry rate (θ₀):                 %.4f (%.1f%%)\n",
  obs_theta0, obs_theta0 * 100
))

# --- 2c. Never-worked prevalence --------------------------------------------
cat("\n--- 2c. Never-worked prevalence (post-ingest data) ---\n")
cat("NOTE: neverworked is NA for employed individuals (not asked in survey)\n\n")

nw_tab <- df[, .(
  count = .N,
  pct   = round(.N / nrow(df) * 100, 2)
), by = .(y1, neverworked1)][order(y1, neverworked1)]
print(nw_tab)

nw_nonemployed <- df[y1 == 0, .(
  total         = .N,
  never_worked  = sum(as.numeric(neverworked1) == 1, na.rm = TRUE),
  pct_nw        = round(mean(as.numeric(neverworked1) == 1, na.rm = TRUE) * 100, 2)
)]
cat("\nAmong wave-1 nonemployed:\n")
print(nw_nonemployed)

# --- 2d. Duration distributions (post-ingest) --------------------------------
cat("\n--- 2d. Duration distributions (post-ingest, employed subsample) ---\n")
cat("Post-ingest tenure for y1=1 (years):\n")
print(summary(df[y1 == 1, tenure1]))

cat("\nPost-ingest timegap for y1=0 (years):\n")
print(summary(df[y1 == 0, timegap1]))

cat("\nTimegap for never-worked vs ever-worked nonemployed:\n")
df_nonemp <- df[y1 == 0]
nw_dur <- df_nonemp[, .(
  n        = .N,
  mean_tg  = round(mean(timegap1, na.rm = TRUE), 2),
  median_tg = round(median(timegap1, na.rm = TRUE), 2),
  max_tg   = round(max(timegap1, na.rm = TRUE), 2)
), by = neverworked1][order(neverworked1)]
print(nw_dur)

# --- 2e. Zero-duration counts -----------------------------------------------
cat("\n--- 2e. Zero-duration counts (post-ingest data) ---\n")
cat("These zeros cause log_emg(0) = -Inf for misclassifying latent histories\n\n")

for (t in 1:3) {
  g_col <- paste0("tenure",  t)
  d_col <- paste0("timegap", t)
  y_col <- paste0("y", t)

  # How the EM model actually sees the data:
  # tenure = 0 for nonemployed  →  breaks EMG for h_t=1 latent states
  # timegap = 0 for employed    →  breaks EMG for h_t=0 latent states
  n_zero_g_nonemp <- df[get(y_col) == 0, sum(get(g_col) == 0, na.rm = TRUE)]
  n_zero_d_emp    <- df[get(y_col) == 1, sum(get(d_col) == 0, na.rm = TRUE)]
  n_zero_g_emp    <- df[get(y_col) == 1, sum(get(g_col) == 0, na.rm = TRUE)]
  n_zero_d_nonemp <- df[get(y_col) == 0, sum(get(d_col) == 0, na.rm = TRUE)]

  total_nonemp <- df[get(y_col) == 0, .N]
  total_emp    <- df[get(y_col) == 1, .N]

  cat(sprintf(
    "Wave %d: tenure=0 & y=0:  %6d / %6d (%.1f%%) — breaks misclassified-as-employed EMG\n",
    t, n_zero_g_nonemp, total_nonemp, 100 * n_zero_g_nonemp / total_nonemp
  ))
  cat(sprintf(
    "Wave %d: timegap=0 & y=1: %6d / %6d (%.1f%%) — breaks misclassified-as-nonemployed EMG\n",
    t, n_zero_d_emp, total_emp, 100 * n_zero_d_emp / total_emp
  ))
  cat(sprintf(
    "Wave %d: tenure>0 & y=1:  %6d / %6d (%.1f%%) — valid for matched employed\n",
    t, total_emp - n_zero_g_emp, total_emp,
    100 * (total_emp - n_zero_g_emp) / total_emp
  ))
}

# --- 2f. Tenure increment distribution (employed continuations) --------------
cat("\n--- 2f. Tenure increment distribution (employed continuations) ---\n")
cat("Model expects: (tenure_t - tenure_{t-1} - 0.25) ~ N(0, 2*sigma2_g)\n\n")

# Waves 2 and 3: both observed as employed, previous also employed
inc_g2 <- df[y1 == 1 & y2 == 1 & tenure1 > 0 & tenure2 > 0,
             tenure2 - tenure1 - 0.25]
inc_g3 <- df[y2 == 1 & y3 == 1 & tenure2 > 0 & tenure3 > 0,
             tenure3 - tenure2 - 0.25]
inc_g  <- c(inc_g2, inc_g3)

cat("Summary of (tenure_t - tenure_{t-1} - 0.25) for emp-emp pairs:\n")
print(summary(inc_g))
cat(sprintf("N observations: %d\n", length(inc_g)))
cat(sprintf("Mean (expect 0): %.4f\n", mean(inc_g, na.rm = TRUE)))
cat(sprintf("SD (= sqrt(2)*sigma_g): %.4f\n", sd(inc_g, na.rm = TRUE)))
cat(sprintf("Implied sigma_g: %.4f years\n",
            sd(inc_g, na.rm = TRUE) / sqrt(2)))

# --- 2g. Timegap increment distribution (nonemployed continuations) ----------
cat("\n--- 2g. Timegap increment distribution (nonemployed continuations) ---\n")
cat("NOTE: In discrete_timegap mode (default), the model uses interval-censored\n")
cat("  Exp(lambda_d) — NOT a Gaussian increment model with sigma2_d.\n")
cat("  The distribution below will show discrete lumps (category transitions),\n")
cat("  which is precisely why the continuous EMG / sigma2_d model was replaced.\n\n")

inc_d2 <- df[y1 == 0 & y2 == 0 & timegap1 > 0 & timegap2 > 0,
             timegap2 - timegap1 - 0.25]
inc_d3 <- df[y2 == 0 & y3 == 0 & timegap2 > 0 & timegap3 > 0,
             timegap3 - timegap2 - 0.25]
inc_d  <- c(inc_d2, inc_d3)

cat("Summary of (timegap_t - timegap_{t-1} - 0.25) for nonemp-nonemp pairs:\n")
print(summary(inc_d))
cat(sprintf("N observations: %d\n", length(inc_d)))
cat(sprintf("Mean (expect 0): %.4f\n", mean(inc_d, na.rm = TRUE)))
cat(sprintf("SD: %.4f years\n", sd(inc_d, na.rm = TRUE)))
cat(sprintf("Implied sigma_d: %.4f years\n",
            sd(inc_d, na.rm = TRUE) / sqrt(2)))

# Unique values reveal the discrete grid
cat("\nUnique timegap increment values (shows discreteness):\n")
cat(round(sort(unique(inc_d))[1:min(20, length(unique(inc_d)))], 4), "\n")

# --- 2h. Plots ---------------------------------------------------------------
cat("\n--- 2h. Saving diagnostic plots to troubleshoot/ ---\n")

# Plot 1: Timegap distribution for nonemployed, by never-worked
p_timegap <- ggplot(
  df[y1 == 0 & !is.na(neverworked1)],
  aes(x = timegap1, fill = factor(neverworked1, labels = c("Ever-worked", "Never-worked")))
) +
  geom_histogram(bins = 40, position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c("#2166ac", "#d6604d"), name = "Type") +
  labs(
    title    = "Timegap distribution for nonemployed (wave 1)",
    subtitle = "Post-ingest data. Never-worked timegap = (age - educ - 6) / 12",
    x        = "Timegap (years)",
    y        = "Count"
  ) +
  theme_minimal()
ggsave("troubleshoot/plot_timegap_dist.png", p_timegap, width = 10, height = 6, dpi = 150)
cat("Saved: troubleshoot/plot_timegap_dist.png\n")

# Plot 2: Tenure increment distribution vs Normal expectation
inc_dt <- data.table(increment = inc_g[!is.na(inc_g)])
sigma_g_implied <- sd(inc_g, na.rm = TRUE) / sqrt(2)
p_tenure_inc <- ggplot(inc_dt[abs(increment) < 2], aes(x = increment)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80,
                 fill = "#2166ac", alpha = 0.7) +
  stat_function(
    fun  = function(x) dnorm(x, 0, sqrt(2) * sigma_g_implied),
    colour = "#d6604d", linewidth = 1
  ) +
  labs(
    title    = "Tenure increment distribution (employed continuations)",
    subtitle = "Blue = data; red = fitted N(0, 2*sigma_g^2)",
    x        = "Increment - 0.25 (years)",
    y        = "Density"
  ) +
  theme_minimal()
ggsave("troubleshoot/plot_tenure_increment.png", p_tenure_inc,
       width = 10, height = 6, dpi = 150)
cat("Saved: troubleshoot/plot_tenure_increment.png\n")

# Plot 3: Timegap increment — discrete structure
inc_d_dt <- data.table(increment = inc_d[!is.na(inc_d)])
p_timegap_inc <- ggplot(inc_d_dt[abs(increment) < 5], aes(x = increment)) +
  geom_histogram(bins = 120, fill = "#2166ac", alpha = 0.7) +
  labs(
    title    = "Timegap increment distribution (nonemployed continuations)",
    subtitle = "Discrete lumps reveal categorical midpoint structure",
    x        = "Increment - 0.25 (years)",
    y        = "Count"
  ) +
  theme_minimal()
ggsave("troubleshoot/plot_timegap_increment.png", p_timegap_inc,
       width = 10, height = 6, dpi = 150)
cat("Saved: troubleshoot/plot_timegap_increment.png\n")


# ==============================================================================
# §3  PROGRESSIVE STRESS TESTS
# ==============================================================================
# Purpose: starting from clean synthetic data, introduce real-data pathologies
# one at a time to isolate which feature drives the implausible estimates.
#
# All tests use n=2000 for speed. True parameters match §1.
# ==============================================================================

cat("\n\n", strrep("#", 70), "\n",
    "§3  PROGRESSIVE STRESS TESTS\n",
    strrep("#", 70), "\n", sep = "")

set.seed(999)
df_base <- simulate_panel(
  n        = 2000,
  alpha    = true_params$alpha,
  theta1   = true_params$theta1,
  theta0   = true_params$theta0,
  pi       = true_params$pi,
  sigma2_g = true_params$sigma2_g,
  # sigma2_d omitted: discrete_timegap = TRUE (default) sets sigma_d = NA
  seed     = 999
)

# Fit baseline: clean synthetic data
cat("\n--- Test 0 (Baseline): Clean synthetic data, n=2000 ---\n")
fit_base <- em_fit_tenure(df_base, misclassification = TRUE, verbose = 1L)
.print_recovery(true_params, fit_base$params, "Baseline (clean synthetic)")

# --- Test A: Zero inappropriate durations ------------------------------------
cat("\n--- Test A: Zero inappropriate durations ---\n")
cat("Mimics ingest encoding: tenure=0 where y=0, timegap=0 where y=1\n")

df_A <- copy(as.data.table(df_base))
for (t in 1:3) {
  y_col <- paste0("y",      t)
  g_col <- paste0("tenure",  t)
  d_col <- paste0("timegap", t)
  df_A[get(y_col) == 0, (g_col) := 0]
  df_A[get(y_col) == 1, (d_col) := 0]
}
df_A <- as.data.frame(df_A)

# Apply the same ad-hoc filter as estimate_pipeline.R
for (t in 1:3) {
  y <- df_A[[paste0("y",      t)]]
  g <- df_A[[paste0("tenure",  t)]]
  d <- df_A[[paste0("timegap", t)]]
  keep <- !((y == 1 & g <= 0) | (y == 0 & d <= 0))
  df_A  <- df_A[keep, ]
}
cat(sprintf("Rows after ad-hoc filter: %d (from %d)\n", nrow(df_A), nrow(df_base)))

fit_A <- em_fit_tenure(df_A, misclassification = TRUE, verbose = 1L)
.print_recovery(true_params, fit_A$params, "Test A: zero inappropriate durations")

# --- Test B: Categorical timegaps -------------------------------------------
cat("\n--- Test B: Categorical timegaps ---\n")
cat("Buckets timegap into QLFS categories and replaces with midpoints\n")

df_B <- copy(as.data.table(df_base))
for (t in 1:3) {
  d_col <- paste0("timegap", t)
  df_B[, (d_col) := .bucket_timegap(get(d_col))]
}
df_B <- as.data.frame(df_B)

fit_B <- em_fit_tenure(df_B, misclassification = TRUE, verbose = 1L)
.print_recovery(true_params, fit_B$params, "Test B: categorical timegaps")

# --- Test C: Never-worked injection -----------------------------------------
cat("\n--- Test C: Never-worked injection ---\n")
cat("Replaces ~40% of nonemployed timegaps with (10-20) years\n")

df_C <- copy(as.data.table(df_base))
set.seed(777)
for (t in 1:3) {
  y_col <- paste0("y",      t)
  d_col <- paste0("timegap", t)
  rows_nonemp <- which(df_C[[y_col]] == 0)
  n_nw        <- round(0.40 * length(rows_nonemp))
  nw_idx      <- sample(rows_nonemp, n_nw)
  # Imputed duration = uniform(10, 20) years, mimicking age - educ - 6
  df_C[nw_idx, (d_col) := runif(.N, 10, 20)]
}
df_C <- as.data.frame(df_C)

fit_C <- em_fit_tenure(df_C, misclassification = TRUE, verbose = 1L)
.print_recovery(true_params, fit_C$params, "Test C: never-worked injection")

# --- Test D: All pathologies combined ----------------------------------------
cat("\n--- Test D: All pathologies combined (A + B + C) ---\n")

df_D <- copy(as.data.table(df_base))
# Apply C: never-worked injection
set.seed(888)
for (t in 1:3) {
  y_col <- paste0("y", t)
  d_col <- paste0("timegap", t)
  rows_nonemp <- which(df_D[[y_col]] == 0)
  n_nw        <- round(0.40 * length(rows_nonemp))
  nw_idx      <- sample(rows_nonemp, n_nw)
  df_D[nw_idx, (d_col) := runif(.N, 10, 20)]
}
# Apply B: categorical timegaps (after C, so never-worked stay at raw values)
for (t in 1:3) {
  d_col <- paste0("timegap", t)
  y_col <- paste0("y", t)
  # Only bucket the non-never-worked nonemployed
  df_D[get(y_col) == 0 & get(d_col) < 10, (d_col) := .bucket_timegap(get(d_col))]
}
# Apply A: zero inappropriate durations
for (t in 1:3) {
  y_col <- paste0("y", t)
  g_col <- paste0("tenure",  t)
  d_col <- paste0("timegap", t)
  df_D[get(y_col) == 0, (g_col) := 0]
  df_D[get(y_col) == 1, (d_col) := 0]
}
df_D <- as.data.frame(df_D)
# Apply ad-hoc filter
for (t in 1:3) {
  y <- df_D[[paste0("y",      t)]]
  g <- df_D[[paste0("tenure",  t)]]
  d <- df_D[[paste0("timegap", t)]]
  keep <- !((y == 1 & g <= 0) | (y == 0 & d <= 0))
  df_D  <- df_D[keep, ]
}
cat(sprintf("Rows after ad-hoc filter: %d\n", nrow(df_D)))

fit_D <- em_fit_tenure(df_D, misclassification = TRUE, verbose = 1L)
.print_recovery(true_params, fit_D$params, "Test D: all pathologies combined")

# --- Test E: Remedy — duration floor + exclude never-worked ------------------
cat("\n--- Test E: Remedy (duration floor on Test D data) ---\n")
cat("Apply floor 0.125 years to all zero durations in df_D\n")
cat("(Cannot 'undo' never-worked injection in synthetic data, but floor helps)\n")

DURATION_FLOOR <- 0.125

df_E <- as.data.table(df_D)
for (t in 1:3) {
  g_col <- paste0("tenure",  t)
  d_col <- paste0("timegap", t)
  df_E[get(g_col) <= 0, (g_col) := DURATION_FLOOR]
  df_E[get(d_col) <= 0, (d_col) := DURATION_FLOOR]
}
df_E <- as.data.frame(df_E)

fit_E <- em_fit_tenure(df_E, misclassification = TRUE, verbose = 1L)
.print_recovery(true_params, fit_E$params,
                "Test E: remedy (duration floor 0.125yr)")

# Stress test summary table
cat("\n", strrep("=", 70), "\n")
cat("§3 SUMMARY: parameter estimates across stress tests\n")
cat(strrep("=", 70), "\n\n")

stress_summary <- rbindlist(list(
  data.table(test = "True params",  theta1 = true_params$theta1, theta0 = true_params$theta0, pi = true_params$pi),
  data.table(test = "0. Baseline",  theta1 = fit_base$params$theta1, theta0 = fit_base$params$theta0, pi = fit_base$params$pi),
  data.table(test = "A. Zero dur.", theta1 = fit_A$params$theta1,    theta0 = fit_A$params$theta0,    pi = fit_A$params$pi),
  data.table(test = "B. Cat. tgap", theta1 = fit_B$params$theta1,    theta0 = fit_B$params$theta0,    pi = fit_B$params$pi),
  data.table(test = "C. Never-wkd", theta1 = fit_C$params$theta1,    theta0 = fit_C$params$theta0,    pi = fit_C$params$pi),
  data.table(test = "D. All path.", theta1 = fit_D$params$theta1,    theta0 = fit_D$params$theta0,    pi = fit_D$params$pi),
  data.table(test = "E. Remedy",    theta1 = fit_E$params$theta1,    theta0 = fit_E$params$theta0,    pi = fit_E$params$pi)
))
stress_summary[, `:=`(
  theta1 = round(theta1, 4),
  theta0 = round(theta0, 4),
  pi     = round(pi,     4)
)]
print(stress_summary)


# ==============================================================================
# §4  TARGETED RE-ESTIMATION ON REAL QLFS DATA
# ==============================================================================
# Purpose: apply remedies R1 and R2 (from DIAGNOSIS.md) to the real QLFS data
# and check whether estimates become consistent with the simpler model baseline
# (theta0 ~ 3%, theta1 ~ 97%, pi ~ 3%).
#
# Requires: df_qlfs from ingest script
# ==============================================================================

cat("\n\n", strrep("#", 70), "\n",
    "§4  TARGETED RE-ESTIMATION ON REAL QLFS DATA\n",
    strrep("#", 70), "\n", sep = "")

if (!exists("df_qlfs")) {
  message("df_qlfs not found — sourcing ingest script ...")
  source("scripts/ingest_data_3waves_SA.R")
}

# Reference: original estimate_pipeline.R results
cat("\n--- Simpler model baseline (from estimate_models_3waves_SA.R) ---\n")
cat("theta0 (entry) ≈ 0.03,  theta1 (stay-employed) ≈ 0.97,  pi ≈ 0.03\n")
cat("These come from binary-sequence MLE, no duration data\n")

cat("\n--- Reference: EM-tenure on full post-filter data (original result) ---\n")
cat("theta1 ≈ 0.465,  theta0 ≈ 0.418,  pi ≈ 0.266 — known implausible values\n")

DURATION_FLOOR <- 0.125  # 1.5 months

# Helper: apply the estimate_pipeline.R ad-hoc filter
.apply_adhoc_filter <- function(df) {
  for (t in 1:3) {
    y <- df[[paste0("y",      t)]]
    g <- df[[paste0("tenure",  t)]]
    d <- df[[paste0("timegap", t)]]
    keep <- !((y == 1 & g <= 0) | (y == 0 & d <= 0))
    df  <- df[keep, ]
  }
  df
}

# Helper: apply duration floor to all zero/negative durations
.apply_duration_floor <- function(df, floor = DURATION_FLOOR) {
  for (t in 1:3) {
    g_col <- paste0("tenure",  t)
    d_col <- paste0("timegap", t)
    df[[g_col]] <- ifelse(df[[g_col]] <= 0, floor, df[[g_col]])
    df[[d_col]] <- ifelse(df[[d_col]] <= 0, floor, df[[d_col]])
  }
  df
}

# Shared initial params for fair comparison
# sigma2_d is excluded: discrete_timegap = TRUE (default) does not use it.
# lambda_d is initialised directly from theta0 via the CTMC link.
custom_init <- list(
  alpha    = 0.47,
  theta1   = 0.95,
  theta0   = 0.05,
  pi       = 0.03,
  sigma2_g = 0.5,
  lambda_g = ctmc_lambda_from_theta(0.95),
  lambda_d = ctmc_lambda_from_theta(0.05)
)

# --- Subgroup 1: Exclude never-worked ----------------------------------------
cat("\n--- Subgroup 1: Exclude never-worked (neverworked1/2/3 == 1) ---\n")
df_sub1 <- df_qlfs |>
  dplyr::filter(is.na(neverworked1) | neverworked1 == 0) |>
  dplyr::filter(is.na(neverworked2) | neverworked2 == 0) |>
  dplyr::filter(is.na(neverworked3) | neverworked3 == 0)
df_sub1 <- .apply_adhoc_filter(df_sub1)
cat(sprintf("n = %d (from %d; excluded %d never-worked)\n",
            nrow(df_sub1), nrow(df_qlfs), nrow(df_qlfs) - nrow(df_sub1)))
cat("Employment rate wave 1:", round(mean(df_sub1$y1, na.rm = TRUE), 4), "\n")

fit_sub1 <- em_fit_tenure(
  df_sub1,
  params0           = custom_init,
  misclassification = TRUE,
  verbose           = 2L
)
cat(sprintf(
  "\nSubgroup 1 estimates: theta1=%.4f  theta0=%.4f  pi=%.4f  alpha=%.4f\n",
  fit_sub1$params$theta1, fit_sub1$params$theta0,
  fit_sub1$params$pi,     fit_sub1$params$alpha
))

# --- Subgroup 2: Exclude never-worked + apply duration floor -----------------
cat("\n--- Subgroup 2: Exclude never-worked + duration floor (0.125 yr) ---\n")
df_sub2 <- df_sub1 |>
  (\(d) .apply_duration_floor(d, DURATION_FLOOR))()
cat(sprintf("n = %d\n", nrow(df_sub2)))

fit_sub2 <- em_fit_tenure(
  df_sub2,
  params0           = custom_init,
  misclassification = TRUE,
  verbose           = 2L
)
cat(sprintf(
  "\nSubgroup 2 estimates: theta1=%.4f  theta0=%.4f  pi=%.4f  alpha=%.4f\n",
  fit_sub2$params$theta1, fit_sub2$params$theta0,
  fit_sub2$params$pi,     fit_sub2$params$alpha
))

# --- Subgroup 3: Exclude never-worked + floor + short timegap only -----------
cat("\n--- Subgroup 3: Subgroup 2 + restrict to timegap < 1 year ---\n")
cat("(Tests whether the categorical midpoint problem is driving remaining bias)\n")
df_sub3 <- df_sub2 |>
  dplyr::filter((timegap1 < 1 | y1 == 1)) |>
  dplyr::filter((timegap2 < 1 | y2 == 1)) |>
  dplyr::filter((timegap3 < 1 | y3 == 1))
cat(sprintf("n = %d (from %d after further restriction)\n",
            nrow(df_sub3), nrow(df_sub2)))
cat("Employment rate wave 1:", round(mean(df_sub3$y1, na.rm = TRUE), 4), "\n")

fit_sub3 <- em_fit_tenure(
  df_sub3,
  params0           = custom_init,
  misclassification = TRUE,
  verbose           = 2L
)
cat(sprintf(
  "\nSubgroup 3 estimates: theta1=%.4f  theta0=%.4f  pi=%.4f  alpha=%.4f\n",
  fit_sub3$params$theta1, fit_sub3$params$theta0,
  fit_sub3$params$pi,     fit_sub3$params$alpha
))

# Final comparison table
cat("\n", strrep("=", 70), "\n")
cat("§4 FINAL COMPARISON TABLE\n")
cat(strrep("=", 70), "\n\n")

comparison <- rbindlist(list(
  data.table(model = "Baseline (simpler MLE)", theta1 = 0.97, theta0 = 0.03, pi = 0.03, n = NA_integer_),
  data.table(model = "EM original (full data)", theta1 = 0.465, theta0 = 0.418, pi = 0.266, n = nrow(df_qlfs)),
  data.table(model = "Sub1: no never-worked",   theta1 = round(fit_sub1$params$theta1, 4), theta0 = round(fit_sub1$params$theta0, 4), pi = round(fit_sub1$params$pi, 4), n = nrow(df_sub1)),
  data.table(model = "Sub2: + duration floor",  theta1 = round(fit_sub2$params$theta1, 4), theta0 = round(fit_sub2$params$theta0, 4), pi = round(fit_sub2$params$pi, 4), n = nrow(df_sub2)),
  data.table(model = "Sub3: + short timegap",   theta1 = round(fit_sub3$params$theta1, 4), theta0 = round(fit_sub3$params$theta0, 4), pi = round(fit_sub3$params$pi, 4), n = nrow(df_sub3))
))
print(comparison)

cat("\n[INTERPRETATION]\n")
cat("If Sub2 or Sub3 gives theta1 ≈ 0.97, theta0 ≈ 0.03, pi ≈ 0.03:\n")
cat("  → Issues 1 (zero durations) and 2 (never-worked) are confirmed as root causes.\n")
cat("  → Fix the ingest script to apply the duration floor by default.\n")
cat("  → Decide whether to exclude never-worked from the EM sample.\n\n")
cat("If estimates remain distorted even in Sub3:\n")
cat("  → The categorical timegap structure (Issue 3) is also a major driver.\n")
cat("  → Model redesign is needed (discrete/censored duration model).\n\n")
cat("See troubleshoot/DIAGNOSIS.md for the proposed discrete/censored model.\n")


# ==============================================================================
# §5  DISCRETE/CENSORED TIMEGAP MODEL: PRE-IMPLEMENTATION DIAGNOSTICS
# ==============================================================================
# Purpose: validate the proposed discrete interval-censored Exp(λ_d) model
# BEFORE modifying any estimation code. This section verifies:
#   5a. Category distribution of timegap in QLFS data (excluding never-worked)
#   5b. Observed category-to-category transitions vs theoretical predictions
#   5c. Interval probabilities vs EMG midpoint densities (distortion analysis)
#   5d. Within-panel start category check (expect ~100% in category 1)
#   5e. Nearest-non-zero imputation: distribution preview
#
# Requires: df_qlfs from ingest script
# ==============================================================================

cat("\n\n", strrep("#", 70), "\n",
    "§5  DISCRETE/CENSORED TIMEGAP MODEL: PRE-IMPLEMENTATION DIAGNOSTICS\n",
    strrep("#", 70), "\n", sep = "")

if (!exists("df_qlfs")) {
  message("df_qlfs not found — sourcing ingest script ...")
  source("scripts/ingest_data_3waves_SA.R")
}

df5 <- as.data.table(df_qlfs)

# Exclude never-worked (Issue 2 remedy — standard for all downstream analysis)
df5 <- df5[is.na(neverworked1) | neverworked1 != 1]
df5 <- df5[is.na(neverworked2) | neverworked2 != 1]
df5 <- df5[is.na(neverworked3) | neverworked3 != 1]
cat(sprintf("\nSample after excluding never-worked: n = %d\n", nrow(df5)))

# --- Category boundaries and mapping -----------------------------------------
# These are the QLFS timegap category intervals (in years)
.BOUNDS_YEARS <- c(0, 0.25, 0.5, 0.75, 1.0, 3.0, 5.0, Inf)
.N_CATS <- 7L

# The ingest script maps codes 1-7 to midpoints (in months) then /12.
# Recover the category code from the midpoint value.
.MIDPOINTS_YEARS <- c(1.5, 4.5, 7.5, 10.5, 24, 48, 90) / 12

.midpoint_to_cat <- function(midpoint_years) {
  # Map midpoints back to category codes 1-7
  # Midpoints: 0.125, 0.375, 0.625, 0.875, 2.0, 4.0, 7.5 (years)
  fcase(
    abs(midpoint_years - .MIDPOINTS_YEARS[1]) < 0.001, 1L,
    abs(midpoint_years - .MIDPOINTS_YEARS[2]) < 0.001, 2L,
    abs(midpoint_years - .MIDPOINTS_YEARS[3]) < 0.001, 3L,
    abs(midpoint_years - .MIDPOINTS_YEARS[4]) < 0.001, 4L,
    abs(midpoint_years - .MIDPOINTS_YEARS[5]) < 0.01,  5L,
    abs(midpoint_years - .MIDPOINTS_YEARS[6]) < 0.01,  6L,
    abs(midpoint_years - .MIDPOINTS_YEARS[7]) < 0.1,   7L,
    default = NA_integer_
  )
}

# Add category codes to the data
for (t in 1:3) {
  d_col   <- paste0("timegap", t)
  cat_col <- paste0("tg_cat",  t)
  y_col   <- paste0("y", t)
  # Only meaningful for nonemployed (y = 0); employed have timegap = 0
  df5[get(y_col) == 0, (cat_col) := .midpoint_to_cat(get(d_col))]
  df5[get(y_col) == 1, (cat_col) := NA_integer_]
}

# --- 5a. Category distribution ------------------------------------------------
cat("\n--- 5a. Timegap category distribution (nonemployed, excluding never-worked) ---\n\n")

cat_dist <- rbindlist(lapply(1:3, function(t) {
  cat_col <- paste0("tg_cat", t)
  y_col   <- paste0("y", t)
  df5[get(y_col) == 0 & !is.na(get(cat_col)), .(
    wave = t,
    n    = .N,
    pct  = round(.N / df5[get(y_col) == 0, .N] * 100, 1)
  ), by = .(cat = get(cat_col))]
}))[order(wave, cat)]

cat("Category distribution by wave:\n")
print(dcast(cat_dist, cat ~ paste0("wave", wave), value.var = "pct", fill = 0))

# Marginal distribution across all waves
cat_marginal <- rbindlist(lapply(1:3, function(t) {
  cat_col <- paste0("tg_cat", t)
  y_col   <- paste0("y", t)
  df5[get(y_col) == 0 & !is.na(get(cat_col)), .(cat = get(cat_col))]
}))
cat_marginal_tab <- cat_marginal[, .(n = .N, pct = round(.N / .N * 100, 1)), by = cat][order(cat)]
cat_marginal_tab[, pct := round(n / sum(n) * 100, 1)]
cat("\nMarginal category distribution (all waves pooled):\n")
print(cat_marginal_tab)

# --- 5b. Observed category-to-category transitions ---------------------------
cat("\n--- 5b. Observed category transitions (nonemployed continuations) ---\n")
cat("These are pairs (c_{t-1}, c_t) where y_{t-1}=0 and y_t=0\n\n")

trans_obs <- rbindlist(lapply(2:3, function(t) {
  cat_prev <- paste0("tg_cat", t - 1)
  cat_curr <- paste0("tg_cat", t)
  y_prev   <- paste0("y",      t - 1)
  y_curr   <- paste0("y",      t)
  df5[get(y_prev) == 0 & get(y_curr) == 0 &
      !is.na(get(cat_prev)) & !is.na(get(cat_curr)),
      .(cat_prev = get(cat_prev), cat_curr = get(cat_curr))]
}))

cat(sprintf("Total nonemployed continuation pairs: %d\n", nrow(trans_obs)))

# Observed transition matrix (counts)
trans_mat_obs <- dcast(
  trans_obs[, .N, by = .(cat_prev, cat_curr)],
  cat_prev ~ cat_curr,
  value.var = "N", fill = 0
)
cat("\nObserved transition counts (rows = c_{t-1}, cols = c_t):\n")
print(trans_mat_obs)

# Observed transition probabilities (row-normalised)
trans_mat_pct <- copy(trans_mat_obs)
row_sums <- rowSums(trans_mat_pct[, -1])
for (col in names(trans_mat_pct)[-1]) {
  set(trans_mat_pct, j = col, value = round(trans_mat_pct[[col]] / row_sums * 100, 1))
}
cat("\nObserved transition probabilities (row %):\n")
print(trans_mat_pct)

# Theoretical transition matrix from the sparse model
cat("\n--- Theoretical transition matrix (sparse Exp(λ_d) model) ---\n")
cat("For categories 1-4: deterministic (one destination).\n")
cat("For categories 5-6: probabilistic (two destinations).\n")
cat("For category 7: absorbing (stays at 7).\n\n")

# Using θ₀ = 0.03 (literature value) → λ_d = -log(0.97)/3
theta0_ref <- 0.03
lambda_d_ref <- -log(1 - theta0_ref) / 3

.interval_prob <- function(lambda, a, b) {
  if (is.infinite(b)) {
    return(exp(-lambda * a))
  }
  exp(-lambda * a) - exp(-lambda * b)
}

.transition_prob <- function(lambda, j, k, bounds = .BOUNDS_YEARS) {
  a_j <- bounds[j]
  b_j <- bounds[j + 1]
  a_k <- bounds[k]
  b_k <- bounds[k + 1]

  # Intersection of [a_j, b_j) and [a_k - 0.25, b_k - 0.25)
  L <- max(a_j, a_k - 0.25)
  U <- min(b_j, if (is.infinite(b_k)) Inf else b_k - 0.25)

  if (is.finite(L) && is.finite(U) && L >= U) return(0)
  if (is.infinite(U) && is.infinite(L)) return(0)

  numer <- .interval_prob(lambda, L, U)
  denom <- .interval_prob(lambda, a_j, b_j)

  if (denom < 1e-300) return(0)
  return(numer / denom)
}

# Build 7x7 theoretical transition matrix
theo_mat <- matrix(0, nrow = .N_CATS, ncol = .N_CATS)
for (j in 1:.N_CATS) {
  for (k in 1:.N_CATS) {
    theo_mat[j, k] <- .transition_prob(lambda_d_ref, j, k)
  }
}
theo_df <- data.table(
  from = 1:.N_CATS,
  as.data.table(round(theo_mat * 100, 1))
)
setnames(theo_df, c("from", paste0("to_", 1:.N_CATS)))
cat(sprintf("Theoretical transition matrix (row %%, λ_d = %.6f, θ₀ = %.4f):\n",
            lambda_d_ref, theta0_ref))
print(theo_df)

cat("\nKey predictions:\n")
cat("  - Categories 1→2, 2→3, 3→4, 4→5: should be ~100% (deterministic)\n")
cat("  - Category 5→5 vs 5→6: split depends on λ_d\n")
cat("  - Category 6→6 vs 6→7: split depends on λ_d\n")
cat("  - Category 7→7: should be ~100% (absorbing)\n")

# --- 5c. Interval probabilities vs EMG midpoint densities --------------------
cat("\n--- 5c. Interval probabilities vs EMG midpoint densities ---\n")
cat("Shows how midpoint-based EMG distorts the likelihood relative to the\n")
cat("correct interval probability.\n\n")

# Use reference λ_d from θ₀ = 0.03 and a moderate σ²_d = 0.01
sigma2_d_ref <- 0.01

cat(sprintf("Reference parameters: θ₀ = %.4f, λ_d = %.6f, σ²_d = %.4f\n\n",
            theta0_ref, lambda_d_ref, sigma2_d_ref))

comparison_5c <- data.table(
  cat       = 1:.N_CATS,
  a_years   = .BOUNDS_YEARS[1:.N_CATS],
  b_years   = .BOUNDS_YEARS[2:(.N_CATS + 1)],
  midpoint  = .MIDPOINTS_YEARS
)

comparison_5c[, interval_prob := mapply(function(a, b) {
  .interval_prob(lambda_d_ref, a, b)
}, a_years, b_years)]

comparison_5c[, emg_at_midpoint := exp(log_emg(midpoint, lambda_d_ref, sigma2_d_ref))]

# Normalise EMG densities to sum to 1 for comparison (multiply by bin width)
comparison_5c[, bin_width := fifelse(is.infinite(b_years), 5.0, b_years - a_years)]
comparison_5c[, emg_approx_prob := emg_at_midpoint * bin_width]
comparison_5c[, emg_approx_prob := emg_approx_prob / sum(emg_approx_prob)]

comparison_5c[, ratio := round(emg_approx_prob / interval_prob, 3)]
comparison_5c[, interval_prob := round(interval_prob, 6)]
comparison_5c[, emg_approx_prob := round(emg_approx_prob, 6)]

cat("Comparison: interval probability vs EMG midpoint approximation\n")
print(comparison_5c[, .(cat, a_years, b_years, midpoint,
                        interval_prob, emg_approx_prob, ratio)])
cat("\nRatio > 1 means EMG over-weights that category; < 1 means under-weights.\n")
cat("Large deviations (especially for wide categories 5-7) show the midpoint\n")
cat("approximation is poor — this is what the discrete model fixes.\n")

# --- 5d. Within-panel start category check -----------------------------------
cat("\n--- 5d. Within-panel start category check ---\n")
cat("For new nonemployment spells (y_{t-1}=1, y_t=0), the true duration\n")
cat("is at most 0.25 years → should always be category 1.\n\n")

starts <- rbindlist(lapply(2:3, function(t) {
  y_prev  <- paste0("y", t - 1)
  y_curr  <- paste0("y", t)
  cat_col <- paste0("tg_cat", t)
  df5[get(y_prev) == 1 & get(y_curr) == 0 & !is.na(get(cat_col)),
      .(wave = t, cat = get(cat_col))]
}))

if (nrow(starts) > 0) {
  cat(sprintf("Total within-panel nonemployment starts: %d\n", nrow(starts)))
  start_tab <- starts[, .(n = .N, pct = round(.N / nrow(starts) * 100, 1)), by = cat][order(cat)]
  print(start_tab)
  cat(sprintf("\nFraction in category 1: %.1f%%\n", start_tab[cat == 1, pct]))
  cat("If this is close to 100%, the deterministic start emission is well-justified.\n")
  cat("If substantially less, some starts report longer durations — investigate.\n")
} else {
  cat("No within-panel starts found (all individuals same state across waves?)\n")
}

# --- 5e. Nearest-non-zero imputation preview ----------------------------------
cat("\n--- 5e. Nearest-non-zero imputation preview ---\n")
cat("For employed individuals (y_t=1), timegap_t = 0 (set by ingest).\n")
cat("The proposed imputation uses the nearest non-zero timegap from another wave.\n")
cat("Preview: what values would be imputed?\n\n")

# For each employed wave, find the nearest non-zero timegap from another wave
impute_preview <- rbindlist(lapply(1:3, function(t) {
  y_col <- paste0("y", t)
  other_waves <- setdiff(1:3, t)

  # Individuals employed at wave t
  idx <- which(df5[[y_col]] == 1)
  if (length(idx) == 0) return(data.table())

  # For each, find nearest non-zero timegap from other waves
  nearest_val <- rep(NA_real_, length(idx))
  nearest_cat <- rep(NA_integer_, length(idx))
  for (t2 in other_waves) {
    y2_col <- paste0("y", t2)
    d2_col <- paste0("timegap", t2)
    cat2_col <- paste0("tg_cat", t2)
    # Prefer closer wave
    can_use <- (df5[[y2_col]][idx] == 0) & (df5[[d2_col]][idx] > 0) & is.na(nearest_val)
    nearest_val[can_use] <- df5[[d2_col]][idx[can_use]]
    nearest_cat[can_use] <- df5[[cat2_col]][idx[can_use]]
  }

  data.table(
    wave = t,
    n_employed = length(idx),
    n_has_donor = sum(!is.na(nearest_val)),
    n_needs_floor = sum(is.na(nearest_val)),
    pct_has_donor = round(sum(!is.na(nearest_val)) / length(idx) * 100, 1),
    pct_needs_floor = round(sum(is.na(nearest_val)) / length(idx) * 100, 1)
  )
}))

cat("Imputation source availability (employed observations):\n")
print(impute_preview)

cat(sprintf(
  "\nTotal employed-wave observations: %d\n",
  sum(impute_preview$n_employed)
))
cat(sprintf(
  "  With donor from another wave: %d (%.1f%%)\n",
  sum(impute_preview$n_has_donor),
  sum(impute_preview$n_has_donor) / sum(impute_preview$n_employed) * 100
))
cat(sprintf(
  "  Needs floor (all 3 waves employed): %d (%.1f%%)\n",
  sum(impute_preview$n_needs_floor),
  sum(impute_preview$n_needs_floor) / sum(impute_preview$n_employed) * 100
))

# Distribution of imputed categories (from donors)
imputed_cats <- rbindlist(lapply(1:3, function(t) {
  y_col <- paste0("y", t)
  other_waves <- setdiff(1:3, t)
  idx <- which(df5[[y_col]] == 1)
  if (length(idx) == 0) return(data.table())

  nearest_cat <- rep(NA_integer_, length(idx))
  for (t2 in other_waves) {
    y2_col <- paste0("y", t2)
    cat2_col <- paste0("tg_cat", t2)
    can_use <- (df5[[y2_col]][idx] == 0) & !is.na(df5[[cat2_col]][idx]) & is.na(nearest_cat)
    nearest_cat[can_use] <- df5[[cat2_col]][idx[can_use]]
  }
  data.table(imputed_cat = nearest_cat[!is.na(nearest_cat)])
}))

if (nrow(imputed_cats) > 0) {
  cat("\nDistribution of imputed timegap categories (from nearest non-zero donor):\n")
  imp_tab <- imputed_cats[, .(n = .N, pct = round(.N / nrow(imputed_cats) * 100, 1)),
                           by = imputed_cat][order(imputed_cat)]
  print(imp_tab)
}

# For tenure: symmetric check (nonemployed with imputed tenure from nearest wave)
tenure_preview <- rbindlist(lapply(1:3, function(t) {
  y_col <- paste0("y", t)
  other_waves <- setdiff(1:3, t)
  idx <- which(df5[[y_col]] == 0)
  if (length(idx) == 0) return(data.table())

  nearest_val <- rep(NA_real_, length(idx))
  for (t2 in other_waves) {
    y2_col <- paste0("y", t2)
    g2_col <- paste0("tenure", t2)
    can_use <- (df5[[y2_col]][idx] == 1) & (df5[[g2_col]][idx] > 0) & is.na(nearest_val)
    nearest_val[can_use] <- df5[[g2_col]][idx[can_use]]
  }

  data.table(
    wave = t,
    n_nonemployed = length(idx),
    n_has_donor = sum(!is.na(nearest_val)),
    n_needs_floor = sum(is.na(nearest_val))
  )
}))

cat("\nTenure imputation preview (nonemployed observations):\n")
print(tenure_preview)
cat(sprintf(
  "  Needs floor (all 3 waves nonemployed): %d (%.1f%% of nonemployed-wave obs)\n",
  sum(tenure_preview$n_needs_floor),
  sum(tenure_preview$n_needs_floor) / sum(tenure_preview$n_nonemployed) * 100
))

# --- 5f. Diagnostic plot: interval probs vs EMG at multiple λ_d --------------
cat("\n--- 5f. Saving §5 diagnostic plots to troubleshoot/ ---\n")

# Plot: interval probability vs category for multiple λ_d values
lambda_vals <- c(0.005, 0.01, 0.02, 0.05, 0.1)
plot_data_5f <- rbindlist(lapply(lambda_vals, function(lam) {
  theta0_val <- 1 - exp(-lam * 3)
  data.table(
    cat = 1:.N_CATS,
    midpoint = .MIDPOINTS_YEARS,
    prob = sapply(1:.N_CATS, function(k) .interval_prob(lam, .BOUNDS_YEARS[k], .BOUNDS_YEARS[k + 1])),
    lambda_d = lam,
    theta0 = theta0_val,
    label = sprintf("λ=%.3f (θ₀=%.3f)", lam, theta0_val)
  )
}))

p_interval <- ggplot(plot_data_5f, aes(x = factor(cat), y = prob, fill = label)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_viridis_d(name = "Parameters") +
  labs(
    title    = "Interval probabilities P(D ∈ [aₖ, bₖ)) under Exp(λ_d)",
    subtitle = "Category 1 = [0, 3mo), ..., Category 7 = [60mo, ∞)",
    x        = "Timegap category",
    y        = "Probability"
  ) +
  theme_minimal()
ggsave("troubleshoot/plot_interval_probs.png", p_interval,
       width = 10, height = 6, dpi = 150)
cat("Saved: troubleshoot/plot_interval_probs.png\n")

# Plot: observed vs theoretical transition matrix for categories 5 and 6
if (nrow(trans_obs) > 0) {
  # Focus on categories 5 and 6 (the only probabilistic transitions)
  obs_56 <- trans_obs[cat_prev %in% c(5, 6)]
  if (nrow(obs_56) > 0) {
    obs_56_tab <- obs_56[, .(n = .N), by = .(cat_prev, cat_curr)]
    obs_56_tab[, pct := round(n / sum(n) * 100, 1), by = cat_prev]
    obs_56_tab[, type := "Observed"]

    # Theoretical predictions at multiple λ_d values
    theo_56 <- rbindlist(lapply(c(0.005, 0.01, 0.02), function(lam) {
      rbindlist(lapply(5:6, function(j) {
        probs <- sapply(1:.N_CATS, function(k) .transition_prob(lam, j, k))
        nonzero <- which(probs > 0.001)
        data.table(
          cat_prev = j,
          cat_curr = nonzero,
          pct = round(probs[nonzero] * 100, 1),
          type = sprintf("Theory λ=%.3f", lam)
        )
      }))
    }))

    cat("\nCategory 5-6 transitions (observed vs theoretical):\n")
    print(obs_56_tab[, .(cat_prev, cat_curr, n, pct, type)])
    cat("\n")
    print(theo_56)
  }
}

cat("\n", strrep("=", 70), "\n")
cat("§5 SUMMARY\n")
cat(strrep("=", 70), "\n\n")
cat("Key questions answered by this section:\n")
cat("  1. What is the empirical category distribution? (5a)\n")
cat("  2. Do observed transitions match the sparse-matrix prediction? (5b)\n")
cat("  3. How much does the midpoint EMG distort vs interval probs? (5c)\n")
cat("  4. Are within-panel starts always in category 1? (5d)\n")
cat("  5. How many observations need floor vs donor imputation? (5e)\n")
cat("  6. Visual: how interval probs vary with λ_d (5f)\n")
cat("\nIf the observed transitions broadly match the theoretical sparse matrix,\n")
cat("the discrete/censored Exp(λ_d) model is appropriate for this data.\n")
cat("Proceed to implementation per troubleshoot/DIAGNOSIS.md and the plan\n")
cat("in .cg-docs/plans/2026-03-19-discrete-timegap-model.md\n")
