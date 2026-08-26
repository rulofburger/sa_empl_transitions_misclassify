# ==============================================================================
# analyse_ee_non_clock_tenure.R
# Created: 2026-05-12
#
# PURPOSE
#   Approximately 65 % of Employed-to-Employed (E->E) continuation pairs in
#   the QLFS panel satisfy the deterministic tenure clock: reported tenure at
#   wave t+1 is exactly 3 months higher than at wave t.  This script analyses
#   the *level* distribution of tenure for the remaining ~35 % ("non-clock"
#   pairs), asking two questions:
#
#     Q1. Does the destination tenure of non-clock E->E pairs follow an
#         exponential distribution?
#
#     Q2. If so, is the exponential rate equal to that of the unconditional
#         wave-1 employed tenure distribution?
#
#   A positive answer to both questions supports the model assumption that
#   contaminated tenure observations are i.i.d. draws from the stationary
#   population tenure distribution, independent of the worker's true state.
#
# GROUPS
#   A — Non-clock E->E destination tenure
#       Tenure at wave t+1 for E->E pairs where |delta - 3 months| >= CLOCK_TOL.
#       This is the "contaminated" observation in the epsilon model.
#   B — Unconditional wave-1 employed tenure
#       Tenure for all workers employed in wave 1 (t1 > 0).
#
# OUTPUTS
#   Figures: output/figures/tenure_ee_distribution/
#   Tables:  output/tables/tenure_ee_distribution/
#
# USAGE
#   Rscript scripts/analyse_ee_non_clock_tenure.R
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(tidyverse)

})

# --- Output directories -------------------------------------------------------
fig_dir <- here("output", "figures", "tenure_ee_distribution")
tab_dir <- here("output", "tables", "tenure_ee_distribution")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

.pw  <- 8L    # plot width (inches)
.ph  <- 5L    # plot height (inches)
.dpi <- 150L  # plot resolution

# Tolerance for "exact +3-month clock" (months); mirrors diagnose script.
CLOCK_TOL_MO <- 0.01

# ==============================================================================
# 1. Load data
# ==============================================================================
source(here("scripts", "ingest_data_3waves_SA.R"))
if (!exists("df_qlfs")) stop("df_qlfs not found after sourcing ingest script.")
setDT(df_qlfs)

# Tenure in months (ingest divides raw months by 12, so we reverse here).
df_qlfs[, t1_mo := tenure1 * 12]
df_qlfs[, t2_mo := tenure2 * 12]
df_qlfs[, t3_mo := tenure3 * 12]

cat(sprintf("\nPanel loaded: N = %s individuals\n\n", format(nrow(df_qlfs), big.mark = ",")))




df_qlfs |> 
    select(y1, y2, y3, t1_mo, t2_mo, t3_mo, weight) |> 
    filter(y1 == 1L & y2 == 1L & y3 == 1L) |>
    mutate(diff12 = t2_mo - t1_mo, diff23 = t3_mo - t2_mo) |>
    filter(diff12 == 3 & diff23 == 3) |>
    head(10)

df_qlfs |> 
    select(y1, y2, y3, t1_mo, t2_mo, t3_mo, weight) |> 
    filter(y1 == 1L & y2 == 1L & y3 == 1L) |>
    mutate(diff12 = t2_mo - t1_mo, diff23 = t3_mo - t2_mo) |>
    filter(diff12 != 3 & diff23 != 3) |>
    head(30)

x <- df_qlfs |> 
    select(y1, y2, y3, t1_mo, t2_mo, t3_mo, weight) |> 
    filter(y1 == 1L & y2 == 1L & y3 == 1L) |>
    mutate(diff12 = t2_mo - t1_mo, diff23 = t3_mo - t2_mo, diff13 = t3_mo - t1_mo) |>
    filter(diff13 == 6 & diff12 != 3) |>
    pull(t2_mo)

rate_est <- 1 / mean(x)
n_x      <- length(x)
ps       <- (seq_len(n_x) - 0.5) / n_x          # Hazen plotting positions
qq_df    <- data.table(
  empirical   = sort(x),
  theoretical = qexp(ps, rate = rate_est)
)

ggplot(qq_df, aes(x = theoretical, y = empirical)) +
  geom_point(alpha = 0.25, size = 0.6, colour = "#4477AA") +
  geom_abline(slope = 1, intercept = 0, colour = "#E41A1C", linetype = "dashed") +
  labs(
    title    = "Exponential Q-Q plot: non-clock E->E (diff13 == 6, diff12 != 3)",
    subtitle = sprintf("lambda_hat = 1/mean(x) = %.5f  (mean = %.1f months)", rate_est, 1/rate_est),
    x        = "Theoretical quantile [Exp(lambda_hat)]",
    y        = "Observed tenure wave 2 (months)"
  ) +
  theme_minimal(base_size = 12)

# Histogram with theoretical Exp(rate_est) density overlay (base R)
hist(x,
     breaks   = "FD",                   # Freedman-Diaconis bin width
     freq     = FALSE,                  # density scale so curve() is comparable
     col      = "#4477AA",
     border   = "white",
     main     = sprintf(
       "Non-clock E->E: histogram vs Exp(lambda_hat)\nlambda_hat = %.5f  (mean = %.1f months)",
       rate_est, 1 / rate_est
     ),
     xlab     = "Tenure wave 2 (months)",
     ylab     = "Density",
     xlim     = c(0, quantile(x, 0.99)))  # truncate at 99th pct for readability
curve(dexp(x, rate = rate_est),
      add   = TRUE,
      col   = "#E41A1C",
      lwd   = 2)
legend("topright",
       legend = c("Observed", sprintf("Exp(lambda = %.5f)", rate_est)),
       fill   = c("#4477AA", NA),
       border = c("white", NA),
       lty    = c(NA, 1),
       lwd    = c(NA, 2),
       col    = c(NA, "#E41A1C"),
       bty    = "n")

# ==============================================================================
# 2. Build E->E continuation pairs and classify clock vs non-clock
# ==============================================================================
# Pool waves 1->2 and 2->3.  Record source tenure (from_mo), destination
# tenure (to_mo), the increment (delta_mo), and the survey weight.
ee12 <- df_qlfs[
  y1 == 1L & y2 == 1L & !is.na(t1_mo) & !is.na(t2_mo),
  .(from_mo = t1_mo, to_mo = t2_mo, delta_mo = t2_mo - t1_mo, w = weight)
]
ee23 <- df_qlfs[
  y2 == 1L & y3 == 1L & !is.na(t2_mo) & !is.na(t3_mo),
  .(from_mo = t2_mo, to_mo = t3_mo, delta_mo = t3_mo - t2_mo, w = weight)
]
ee <- rbind(ee12, ee23)

ee[, is_clock := abs(delta_mo - 3) < CLOCK_TOL_MO]

n_ee       <- nrow(ee)
n_clock    <- sum(ee$is_clock, na.rm = TRUE)
n_nonclock <- sum(!ee$is_clock, na.rm = TRUE)
pct_clock  <- 100 * n_clock / n_ee

cat(sprintf("E->E pairs (pooled waves 1-2 and 2-3):  %s\n", format(n_ee, big.mark = ",")))
cat(sprintf("  Clock (|delta - 3| < %.2f mo):  %s  (%.1f%%)\n",
            CLOCK_TOL_MO, format(n_clock, big.mark = ","), pct_clock))
cat(sprintf("  Non-clock:                      %s  (%.1f%%)\n\n",
            format(n_nonclock, big.mark = ","), 100 - pct_clock))

# ==============================================================================
# 3. Define analysis groups
# ==============================================================================
# Group A: destination tenure for non-clock pairs (the contaminated observation).
# Group B: wave-1 employed tenure (the unconditional reference distribution).
# Both restricted to strictly positive tenure (well-defined exponential support).
grp_A <- ee[is_clock == FALSE & to_mo > 0, .(tenure_mo = to_mo, w)]
grp_B <- df_qlfs[y1 == 1L & !is.na(t1_mo) & t1_mo > 0, .(tenure_mo = t1_mo, w = weight)]

cat(sprintf("Group A — non-clock E->E destination tenure:   N = %s\n",
            format(nrow(grp_A), big.mark = ",")))
cat(sprintf("Group B — wave-1 unconditional employed:        N = %s\n\n",
            format(nrow(grp_B), big.mark = ",")))

# ==============================================================================
# 4. Descriptive statistics
# ==============================================================================
describe_group <- function(dt, label) {
  x  <- dt$tenure_mo
  wt <- dt$w / sum(dt$w, na.rm = TRUE)
  data.table(
    group     = label,
    N         = nrow(dt),
    mean_w    = round(sum(x * wt, na.rm = TRUE), 1),
    mean_uw   = round(mean(x, na.rm = TRUE), 1),
    median_uw = round(median(x, na.rm = TRUE), 1),
    sd_uw     = round(sd(x, na.rm = TRUE), 1),
    pct_lt6   = round(100 * mean(x < 6,   na.rm = TRUE), 1),
    pct_lt12  = round(100 * mean(x < 12,  na.rm = TRUE), 1),
    pct_ge36  = round(100 * mean(x >= 36, na.rm = TRUE), 1)
  )
}

desc <- rbind(
  describe_group(grp_A, "Non-clock E->E"),
  describe_group(grp_B, "Wave-1 employed")
)
fwrite(desc, file.path(tab_dir, "01_descriptive_stats.csv"))

cat("---- Descriptive statistics (months) ----\n")
print(desc)
cat("\n")

# ==============================================================================
# 5. Exponential MLE fits
# ==============================================================================
# For Exp(lambda): MLE lambda_hat = 1 / mean(x).
# Weighted MLE:    lambda_hat_w  = 1 / weighted_mean(x, w).
# Delta-method SE: SE(lambda_hat) = lambda_hat / sqrt(n)  (unweighted only).
fit_exp <- function(dt, label) {
  x  <- dt$tenure_mo
  w  <- dt$w
  n  <- length(x)
  lambda_uw <- 1 / mean(x, na.rm = TRUE)
  se_uw     <- lambda_uw / sqrt(n)
  lambda_wt <- 1 / (sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE))
  data.table(
    group     = label,
    N         = n,
    lambda_uw = lambda_uw,
    se_uw     = se_uw,
    ci95_lo   = lambda_uw - 1.96 * se_uw,
    ci95_hi   = lambda_uw + 1.96 * se_uw,
    mean_uw   = 1 / lambda_uw,
    lambda_wt = lambda_wt,
    mean_wt   = 1 / lambda_wt
  )
}

fit_A <- fit_exp(grp_A, "Non-clock E->E")
fit_B <- fit_exp(grp_B, "Wave-1 employed")
fits  <- rbind(fit_A, fit_B)
fwrite(fits, file.path(tab_dir, "02_exponential_fits.csv"))

cat("---- Exponential MLE fits ----\n")
cat("  Unweighted MLE (delta-method SE):\n")
cat(sprintf("    Non-clock:   lambda = %.5f  mean = %.1f mo  95%% CI [%.5f, %.5f]\n",
            fit_A$lambda_uw, fit_A$mean_uw, fit_A$ci95_lo, fit_A$ci95_hi))
cat(sprintf("    Wave-1:      lambda = %.5f  mean = %.1f mo  95%% CI [%.5f, %.5f]\n",
            fit_B$lambda_uw, fit_B$mean_uw, fit_B$ci95_lo, fit_B$ci95_hi))
cat(sprintf("    Rate ratio (A/B):  %.4f\n\n", fit_A$lambda_uw / fit_B$lambda_uw))

cat("  Survey-weighted MLE:\n")
cat(sprintf("    Non-clock:   lambda_w = %.5f  mean_w = %.1f mo\n",
            fit_A$lambda_wt, fit_A$mean_wt))
cat(sprintf("    Wave-1:      lambda_w = %.5f  mean_w = %.1f mo\n",
            fit_B$lambda_wt, fit_B$mean_wt))
cat(sprintf("    Rate ratio (A/B):  %.4f\n\n", fit_A$lambda_wt / fit_B$lambda_wt))

# ==============================================================================
# 6. Likelihood-ratio test for equal exponential rates (H0: lambda_A = lambda_B)
# ==============================================================================
# Under H0 the pooled MLE is lambda_0 = (n_A + n_B) / (sum_A x + sum_B x).
# LRT = 2 * (logL_unrestricted - logL_restricted) ~ chi2(1) under H0.
lrt_exp_equal_rates <- function(xA, xB) {
  nA <- length(xA);  nB <- length(xB)
  lA <- 1 / mean(xA);  lB <- 1 / mean(xB)
  l0 <- (nA + nB) / (sum(xA) + sum(xB))
  logL_unr <- nA * log(lA) - lA * sum(xA) + nB * log(lB) - lB * sum(xB)
  logL_res  <- (nA + nB) * log(l0) - l0 * (sum(xA) + sum(xB))
  LRT  <- 2 * (logL_unr - logL_res)
  pval <- pchisq(LRT, df = 1, lower.tail = FALSE)
  list(lambda_A = lA, lambda_B = lB, lambda_pool = l0,
       rate_ratio = lA / lB, LRT = LRT, df = 1, pval = pval)
}

lrt <- lrt_exp_equal_rates(grp_A$tenure_mo, grp_B$tenure_mo)

cat("---- LR test: H0: lambda_A = lambda_B ----\n")
cat(sprintf("  lambda_A (non-clock):  %.5f\n", lrt$lambda_A))
cat(sprintf("  lambda_B (wave-1):     %.5f\n", lrt$lambda_B))
cat(sprintf("  Rate ratio (A/B):      %.4f\n", lrt$rate_ratio))
cat(sprintf("  LRT = %.2f  df = 1  p = %.4g\n", lrt$LRT, lrt$pval))
cat("  Note: with N > 100,000 the LRT has power to reject trivially small\n")
cat("  differences.  The rate ratio and visual diagnostics are more\n")
cat("  substantively informative than the p-value.\n\n")

lrt_tbl <- data.table(
  lambda_A    = lrt$lambda_A,   lambda_B    = lrt$lambda_B,
  lambda_pool = lrt$lambda_pool, rate_ratio  = lrt$rate_ratio,
  LRT         = lrt$LRT,        df          = lrt$df,
  pval        = lrt$pval
)
fwrite(lrt_tbl, file.path(tab_dir, "03_lrt_equal_rates.csv"))

# ==============================================================================
# 7. Kolmogorov-Smirnov tests
# ==============================================================================
# KS test 1: Is Group A exponential (with its own fitted rate)?
ks_A  <- ks.test(grp_A$tenure_mo, "pexp", rate = fit_A$lambda_uw)
# KS test 2: Is Group B exponential (with its own fitted rate)?
ks_B  <- ks.test(grp_B$tenure_mo, "pexp", rate = fit_B$lambda_uw)
# KS test 3: Do A and B share the same distribution?
ks_AB <- ks.test(grp_A$tenure_mo, grp_B$tenure_mo)

cat("---- Kolmogorov-Smirnov tests (unweighted) ----\n")
cat("  Note: KS p-values for large N will reject tiny deviations.  The D\n")
cat("  statistic (max absolute CDF distance) is the key effect-size metric.\n")
cat(sprintf("  A vs Exp(lambda_A):   D = %.4f  p = %.4g\n", ks_A$statistic,  ks_A$p.value))
cat(sprintf("  B vs Exp(lambda_B):   D = %.4f  p = %.4g\n", ks_B$statistic,  ks_B$p.value))
cat(sprintf("  A vs B (two-sample):  D = %.4f  p = %.4g\n\n", ks_AB$statistic, ks_AB$p.value))

ks_tbl <- data.table(
  test    = c("Non-clock vs Exp(lambda_A)", "Wave-1 vs Exp(lambda_B)", "Non-clock vs Wave-1"),
  D       = c(ks_A$statistic, ks_B$statistic, ks_AB$statistic),
  p_value = c(ks_A$p.value,   ks_B$p.value,   ks_AB$p.value)
)
fwrite(ks_tbl, file.path(tab_dir, "04_ks_tests.csv"))

# ==============================================================================
# 8. Robustness: source tenure of non-clock pairs
# ==============================================================================
# If contamination enters at the destination wave (wave t+1 is a random draw),
# the SOURCE tenure (wave t) should behave like genuine employed tenure and
# should NOT look like the unconditional distribution.  Comparing source vs
# destination rate provides a rough check on this asymmetry.
grp_A_src <- ee[is_clock == FALSE & from_mo > 0, .(tenure_mo = from_mo, w)]
fit_A_src  <- fit_exp(grp_A_src, "Non-clock E->E source tenure")
ks_src_vs_B <- ks.test(grp_A_src$tenure_mo, grp_B$tenure_mo)

cat("---- Robustness: source (wave t) tenure of non-clock pairs ----\n")
cat(sprintf("  N = %s\n", format(nrow(grp_A_src), big.mark = ",")))
cat(sprintf("  lambda (source) = %.5f  (mean = %.1f mo)\n",
            fit_A_src$lambda_uw, fit_A_src$mean_uw))
cat(sprintf("  lambda (dest)   = %.5f  (mean = %.1f mo)\n",
            fit_A$lambda_uw, fit_A$mean_uw))
cat(sprintf("  KS D (source vs wave-1 unconditional): %.4f  p = %.4g\n\n",
            ks_src_vs_B$statistic, ks_src_vs_B$p.value))

rob_tbl <- data.table(
  measure    = c("Source tenure lambda (uw)", "Dest tenure lambda (uw)", "Wave-1 lambda (uw)"),
  lambda     = c(fit_A_src$lambda_uw, fit_A$lambda_uw, fit_B$lambda_uw),
  mean_mo    = c(fit_A_src$mean_uw,   fit_A$mean_uw,   fit_B$mean_uw),
  ks_vs_B_D  = c(ks_src_vs_B$statistic, ks_AB$statistic, NA_real_),
  ks_vs_B_p  = c(ks_src_vs_B$p.value,   ks_AB$p.value,   NA_real_)
)
fwrite(rob_tbl, file.path(tab_dir, "05_robustness_source_vs_dest.csv"))

# ==============================================================================
# 9. Visual diagnostics
# ==============================================================================
# Truncate at the 99th percentile for readability.
x_max  <- quantile(c(grp_A$tenure_mo, grp_B$tenure_mo), 0.99, na.rm = TRUE)
x_grid <- seq(0.01, x_max, length.out = 600)

# ---- 9a. ECDF vs fitted exponential CDF ----
ecdf_A <- ecdf(grp_A$tenure_mo)
ecdf_B <- ecdf(grp_B$tenure_mo)

cdf_df <- rbind(
  data.table(x = x_grid, y = ecdf_A(x_grid),
             group = "Non-clock E->E", series = "Empirical CDF"),
  data.table(x = x_grid, y = ecdf_B(x_grid),
             group = "Wave-1 employed", series = "Empirical CDF"),
  data.table(x = x_grid, y = pexp(x_grid, rate = fit_A$lambda_uw),
             group = "Non-clock E->E", series = "Fitted Exp(lambda)"),
  data.table(x = x_grid, y = pexp(x_grid, rate = fit_B$lambda_uw),
             group = "Wave-1 employed", series = "Fitted Exp(lambda)")
)

p_cdf <- ggplot(cdf_df, aes(x = x, y = y, colour = group, linetype = series)) +
  geom_line(linewidth = 0.85) +
  scale_colour_manual(
    values = c("Non-clock E->E" = "#E41A1C", "Wave-1 employed" = "#377EB8")
  ) +
  scale_linetype_manual(
    values = c("Empirical CDF" = "solid", "Fitted Exp(lambda)" = "dashed")
  ) +
  labs(
    title    = "Empirical CDFs vs. MLE exponential fits",
    subtitle = sprintf(
      "Non-clock: lambda=%.4f (mean %.0f mo) | Wave-1: lambda=%.4f (mean %.0f mo)",
      fit_A$lambda_uw, fit_A$mean_uw, fit_B$lambda_uw, fit_B$mean_uw
    ),
    x = "Tenure (months)", y = "Cumulative probability",
    colour = NULL, linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "01_ecdf_vs_fitted_exp.png"),
       p_cdf, width = .pw, height = .ph, dpi = .dpi)

# ---- 9b. Histograms with exponential density overlay ----
exp_dens <- rbind(
  data.table(x = x_grid, density = dexp(x_grid, rate = fit_A$lambda_uw),
             group = "Non-clock E->E destination tenure"),
  data.table(x = x_grid, density = dexp(x_grid, rate = fit_B$lambda_uw),
             group = "Wave-1 unconditional employed tenure")
)

hist_df <- rbind(
  data.table(tenure_mo = grp_A[tenure_mo <= x_max]$tenure_mo,
             group = "Non-clock E->E destination tenure"),
  data.table(tenure_mo = grp_B[tenure_mo <= x_max]$tenure_mo,
             group = "Wave-1 unconditional employed tenure")
)

p_hist <- ggplot(hist_df, aes(x = tenure_mo)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 3,
                 fill = "#4477AA", colour = "white", alpha = 0.65) +
  geom_line(data = exp_dens, aes(x = x, y = density),
            colour = "#E41A1C", linewidth = 0.9) +
  facet_wrap(~group, ncol = 1, scales = "free_y") +
  labs(
    title    = "Tenure distribution vs. MLE-fitted exponential",
    subtitle = "Red curve: Exp(lambda_hat). Binwidth = 3 months. Truncated at 99th percentile.",
    x = "Tenure (months)", y = "Density"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(fig_dir, "02_histogram_vs_exp.png"),
       p_hist, width = .pw, height = .ph + 2L, dpi = .dpi)

# ---- 9c. Exponential Q-Q plots ----
# For Exp(lambda), the p-th quantile is Q(p) = -log(1-p) / lambda.
make_qq <- function(x, lambda, label) {
  xs <- sort(x[!is.na(x) & x > 0])
  n  <- length(xs)
  ps <- (seq_len(n) - 0.5) / n
  data.table(empirical = xs, theoretical = qexp(ps, rate = lambda), group = label)
}

qq_A <- make_qq(grp_A$tenure_mo, fit_A$lambda_uw, "Non-clock E->E destination tenure")
qq_B <- make_qq(grp_B$tenure_mo, fit_B$lambda_uw, "Wave-1 unconditional employed tenure")
qq   <- rbind(qq_A, qq_B)
qq_max <- quantile(qq$empirical, 0.99, na.rm = TRUE)

p_qq <- ggplot(qq[empirical <= qq_max], aes(x = theoretical, y = empirical)) +
  geom_point(alpha = 0.15, size = 0.4, colour = "#4477AA") +
  geom_abline(slope = 1, intercept = 0, colour = "#E41A1C", linetype = "dashed") +
  facet_wrap(~group, ncol = 2, scales = "free") +
  labs(
    title    = "Exponential Q-Q plots",
    subtitle = "Dashed diagonal: perfect exponential fit. Truncated at 99th percentile.",
    x = "Theoretical quantile [Exp(lambda_hat)]", y = "Observed tenure (months)"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(fig_dir, "03_qq_plots.png"),
       p_qq, width = .pw + 2L, height = .ph, dpi = .dpi)

# ---- 9d. Log-survival plot: exponential => log(1-F(x)) linear in x ----
make_logsurv <- function(x, lambda_hat, label) {
  xs <- sort(x[!is.na(x) & x > 0 & x <= x_max])
  Fx <- ecdf(xs)(xs)
  data.table(
    x         = xs,
    empirical = log(pmax(1 - Fx, 1e-6)),
    reference = -lambda_hat * xs,
    group     = label
  )
}

ls_A <- make_logsurv(grp_A$tenure_mo, fit_A$lambda_uw, "Non-clock E->E destination tenure")
ls_B <- make_logsurv(grp_B$tenure_mo, fit_B$lambda_uw, "Wave-1 unconditional employed tenure")

ls_long <- rbind(
  melt(ls_A, id.vars = c("x", "group"),
       measure.vars = c("empirical", "reference"), variable.name = "series"),
  melt(ls_B, id.vars = c("x", "group"),
       measure.vars = c("empirical", "reference"), variable.name = "series")
)
ls_long[, series := fifelse(
  series == "reference", "Exp benchmark (-lambda*x)", "Empirical log-survival"
)]

p_ls <- ggplot(ls_long, aes(x = x, y = value, colour = series, linetype = series)) +
  geom_line(linewidth = 0.75) +
  scale_colour_manual(
    values = c("Empirical log-survival" = "#4477AA",
               "Exp benchmark (-lambda*x)" = "#E41A1C")
  ) +
  scale_linetype_manual(
    values = c("Empirical log-survival" = "solid",
               "Exp benchmark (-lambda*x)" = "dashed")
  ) +
  facet_wrap(~group, ncol = 2, scales = "free") +
  labs(
    title    = "Log-survival plot (exponential goodness-of-fit check)",
    subtitle = "Linearity in x implies exponential distribution. Red dashed: -lambda*x benchmark.",
    x = "Tenure (months)", y = "log(1 - F(x))",
    colour = NULL, linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "04_log_survival.png"),
       p_ls, width = .pw + 2L, height = .ph, dpi = .dpi)

# ---- 9e. Overlay: source vs destination tenure for non-clock pairs ----
# Shows whether contamination is asymmetric (destination different from source).
ovl_df <- rbind(
  data.table(tenure_mo = grp_A_src[tenure_mo <= x_max]$tenure_mo,
             series = "Non-clock source (wave t)"),
  data.table(tenure_mo = grp_A[tenure_mo <= x_max]$tenure_mo,
             series = "Non-clock destination (wave t+1)"),
  data.table(tenure_mo = grp_B[tenure_mo <= x_max]$tenure_mo,
             series = "Wave-1 unconditional")
)

p_ovl <- ggplot(ovl_df, aes(x = tenure_mo, colour = series)) +
  stat_ecdf(linewidth = 0.8) +
  scale_colour_manual(
    values = c(
      "Non-clock source (wave t)"        = "#FF7F00",
      "Non-clock destination (wave t+1)" = "#E41A1C",
      "Wave-1 unconditional"             = "#377EB8"
    )
  ) +
  labs(
    title    = "Empirical CDFs: source vs destination tenure for non-clock pairs",
    subtitle = "If contamination is at wave t+1, destination (red) should match wave-1 (blue); source (orange) should not.",
    x = "Tenure (months)", y = "Cumulative probability",
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "05_source_vs_dest_vs_unconditional.png"),
       p_ovl, width = .pw, height = .ph, dpi = .dpi)

# ==============================================================================
# 10. Summary
# ==============================================================================
cat("========================================================\n")
cat(" SUMMARY\n") 
cat("========================================================\n")
cat(sprintf("  Non-clock share:        %.1f%% of all E->E pairs\n", 100 - pct_clock))
cat("\n  Exponential rates (unweighted MLE):\n")
cat(sprintf("    lambda_A (non-clock):  %.5f  (mean = %.1f months)\n",
            fit_A$lambda_uw, fit_A$mean_uw))
cat(sprintf("    lambda_B (wave-1):     %.5f  (mean = %.1f months)\n",
            fit_B$lambda_uw, fit_B$mean_uw))
cat(sprintf("    Rate ratio (A/B):      %.4f\n", lrt$rate_ratio))
cat(sprintf("    LRT p-value:           %.4g\n", lrt$pval))
cat("\n  KS D statistics:\n")
cat(sprintf("    Non-clock vs Exp:      D = %.4f  (p = %.4g)\n",
            ks_A$statistic, ks_A$p.value))
cat(sprintf("    Wave-1 vs Exp:         D = %.4f  (p = %.4g)\n",
            ks_B$statistic, ks_B$p.value))
cat(sprintf("    Non-clock vs Wave-1:   D = %.4f  (p = %.4g)\n",
            ks_AB$statistic, ks_AB$p.value))
cat(sprintf("\n  Figures written to: %s\n", fig_dir))
cat(sprintf("  Tables  written to: %s\n\n", tab_dir))
