# ==============================================================================
# EM-tenure: Simulation-based diagnostic analysis
# ==============================================================================
# Supervisor-recommended analysis to understand model misspecification.
#
# Steps:
#   1. Load real QLFS data (post-ingest, never-worked excluded)
#   2. Simulate data from THREE parameter sets:
#      a) "Sim-EM": parameters from the current EM-tenure estimates
#      b) "Sim-Expected": plausible params from binary MLE on FULL sample
#      c) "Sim-Adjusted": same as (b) but recalibrated for the filtered
#         population (never-worked excluded) using empirical moments
#   3. Compare all simulated datasets against real data across 7 dimensions
#   4. Produce summary tables and diagnostic plots
#
# The goal is to understand WHAT the EM-tenure model is fitting and WHERE
# the model specification breaks down.
#
# Usage from project root:
#   source("EM-tenure/simulation_diagnostic.R")
# ==============================================================================

library(data.table)
library(ggplot2)

# --- Load EM-tenure module ---
source("EM-tenure/R/source_all.R")
source("EM-tenure/R/compare_distributions.R")

# --- Output directory ---
output_dir <- "EM-tenure/output"
fig_dir    <- file.path(output_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# §1  LOAD REAL DATA
# ==============================================================================

message("=== §1 Loading real QLFS data ===")

library(tidyverse)
source("scripts/ingest_data_3waves_SA.R")

df_qlfs <- df_qlfs |>
  mutate(across(where(haven::is.labelled), \(x) as.numeric(x)))

# Exclude never-worked (same filter as estimate_pipeline.R)
df_real <- df_qlfs |>
  filter(!(!is.na(neverworked1) & as.numeric(neverworked1) == 1)) |>
  filter(!(!is.na(neverworked2) & as.numeric(neverworked2) == 1)) |>
  filter(!(!is.na(neverworked3) & as.numeric(neverworked3) == 1))

# Ensure weight column
if (!"weight" %in% names(df_real)) {
  df_real$weight <- df_real$weight1
}

message(sprintf("Real data: n = %d", nrow(df_real)))


# ==============================================================================
# §2  DEFINE PARAMETER SETS
# ==============================================================================

message("\n=== §2 Defining parameter sets ===")

# --- Sim-EM: current EM-tenure estimates ---
# UPDATE THESE after the estimation pipeline finishes. The values below are
# from the iteration trace (EM has not converged yet — use best available).
# If you have a saved fit object, load it:
#   fit <- readRDS("output/results/fit_miscl_YYYYMMDD_HHMMSS.rds")
#   params_em <- fit$params
params_em <- list(
  alpha    = 0.775,
  theta1   = 0.50,    # implausible — this is the problem
  theta0   = 0.44,    # implausible
  pi       = 0.14,    # implausible
  sigma2_g = 0.12
)

# --- Sim-Expected: plausible parameters from simple binary model ---
# These come from the binary-sequence MLE (no duration data) and are
# consistent with the literature on South African labour markets.
params_expected <- list(
  alpha    = 0.47,
  theta1   = 0.95,
  theta0   = 0.05,
  pi       = 0.03,
  sigma2_g = 0.01    # small measurement error
)

# --- Sim-Adjusted: expected parameters recalibrated for filtered population ---
# The binary-sequence MLE was estimated on the FULL sample (including
# never-worked). The EM-tenure estimation excludes never-worked (~35% of obs,
# virtually all nonemployed). This changes the target population:
#   - alpha jumps from 0.47 to ~0.69 (employed count unchanged, nonemp halved)
#   - theta1 is nearly unchanged (never-worked don't enter theta1 calculation)
#   - theta0 roughly doubles (never-worked were ~55% of nonemployed, all
#     "stuck" in state 0, massively deflating the exit rate)
#   - pi unchanged (never-worked are correctly classified, don't drive pi)
# Values below are derived from empirical filtered transition rates with a
# simple misclassification correction: theta_adj = (naive - pi) / (1 - 2*pi)

# Compute empirical moments from df_real to anchor the adjustment
.alpha_filt <- mean(c(df_real$y1, df_real$y2, df_real$y3))
.tab12 <- table(df_real$y1, df_real$y2)
.tab23 <- table(df_real$y2, df_real$y3)
.naive_t1 <- mean(c(.tab12["1","1"] / sum(.tab12["1",]),
                     .tab23["1","1"] / sum(.tab23["1",])))
.naive_t0 <- mean(c(.tab12["0","1"] / sum(.tab12["0",]),
                     .tab23["0","1"] / sum(.tab23["0",])))
.pi_approx <- 0.03

params_adjusted <- list(
  alpha    = round(.alpha_filt, 3),
  theta1   = round((.naive_t1 - .pi_approx) / (1 - 2 * .pi_approx), 3),
  theta0   = round((.naive_t0 - .pi_approx) / (1 - 2 * .pi_approx), 3),
  pi       = 0.03,
  sigma2_g = 0.01
)

message(sprintf("  Empirical filtered alpha = %.3f", .alpha_filt))
message(sprintf("  Naive theta1 = %.3f, theta0 = %.3f", .naive_t1, .naive_t0))
message(sprintf("  Adjusted theta1 = %.3f, theta0 = %.3f (pi=%.2f correction)",
                params_adjusted$theta1, params_adjusted$theta0, .pi_approx))
rm(.alpha_filt, .tab12, .tab23, .naive_t1, .naive_t0, .pi_approx)

# Print all three
message("\nSim-EM parameters:")
message(paste(capture.output(str(params_em)), collapse = "\n"))
message("\nSim-Expected parameters (full sample):")
message(paste(capture.output(str(params_expected)), collapse = "\n"))
message("\nSim-Adjusted parameters (filtered population):")
message(paste(capture.output(str(params_adjusted)), collapse = "\n"))


# ==============================================================================
# §3  SIMULATE DATA
# ==============================================================================

message("\n=== §3 Simulating data ===")

n_sim <- nrow(df_real)
set.seed(2026)

message(sprintf("Simulating n = %d from Sim-EM parameters...", n_sim))
df_sim_em <- simulate_panel(
  n      = n_sim,
  alpha  = params_em$alpha,
  theta1 = params_em$theta1,
  theta0 = params_em$theta0,
  pi     = params_em$pi,
  sigma2_g = params_em$sigma2_g,
  discrete_timegap = TRUE,
  seed   = 42
)

message(sprintf("Simulating n = %d from Sim-Expected parameters...", n_sim))
df_sim_exp <- simulate_panel(
  n      = n_sim,
  alpha  = params_expected$alpha,
  theta1 = params_expected$theta1,
  theta0 = params_expected$theta0,
  pi     = params_expected$pi,
  sigma2_g = params_expected$sigma2_g,
  discrete_timegap = TRUE,
  seed   = 123
)

message(sprintf("Simulating n = %d from Sim-Adjusted parameters...", n_sim))
df_sim_adj <- simulate_panel(
  n      = n_sim,
  alpha  = params_adjusted$alpha,
  theta1 = params_adjusted$theta1,
  theta0 = params_adjusted$theta0,
  pi     = params_adjusted$pi,
  sigma2_g = params_adjusted$sigma2_g,
  discrete_timegap = TRUE,
  seed   = 456
)


# ==============================================================================
# §4  COMPARE: SUMMARY TABLES
# ==============================================================================

message("\n=== §4 Distributional comparisons ===")

# --- (a) Employment rates ---
cat("\n--- (a) Employment rates ---\n")
er <- rbindlist(list(
  compute_employment_rates(df_real,    "Real"),
  compute_employment_rates(df_sim_em,  "Sim-EM"),
  compute_employment_rates(df_sim_exp, "Sim-Expected"),
  compute_employment_rates(df_sim_adj, "Sim-Adjusted")
))
print(dcast(er, wave ~ dataset, value.var = "emp_rate"))

# --- (b) Transition matrices ---
cat("\n--- (b) Transition matrices (wave 1→2) ---\n")
tr <- rbindlist(list(
  compute_transitions(df_real,    "y1", "y2", "Real"),
  compute_transitions(df_sim_em,  "y1", "y2", "Sim-EM"),
  compute_transitions(df_sim_exp, "y1", "y2", "Sim-Expected"),
  compute_transitions(df_sim_adj, "y1", "y2", "Sim-Adjusted")
))
print(dcast(tr, from_state + to_state ~ dataset, value.var = "row_pct"))

# --- (c) Tenure distribution ---
cat("\n--- (c) Tenure distribution (pooled, employed obs) ---\n")
ten <- rbindlist(list(
  compute_tenure_stats(df_real,    "Real"),
  compute_tenure_stats(df_sim_em,  "Sim-EM"),
  compute_tenure_stats(df_sim_exp, "Sim-Expected"),
  compute_tenure_stats(df_sim_adj, "Sim-Adjusted")
))
print(ten)

# --- (d) Timegap category distribution ---
cat("\n--- (d) Timegap category distribution (pooled, nonemployed obs) ---\n")
tg <- rbindlist(list(
  compute_timegap_cat_dist(df_real,    "Real"),
  compute_timegap_cat_dist(df_sim_em,  "Sim-EM"),
  compute_timegap_cat_dist(df_sim_exp, "Sim-Expected"),
  compute_timegap_cat_dist(df_sim_adj, "Sim-Adjusted")
))
print(dcast(tg, cat ~ dataset, value.var = "pct", fill = 0))

# --- (e) Tenure increments ---
cat("\n--- (e) Tenure increments (emp→emp pairs) ---\n")
ti_real <- compute_tenure_increments(df_real,    "Real")
ti_em   <- compute_tenure_increments(df_sim_em,  "Sim-EM")
ti_exp  <- compute_tenure_increments(df_sim_exp, "Sim-Expected")
ti_adj  <- compute_tenure_increments(df_sim_adj, "Sim-Adjusted")
print(rbindlist(list(ti_real$stats, ti_em$stats, ti_exp$stats, ti_adj$stats)))

# --- (f) Timegap category transitions ---
cat("\n--- (f) Timegap category transitions (nonemp→nonemp) ---\n")
cat("Real:\n")
print(compute_timegap_transitions(df_real, "Real"))
cat("\nSim-EM:\n")
print(compute_timegap_transitions(df_sim_em, "Sim-EM"))
cat("\nSim-Expected:\n")
print(compute_timegap_transitions(df_sim_exp, "Sim-Expected"))
cat("\nSim-Adjusted:\n")
print(compute_timegap_transitions(df_sim_adj, "Sim-Adjusted"))

# --- (g) Three-wave sequence distribution ---
cat("\n--- (g) Three-wave sequence pattern distribution ---\n")
seq_dist <- rbindlist(list(
  compute_sequence_dist(df_real,    "Real"),
  compute_sequence_dist(df_sim_em,  "Sim-EM"),
  compute_sequence_dist(df_sim_exp, "Sim-Expected"),
  compute_sequence_dist(df_sim_adj, "Sim-Adjusted")
))
print(dcast(seq_dist, pattern ~ dataset, value.var = "pct", fill = 0))


# ==============================================================================
# §5  DIVERGENCE SUMMARY
# ==============================================================================

cat("\n--- Divergence summary ---\n")
div_em  <- compute_divergence_summary(df_real, df_sim_em,  "Real", "Sim-EM")
div_exp <- compute_divergence_summary(df_real, df_sim_exp, "Real", "Sim-Expected")
div_adj <- compute_divergence_summary(df_real, df_sim_adj, "Real", "Sim-Adjusted")
div_all <- rbindlist(list(div_em, div_exp, div_adj))
print(dcast(div_all, dimension ~ comparison, value.var = "value"))


# ==============================================================================
# §6  DIAGNOSTIC PLOTS
# ==============================================================================

message("\n=== §6 Saving diagnostic plots ===")

theme_diag <- theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# --- Plot 1: Tenure density (log-x) ---
tenure_pool <- rbindlist(list(
  extract_tenure_pool(df_real,    "Real"),
  extract_tenure_pool(df_sim_em,  "Sim-EM"),
  extract_tenure_pool(df_sim_exp, "Sim-Expected"),
  extract_tenure_pool(df_sim_adj, "Sim-Adjusted")
))

p1 <- ggplot(tenure_pool[tenure > 0.01 & tenure < 50],
             aes(x = tenure, color = dataset)) +
  geom_density(linewidth = 0.8) +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10, 20, 50),
                labels = c("0.1", "0.25", "0.5", "1", "2", "5", "10", "20", "50")) +
  labs(title = "Tenure distribution: Real vs Simulated",
       subtitle = "Log scale — note where the simulated modes fall",
       x = "Tenure (years, log scale)", y = "Density", color = NULL) +
  theme_diag
ggsave(file.path(fig_dir, "tenure_density_comparison.png"), p1,
       width = 10, height = 6, dpi = 150)
message("Saved: ", file.path(fig_dir, "tenure_density_comparison.png"))

# --- Plot 2: Tenure increment histogram ---
inc_pool <- rbindlist(list(
  data.table(increment = ti_real$increments, dataset = "Real"),
  data.table(increment = ti_em$increments,   dataset = "Sim-EM"),
  data.table(increment = ti_exp$increments,  dataset = "Sim-Expected"),
  data.table(increment = ti_adj$increments,  dataset = "Sim-Adjusted")
))
# Trim extreme increments for readability
inc_pool_trim <- inc_pool[abs(increment) < 5]

p2 <- ggplot(inc_pool_trim, aes(x = increment, fill = dataset)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80, alpha = 0.5,
                 position = "identity") +
  labs(title = "Tenure increment distribution (emp→emp pairs)",
       subtitle = "Model expects N(0, 2σ²_g) — check shape and spread",
       x = "tenure_t - tenure_{t-1} - 0.25 (years)",
       y = "Density", fill = NULL) +
  theme_diag
ggsave(file.path(fig_dir, "tenure_increment_comparison.png"), p2,
       width = 10, height = 6, dpi = 150)
message("Saved: ", file.path(fig_dir, "tenure_increment_comparison.png"))

# --- Plot 3: Timegap category bar chart ---
p3 <- ggplot(tg, aes(x = factor(cat), y = pct, fill = dataset)) +
  geom_col(position = "dodge") +
  labs(title = "Timegap category distribution (nonemployed obs)",
       subtitle = "Exp(λ_d) predicts monotonically decreasing — does the data agree?",
       x = "Timegap category", y = "Percent", fill = NULL) +
  theme_diag
ggsave(file.path(fig_dir, "timegap_category_comparison.png"), p3,
       width = 10, height = 6, dpi = 150)
message("Saved: ", file.path(fig_dir, "timegap_category_comparison.png"))

# --- Plot 4: Three-wave sequence pattern bar chart ---
p4 <- ggplot(seq_dist, aes(x = pattern, y = pct, fill = dataset)) +
  geom_col(position = "dodge") +
  labs(title = "Three-wave sequence pattern distribution",
       subtitle = "Which binary patterns does each parameter set match?",
       x = "Pattern (y1 y2 y3)", y = "Percent", fill = NULL) +
  theme_diag
ggsave(file.path(fig_dir, "sequence_pattern_comparison.png"), p4,
       width = 10, height = 6, dpi = 150)
message("Saved: ", file.path(fig_dir, "sequence_pattern_comparison.png"))

# --- Plot 5: Tenure CDF comparison ---
p5 <- ggplot(tenure_pool[tenure > 0.01 & tenure < 50],
             aes(x = tenure, color = dataset)) +
  stat_ecdf(linewidth = 0.7) +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10, 20, 50),
                labels = c("0.1", "0.25", "0.5", "1", "2", "5", "10", "20", "50")) +
  labs(title = "Tenure CDF: Real vs Simulated",
       subtitle = "The gap between curves shows the distributional mismatch",
       x = "Tenure (years, log scale)", y = "Cumulative probability",
       color = NULL) +
  theme_diag
ggsave(file.path(fig_dir, "tenure_cdf_comparison.png"), p5,
       width = 10, height = 6, dpi = 150)
message("Saved: ", file.path(fig_dir, "tenure_cdf_comparison.png"))


# ==============================================================================
# §7  INTERPRETATION GUIDE
# ==============================================================================

cat("\n", strrep("=", 70), "\n")
cat("INTERPRETATION GUIDE\n")
cat(strrep("=", 70), "\n\n")

cat("0. POPULATION MISMATCH NOTE:\n")
cat("   The binary MLE (Sim-Expected) was estimated on the FULL sample\n")
cat("   including never-worked individuals. The EM-tenure (Sim-EM) and the\n")
cat("   real data (df_real) EXCLUDE never-worked. Sim-Adjusted corrects\n")
cat("   Sim-Expected for this population difference using empirical moments\n")
cat("   from df_real. Key shifts: alpha 0.47→0.69, theta0 0.05→0.12.\n")
cat("   Compare Sim-Adjusted to Real (same population) for a fair benchmark.\n\n")

cat("1. TENURE DISTRIBUTION (Plots 1, 5):\n")
cat("   - Real data has median ~3.8 yrs, heavy right tail to 70 yrs\n")
cat("   - Sim-Expected (theta1=0.95, alpha=0.47): wrong population, very short spells\n")
cat("   - Sim-Adjusted (theta1~0.96, alpha~0.69): right population, still short spells\n")
cat("   - Sim-EM (theta1=0.50): slightly longer but still short spells\n")
cat("   → Neither matches. The exponential duration model cannot produce\n")
cat("     the observed tenure distribution. This is the core misspecification.\n")
cat("   → Sim-Adjusted is the fair comparison — same population as Real.\n\n")

cat("2. TRANSITION MATRICES (Table b):\n")
cat("   - Real (filtered): naive theta1~0.93, theta0~0.14\n")
cat("   - Sim-Adjusted should match these closely (same-population benchmark)\n")
cat("   - Sim-Expected understates theta0 (never-worked deflate it)\n")
cat("   - Sim-EM will show ~50/50 transitions (wrong)\n")
cat("   → If Sim-Adjusted matches transitions well, the binary model is\n")
cat("     correct for this population but duration data drags EM away.\n\n")

cat("3. SEQUENCE PATTERNS (Plot 4):\n")
cat("   - In filtered data: 111 dominates (~42%), 000 smaller (~18%)\n")
cat("   - Sim-Expected over-predicts 000 (wrong alpha)\n")
cat("   - Sim-Adjusted should match the pattern mix better\n")
cat("   - Sim-EM may redistribute mass toward mixed patterns\n")
cat("   → This reveals WHAT the EM is fitting vs what it should match.\n\n")

cat("4. TENURE INCREMENTS (Plot 2, Table e):\n")
cat("   - Under the model: N(0, 2*sigma2_g)\n")
cat("   - Real data may show heaping, skewness, or heavy tails\n")
cat("   → Increment distribution is less sensitive to the CTMC link\n")
cat("     (it depends on sigma2_g, not lambda_g). If this matches well,\n")
cat("     the measurement error model is reasonable.\n\n")

cat("5. TIMEGAP CATEGORIES (Plot 3):\n")
cat("   - Exp(lambda_d) predicts monotonically decreasing mass across cats\n")
cat("   - Real data may have non-monotone shape (mode at higher categories)\n")
cat("   → If non-monotone, the exponential assumption fails on both sides.\n\n")

cat("6. DIVERGENCE TABLE:\n")
cat("   - Lower = better match to real data\n")
cat("   - Compare Sim-EM vs Sim-Adjusted (not Sim-Expected) for a fair test\n")
cat("   - Sim-Expected vs Sim-Adjusted difference shows the population effect\n")
cat("   - The pattern of which dimensions each simulation wins reveals\n")
cat("     exactly what the EM is trading off.\n\n")

# Save the divergence table as CSV
fwrite(dcast(div_all, dimension ~ comparison, value.var = "value"),
       file.path(output_dir, "divergence_summary.csv"))
message("Saved: ", file.path(output_dir, "divergence_summary.csv"))

message("\n=== Simulation diagnostic complete ===")
