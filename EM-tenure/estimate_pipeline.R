# ==============================================================================
# EM-tenure: End-to-end estimation pipeline
# ==============================================================================
# Created: 2026-04-29
# This script:
#   1. Sources the EM-tenure module
#   2. Sources the data ingestion script
#   3. Prepares the data (rename weight column if needed)
#   4. Runs the EM (non-stationary and stationary)
#   5. Saves each fit as a timestamped .rds in output/results/
#   6. Appends a one-row summary per model to output/results/run_summary.csv
#
# Usage from project root:
#   source("EM-tenure/estimate_pipeline.R")
# ==============================================================================

# Verify working directory: all paths are relative to the project root.
if (!file.exists("EM-tenure/R/source_all.R")) {
  stop("estimate_pipeline.R must be sourced from the project root. ",
       "Expected to find 'EM-tenure/R/source_all.R' relative to cwd.")
}

# --- Load EM-tenure module ---
source("EM-tenure/R/source_all.R")

# --- Ingest data ---
# Always re-ingest to avoid stale data from a previous session.
if (!file.exists("data/raw/df_qlfs_A.rds")) {
  stop("Missing input data: 'data/raw/df_qlfs_A.rds'. ",
       "Obtain the QLFS data and place it at that path before running the pipeline.")
}
library("tidyverse")
source("scripts/ingest_data_3waves_SA.R")

df_qlfs <- df_qlfs |>
  mutate(across(where(haven::is.labelled), \(x) as.numeric(x)))

# --- Never-worked: impute timegap as (age - 16) * 12 months, mapped to category ---
# Never-worked individuals have no previous employment. Their nonemployment
# duration is at least (age - 16) years. We convert to months and bin into
# the standard 7 timegap categories. This avoids excluding them and biasing
# the sample toward shorter nonemployment spells.
.nw_impute_timegap_cat <- function(age_years) {
  dur_months <- (age_years - 16) * 12
  dur_months <- pmax(dur_months, 0)  # safety: age < 16 → 0 months
  # Standard QLFS timegap bins (in months):
  # 1: [0,3), 2: [3,6), 3: [6,9), 4: [9,12), 5: [12,36), 6: [36,60), 7: [60,Inf)
  cut_vals <- c(0, 3, 6, 9, 12, 36, 60, Inf)
  as.integer(cut(dur_months, breaks = cut_vals, right = FALSE, include.lowest = TRUE))
}

for (.wave in 1:3) {
  nw_col  <- paste0("neverworked", .wave)
  tc_col  <- paste0("timegap_cat", .wave)
  age_col <- paste0("age", .wave)

  if (nw_col %in% names(df_qlfs) && age_col %in% names(df_qlfs)) {
    is_nw <- !is.na(df_qlfs[[nw_col]]) & as.numeric(df_qlfs[[nw_col]]) == 1
    if (any(is_nw)) {
      df_qlfs[[tc_col]][is_nw] <- .nw_impute_timegap_cat(df_qlfs[[age_col]][is_nw])
    }
  }
}
rm(.wave, .nw_impute_timegap_cat)
# --- Prepare data ---
# Ensure weight column exists (the ingest script uses weight1/weight2/weight3;
# we pick weight1 as the baseline weight).
if (!"weight" %in% names(df_qlfs)) {
  df_qlfs$weight <- df_qlfs$weight1
}

# Durations should already be in years and imputed from the ingest script.
# Verify columns exist (discrete mode: timegap_cat1-3 instead of timegap1-3)
required_cols <- c("y1", "y2", "y3",
                   "tenure1", "tenure2", "tenure3",
                   "timegap_cat1", "timegap_cat2", "timegap_cat3",
                   "weight")
stopifnot(all(required_cols %in% names(df_qlfs)))

# --- Run EM ---
message("=== Estimating EM models ===")

# NOTE: The ad-hoc zero-duration filter has been removed. The ingest script
# now performs nearest-non-zero imputation (Issue 1 fix), so all tenure
# values for employed waves and all timegap_cat values are valid before
# reaching this pipeline. The never-worked filter above removes individuals
# who were never employed, eliminating the source of structural zeros.

# Default initial parameters (used as fallback when no prior run exists for a model).
custom_init      <- init_params(df_qlfs, discrete_timegap = TRUE, linked = FALSE)
custom_init_rho  <- init_params_rho(df_qlfs, linked = FALSE)

# Per-model warm-start: load the most recent saved result for each of the 8
# model variants. Each model gets its own param object so that, e.g., the
# linked models warm-start from linked-converged params and the rho models
# warm-start from rho-converged params (which include the rho parameter).
# Falls back to init_params() / init_params_rho() defaults when no prior
# result exists for a given model.
# load_warm_starts() is defined in utils.R (already sourced via source_all.R above).
message("\n=== Loading per-model warm-start parameters ===")
ws <- load_warm_starts("output/results", df_qlfs, verbose = TRUE)

message("Warm-start status:")
for (.nm in c("miscl", "miscl_stationary", "miscl_linked", "miscl_stationary_linked",
              "rho", "rho_stationary", "rho_linked", "rho_stationary_linked",
              "eps", "eps_stationary", "eps_linked", "eps_stationary_linked")) {
  message(sprintf("  %-30s: %s", .nm, if (!is.null(ws[[.nm]])) "warm-start loaded" else "INIT_PARAMS (no prior run)"))
}
rm(.nm)

# --- Results storage setup ---
# run_id ties together the .rds files and the summary CSV row for this run.
# NOTE: The EM algorithm is fully deterministic given fixed data and starting
# parameters (no random operations in e-step, m-step, or driver). set.seed()
# is not needed here.
run_id     <- format(Sys.time(), "%Y%m%d_%H%M%S")
results_dir <- "output/results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
summary_csv <- file.path(results_dir, "run_summary.csv")

message("=== Estimating non-stationary model (free) ===")
t1_fit_em <- Sys.time()
fit_em <- em_fit_tenure(
  df                = df_qlfs,
  params0           = ws$miscl %||% custom_init,
  stationary        = FALSE,
  linked            = FALSE,
  discrete_timegap  = TRUE,
  verbose           = 2L
)
t2_fit_em <- Sys.time()

message("=== Estimating stationary model (free) ===")
t1_fit_stationary_em <- Sys.time()
fit_stationary_em <- em_fit_tenure(
  df                = df_qlfs,
  params0           = ws$miscl_stationary %||% custom_init,
  stationary        = TRUE,
  linked            = FALSE,
  discrete_timegap  = TRUE,
  verbose           = 2L
)
t2_fit_stationary_em <- Sys.time()
# --- CTMC-linked estimation ---
# The linked specification constrains lambda_g = -log(theta1)/Delta and
# lambda_d = -log(1-theta0)/Delta via the CTMC link, so durations and
# transitions are jointly identified from a single pair (theta1, theta0).
# ws$miscl_linked warm-starts from a previously converged linked result;
# ws$miscl_stationary_linked does the same for the stationary linked model.
# linked_init2 is intentionally data-driven (no warm-start).

message("=== Estimating non-stationary model (CTMC-linked) ===")
t1_fit_linked <- Sys.time()
fit_linked <- em_fit_tenure(
  df                = df_qlfs,
  params0           = ws$miscl_linked %||% custom_init,
  stationary        = FALSE,
  linked            = TRUE,
  discrete_timegap  = TRUE,
  verbose           = 2L
)
t2_fit_linked <- Sys.time()
#t1_fit_linked_data_driven_params <- Sys.time()
#fit_linked_data_driven_params <- em_fit_tenure(
#  df                = df_qlfs,
#  params0           = linked_init2,
#  stationary        = FALSE,
#  linked            = TRUE,
#  discrete_timegap  = TRUE,
#  verbose           = 2L
#)
#t2_fit_linked_data_driven_params <- Sys.time()

message("=== Estimating stationary model (CTMC-linked) ===")
t1_fit_stationary_linked <- Sys.time()
fit_stationary_linked <- em_fit_tenure(
  df                = df_qlfs,
  params0           = ws$miscl_stationary_linked %||% custom_init,
  stationary        = TRUE,
  linked            = TRUE,
  discrete_timegap  = TRUE,
  verbose           = 2L
)
t2_fit_stationary_linked <- Sys.time()
# --- Save results ---
# Each fit is saved as a self-contained .rds (includes gamma, history, params).
# A flat summary row is appended to run_summary.csv for quick comparison
# across runs without loading the full .rds files.
saveRDS(fit_em,
        file.path(results_dir, sprintf("fit_miscl_%s.rds", run_id)))
saveRDS(fit_stationary_em,
        file.path(results_dir, sprintf("fit_miscl_stationary_%s.rds", run_id)))
saveRDS(fit_linked,
        file.path(results_dir, sprintf("fit_miscl_linked_%s.rds", run_id)))
saveRDS(fit_stationary_linked,
        file.path(results_dir, sprintf("fit_miscl_stationary_linked_%s.rds", run_id)))
message(sprintf("Fits saved [run_id = %s]", run_id))

message(sprintf("\nLog-likelihood (free, non-stationary):   %.4f", fit_em$loglik))
message(sprintf("Log-likelihood (free, stationary):       %.4f", fit_stationary_em$loglik))
message(sprintf("Log-likelihood (linked, non-stationary): %.4f", fit_linked$loglik))
message(sprintf("Log-likelihood (linked, stationary):     %.4f", fit_stationary_linked$loglik))

fits <- list(
  miscl                    = fit_em,
  miscl_stationary         = fit_stationary_em,
  miscl_linked             = fit_linked,
  miscl_stationary_linked  = fit_stationary_linked
)

# Flat summary — one row per model per run.
# A single unified function covers base, rho, and eps models: model-specific
# fields (rho, eps, sigma2_g, sigma2_d) fall back to NA_real_ when absent.
.make_summary_row <- function(run_id, model, fit) {
  p <- fit$params
  data.table::data.table(
    run_id      = run_id,
    model       = model,
    converged   = fit$converged,
    iterations  = fit$iterations,
    loglik      = fit$loglik,
    alpha       = p$alpha,
    theta0      = p$theta0,
    theta1      = p$theta1,
    pi          = p$pi,
    rho         = p$rho      %||% NA_real_,
    eps         = p$eps      %||% NA_real_,
    sigma2_g    = p$sigma2_g %||% NA_real_,
    lambda_g    = p$lambda_g,
    lambda_d    = p$lambda_d,
    sigma2_d    = p$sigma2_d %||% NA_real_
  )
}

summary_rows <- data.table::rbindlist(
  lapply(names(fits), function(m) .make_summary_row(run_id, m, fits[[m]]))
)

data.table::fwrite(
  summary_rows,
  summary_csv,
  append    = file.exists(summary_csv),
  col.names = !file.exists(summary_csv)
)
message(sprintf("Summary appended: %s  [run_id = %s]", summary_csv, run_id))


# ##############################################################################
# RHO-AUGMENTED ESTIMATION (duration contamination model)
# ##############################################################################
# The rho model adds a duration contamination probability rho that separates
# tenure reporting error (correct state, wrong duration) from employment
# misclassification (wrong state, correct duration for that state).
#
# Created: 2026-04-29
# ##############################################################################

message("\n=== Estimating rho model (free, non-stationary) ===")
t1_fit_rho <- Sys.time()
fit_rho <- em_fit_tenure_rho(
  df         = df_qlfs,
  params0    = ws$rho %||% custom_init_rho,
  stationary = FALSE,
  linked     = FALSE,
  verbose    = 2L
)
t2_fit_rho <- Sys.time()
message("\n=== Estimating rho model (free, stationary) ===")
t1_fit_rho_stationary <- Sys.time()
fit_rho_stationary <- em_fit_tenure_rho(
  df         = df_qlfs,
  params0    = ws$rho_stationary %||% custom_init_rho,
  stationary = TRUE,
  linked     = FALSE,
  verbose    = 2L
)
t2_fit_rho_stationary <- Sys.time()

message("\n=== Estimating rho model (CTMC-linked, non-stationary) ===")
t1_fit_rho_linked <- Sys.time()
fit_rho_linked <- em_fit_tenure_rho(
  df         = df_qlfs,
  params0    = ws$rho_linked %||% custom_init_rho,
  stationary = FALSE,
  linked     = TRUE,
  verbose    = 2L
)
t2_fit_rho_linked <- Sys.time()

message("\n=== Estimating rho model (CTMC-linked, stationary) ===")
t1_fit_rho_stationary_linked <- Sys.time()
fit_rho_stationary_linked <- em_fit_tenure_rho(
  df         = df_qlfs,
  params0    = ws$rho_stationary_linked %||% custom_init_rho,
  stationary = TRUE,
  linked     = TRUE,
  verbose    = 2L
)
t2_fit_rho_stationary_linked <- Sys.time()

# --- Save rho model results ---
saveRDS(fit_rho,
        file.path(results_dir, sprintf("fit_rho_%s.rds", run_id)))
saveRDS(fit_rho_stationary,
        file.path(results_dir, sprintf("fit_rho_stationary_%s.rds", run_id)))
saveRDS(fit_rho_linked,
        file.path(results_dir, sprintf("fit_rho_linked_%s.rds", run_id)))
saveRDS(fit_rho_stationary_linked,
        file.path(results_dir, sprintf("fit_rho_stationary_linked_%s.rds", run_id)))

message(sprintf("\nLog-likelihood (rho, free, non-stationary): %.4f", fit_rho$loglik))
message(sprintf("Log-likelihood (rho, free, stationary):     %.4f", fit_rho_stationary$loglik))
message(sprintf("Log-likelihood (rho, linked, non-stat):     %.4f", fit_rho_linked$loglik))
message(sprintf("Log-likelihood (rho, linked, stationary):   %.4f", fit_rho_stationary_linked$loglik))

# LR test: rho model vs base model (1 df for rho)
lr_stat <- 2 * (fit_rho$loglik - fit_em$loglik)
lr_pval <- pchisq(lr_stat, df = 1, lower.tail = FALSE)
message(sprintf("\nLR test (rho vs base): stat = %.4f, p = %.6f", lr_stat, lr_pval))

rho_fits <- list(
  rho_free              = fit_rho,
  rho_stationary        = fit_rho_stationary,
  rho_linked            = fit_rho_linked,
  rho_stationary_linked = fit_rho_stationary_linked
)

rho_summary <- data.table::rbindlist(
  lapply(names(rho_fits), function(m) .make_summary_row(run_id, m, rho_fits[[m]]))
)

data.table::fwrite(
  rho_summary,
  summary_csv,
  append    = TRUE,
  col.names = FALSE
)
message(sprintf("Rho summary appended: %s  [run_id = %s]", summary_csv, run_id))

# --- Runtime summary ---
model_runtimes <- list(
  fit_em                     = t2_fit_em - t1_fit_em,
  fit_stationary_em          = t2_fit_stationary_em - t1_fit_stationary_em,
  fit_linked                 = t2_fit_linked - t1_fit_linked,
  fit_stationary_linked      = t2_fit_stationary_linked - t1_fit_stationary_linked,
  fit_rho                    = t2_fit_rho - t1_fit_rho,
  fit_rho_stationary         = t2_fit_rho_stationary - t1_fit_rho_stationary,
  fit_rho_linked             = t2_fit_rho_linked - t1_fit_rho_linked,
  fit_rho_stationary_linked  = t2_fit_rho_stationary_linked - t1_fit_rho_stationary_linked
)

message("\n=== Model runtimes (t2 - t1) ===")
print(model_runtimes)

total_runtime <- sum(vapply(model_runtimes, function(x) {
  as.numeric(x, units = "secs")
}, numeric(1L)))
message(sprintf("Total model runtime: %.2f seconds", total_runtime))


# ##############################################################################
# EPSILON-AUGMENTED ESTIMATION (Spec I: point-mass + Exp contamination)
# ##############################################################################
# The eps model replaces the Gaussian-EMG measurement model with a structurally
# correct point-mass + Exponential contamination mixture, and adds a spell-pair
# joint emission so that tenure flanks identify misclassification (Story A).
# Companion spec: documents/EM tenure epsilon.tex
# Design rationale: EM-tenure/feedback/2026-04-30-epsilon-spec-design.md
# Discussed in 'EM tenure epsilon.tex'
#
# Created: 2026-04-30
# ##############################################################################

message("\n=== Estimating eps model (Spec I) ===")

# Initial parameters for eps model — separate objects for free vs linked
# variants so the fallback parameter structure matches each model's spec.
custom_init_eps        <- init_params_eps(df_qlfs, linked = FALSE)
custom_init_eps_linked <- init_params_eps(df_qlfs, linked = TRUE)

message("=== Estimating eps model (free, non-stationary) ===")
t1_fit_eps <- Sys.time()
fit_eps <- em_fit_tenure_eps(
  df         = df_qlfs,
  params0    = ws$eps %||% custom_init_eps,
  stationary = FALSE,
  linked     = FALSE,
  verbose    = 2L
)
t2_fit_eps <- Sys.time()

message("=== Estimating eps model (free, stationary) ===")
t1_fit_eps_stationary <- Sys.time()
fit_eps_stationary <- em_fit_tenure_eps(
  df         = df_qlfs,
  params0    = ws$eps_stationary %||% custom_init_eps,
  stationary = TRUE,
  linked     = FALSE,
  verbose    = 2L
)
t2_fit_eps_stationary <- Sys.time()

message("=== Estimating eps model (CTMC-linked, non-stationary) ===")
t1_fit_eps_linked <- Sys.time()
fit_eps_linked <- em_fit_tenure_eps(
  df         = df_qlfs,
  params0    = ws$eps_linked %||% custom_init_eps_linked,
  stationary = FALSE,
  linked     = TRUE,
  verbose    = 2L
)
t2_fit_eps_linked <- Sys.time()

message("=== Estimating eps model (CTMC-linked, stationary) ===")
t1_fit_eps_stationary_linked <- Sys.time()
fit_eps_stationary_linked <- em_fit_tenure_eps(
  df         = df_qlfs,
  params0    = ws$eps_stationary_linked %||% custom_init_eps_linked,
  stationary = TRUE,
  linked     = TRUE,
  verbose    = 2L
)
t2_fit_eps_stationary_linked <- Sys.time()

# --- Save eps model results ---
saveRDS(fit_eps,
        file.path(results_dir, sprintf("fit_eps_%s.rds", run_id)))
saveRDS(fit_eps_stationary,
        file.path(results_dir, sprintf("fit_eps_stationary_%s.rds", run_id)))
saveRDS(fit_eps_linked,
        file.path(results_dir, sprintf("fit_eps_linked_%s.rds", run_id)))
saveRDS(fit_eps_stationary_linked,
        file.path(results_dir, sprintf("fit_eps_stationary_linked_%s.rds", run_id)))

message(sprintf("\nLog-likelihood (eps, free, non-stationary): %.4f", fit_eps$loglik))
message(sprintf("Log-likelihood (eps, free, stationary):     %.4f", fit_eps_stationary$loglik))
message(sprintf("Log-likelihood (eps, linked, non-stat):     %.4f", fit_eps_linked$loglik))
message(sprintf("Log-likelihood (eps, linked, stationary):   %.4f", fit_eps_stationary_linked$loglik))

# Delta-LL vs base (non-nested in distribution family but informative)
lr_eps_vs_base <- 2 * (fit_eps$loglik - fit_em$loglik)
message(sprintf("\nDelta-LL (eps vs base, non-stat): %.4f", lr_eps_vs_base))

eps_fits <- list(
  eps_free              = fit_eps,
  eps_stationary        = fit_eps_stationary,
  eps_linked            = fit_eps_linked,
  eps_stationary_linked = fit_eps_stationary_linked
)

eps_summary <- data.table::rbindlist(
  lapply(names(eps_fits), function(m) .make_summary_row(run_id, m, eps_fits[[m]]))
)

data.table::fwrite(
  eps_summary,
  summary_csv,
  append    = TRUE,
  col.names = FALSE
)
message(sprintf("Eps summary appended: %s  [run_id = %s]", summary_csv, run_id))

# Eps model runtimes
eps_runtimes <- list(
  fit_eps                   = t2_fit_eps - t1_fit_eps,
  fit_eps_stationary        = t2_fit_eps_stationary - t1_fit_eps_stationary,
  fit_eps_linked            = t2_fit_eps_linked - t1_fit_eps_linked,
  fit_eps_stationary_linked = t2_fit_eps_stationary_linked - t1_fit_eps_stationary_linked
)
message(sprintf("  fit_eps:                    %.1f s", as.numeric(eps_runtimes$fit_eps,                   units = "secs")))
message(sprintf("  fit_eps_stationary:         %.1f s", as.numeric(eps_runtimes$fit_eps_stationary,        units = "secs")))
message(sprintf("  fit_eps_linked:             %.1f s", as.numeric(eps_runtimes$fit_eps_linked,            units = "secs")))
message(sprintf("  fit_eps_stationary_linked:  %.1f s", as.numeric(eps_runtimes$fit_eps_stationary_linked, units = "secs")))

total_runtime_eps <- sum(vapply(eps_runtimes, function(x) as.numeric(x, units = "secs"), numeric(1L)))
message(sprintf("Total eps model runtime: %.2f seconds", total_runtime_eps))

# ##############################################################################
# LIKELIHOOD RATIO TESTS: fit_eps_stationary_linked vs other eps models
# ##############################################################################
# fit_eps_stationary_linked is the preferred model (most parsimonious: 4 free
# parameters -- theta0, theta1, pi, eps).  The three alternative models are
# strictly less restricted, so all three tests are proper nested LR tests.
#
# Parameter counts:
#   fit_eps                     : 7 (alpha, theta0, theta1, pi, eps, lambda_g, lambda_d)
#   fit_eps_stationary          : 6 (alpha = theta1/(theta0+theta1) imposed; -1 df)
#   fit_eps_linked              : 5 (CTMC link lambda_g, lambda_d imposed; -2 df)
#   fit_eps_stationary_linked   : 4 (both restrictions; -3 df vs free model)
# ##############################################################################

message("\n=== LR tests: fit_eps_stationary_linked (H0) vs less-restricted eps models ===")

# --- H0: stationary_linked  vs  H1: linked (non-stationary) ---
# Tests the stationarity restriction alpha = theta1 / (theta0 + theta1); 1 df.
lr_stat_vs_linked <- 2 * (fit_eps_linked$loglik - fit_eps_stationary_linked$loglik)
lr_pval_vs_linked <- pchisq(lr_stat_vs_linked, df = 1L, lower.tail = FALSE)
message(sprintf(
  "  stationary_linked vs linked (non-stationary) | df=1 | LR=%.4f | p=%.6f",
  lr_stat_vs_linked, lr_pval_vs_linked
))

# --- H0: stationary_linked  vs  H1: stationary (free lambda) ---
# Tests the CTMC-link restrictions lambda_g = -log(theta1)/Delta and
# lambda_d = -log(1-theta0)/Delta; 2 df.
lr_stat_vs_stationary <- 2 * (fit_eps_stationary$loglik - fit_eps_stationary_linked$loglik)
lr_pval_vs_stationary <- pchisq(lr_stat_vs_stationary, df = 2L, lower.tail = FALSE)
message(sprintf(
  "  stationary_linked vs stationary (free lambda) | df=2 | LR=%.4f | p=%.6f",
  lr_stat_vs_stationary, lr_pval_vs_stationary
))

# --- H0: stationary_linked  vs  H1: free (non-stationary, free lambda) ---
# Tests both the stationarity restriction and the CTMC-link jointly; 3 df.
lr_stat_vs_free <- 2 * (fit_eps$loglik - fit_eps_stationary_linked$loglik)
lr_pval_vs_free <- pchisq(lr_stat_vs_free, df = 3L, lower.tail = FALSE)
message(sprintf(
  "  stationary_linked vs free (non-stat, free lambda) | df=3 | LR=%.4f | p=%.6f",
  lr_stat_vs_free, lr_pval_vs_free
))

message(sprintf(
  "\n  Preferred model: fit_eps_stationary_linked  (loglik = %.4f)",
  fit_eps_stationary_linked$loglik
))




