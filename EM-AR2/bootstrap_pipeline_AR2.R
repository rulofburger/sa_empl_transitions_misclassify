# ==============================================================================
# EM-AR2: Bootstrap inference pipeline
# Created: 2026-05-07
#
# Computes bootstrap standard errors for all 3 AR(2) EM model configurations:
#   - nopi  : no misclassification (pi = 0)
#   - sym   : symmetric misclassification (pi)
#   - asym  : asymmetric misclassification (pi0, pi1)
#
# Bootstrap design:
#   - Nonparametric: individuals resampled with replacement (all 4 waves +
#     weight preserved). Preserves intra-individual correlation without
#     assuming a parametric model.
#   - B = 200 reps (labeled in filenames). Sufficient for 2 d.p. SE estimation
#     (Efron & Tibshirani 1993, rule of thumb B ≈ 200 for SEs).
#   - Warm-start from point estimate: EM initialized at θ̂ (MLE from full data)
#     rather than random starts, reducing computation and guarding saddle points.
#   - Seeds: deterministic B × 3 matrix from MASTER_SEED.
#   - Quality flags: non-convergence and reps with loglik more than 50 nats
#     below the point estimate are stored but excluded from SE computation.
#   - Parallelisation: parallel::mclapply (fork-based, macOS/Linux).
#     NOT available on Windows — falls back to single-core automatically.
#
# Point estimates: loaded from the most recent timestamped .rds file per model
#   in EM-AR2/output/results/ via .find_latest_fit().
#
# Outputs (per model):
#   EM-AR2/output/results/bootstrap/boot_{label}_B200.rds
#   Each .rds contains: boot_results, summary, point_params, point_implied,
#   point_loglik, B, n_ok, label, run_ts.
#
# Prerequisites:
#   estimate_pipeline.R must have been run first.
#
# Usage (from project root):
#   source("EM-AR2/bootstrap_pipeline_AR2.R")
# ==============================================================================

# Verify working directory -------------------------------------------------------
if (!file.exists("EM-AR2/R/source_all.R")) {
  stop(
    "bootstrap_pipeline_AR2.R must be sourced from the project root. ",
    "Expected to find 'EM-AR2/R/source_all.R' relative to cwd."
  )
}

# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------

library(tidyverse)   # required by ingest_data_4waves_SA.R
library(parallel)

# Source EM-baseline for bootstrap_resample() and summarise_bootstrap()
source("EM-baseline/R/source_all.R")

# Source EM-AR2 for em_fit_ar2(), implied_ar2(), .find_latest_fit(),
# bootstrap_one_ar2()
source("EM-AR2/R/source_all.R")

# Bootstrap hyperparameters
# *** IMPORTANT: changing B or MASTER_SEED invalidates cached results.
# *** Delete existing boot_*_B200.rds files before re-running.
B           <- 200L
MASTER_SEED <- 42L
N_CORES     <- max(1L, parallel::detectCores() - 1L)
# Fork-based mclapply is not available on Windows; fall back to single-core.
if (.Platform$OS.type == "windows") {
  N_CORES <- 1L
  warning("bootstrap_pipeline_AR2.R: mclapply unavailable on Windows; running single-core.")
}
# Upper bound on misclassification probability for numerical stability.
# Values at or above 0.5 cause identification failure in the EM algorithm.
PI_CAP      <- 0.49

# Directories
results_dir   <- "EM-AR2/output/results"
boot_dir      <- "EM-AR2/output/results/bootstrap"
dir.create(boot_dir, recursive = TRUE, showWarnings = FALSE)

run_ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

cat(sprintf("EM-AR2 Bootstrap pipeline — B=%d | cores=%d | seed=%d\n",
            B, N_CORES, MASTER_SEED))
cat(sprintf("Run timestamp: %s\n\n", run_ts))

# Deterministic per-rep seeds: B × 3 matrix (rows = reps, cols = models)
set.seed(MASTER_SEED)
N_MODELS    <- 3L
seed_matrix <- matrix(
  sample.int(.Machine$integer.max, B * N_MODELS, replace = FALSE),
  nrow = B, ncol = N_MODELS
)

# ------------------------------------------------------------------------------
# 1. Load and validate data
# ------------------------------------------------------------------------------

if (!file.exists("scripts/ingest_data_4waves_SA.R")) {
  stop("Missing ingestion script: 'scripts/ingest_data_4waves_SA.R'.")
}
if (!file.exists("data/raw/df_qlfs_A.rds")) {
  stop("Missing input data: 'data/raw/df_qlfs_A.rds'.")
}

source("scripts/ingest_data_4waves_SA.R")  # loads df_qlfs
if (!exists("df_qlfs") || !is.data.frame(df_qlfs))
  stop("ingest_data_4waves_SA.R did not produce a 'df_qlfs' data frame. ",
       "Check the ingestion script for errors.")

df_raw <- df_qlfs

if (nrow(df_raw) == 0L)
  stop("The complete four-wave estimation sample is empty. Check the input panel.")
if (anyNA(df_raw$weight))
  stop(sprintf("weight column contains %d NA(s). Check ingestion script.",
               sum(is.na(df_raw$weight))))
if (any(df_raw$weight < 0))
  stop("weight column contains negative values.")

n_obs        <- nrow(df_raw)
weight_total <- sum(df_raw$weight)
if (weight_total <= 0)
  stop(sprintf("Total weight is non-positive (%g). Cannot normalize.", weight_total))

df_4w        <- df_raw
df_4w$weight <- n_obs * df_4w$weight / weight_total

required_cols <- c("y1", "y2", "y3", "y4", "weight")
missing <- setdiff(required_cols, names(df_4w))
if (length(missing) > 0L)
  stop("Data missing required columns: ", paste(missing, collapse = ", "))

# Cast y columns to integer (guard against factor output from ingestion)
for (col in c("y1", "y2", "y3", "y4")) {
  if (is.factor(df_4w[[col]]))
    df_4w[[col]] <- as.integer(as.character(df_4w[[col]]))
}
rm(col, df_qlfs, df_raw)

cat(sprintf("Data loaded: N = %d\n\n", n_obs))

# ------------------------------------------------------------------------------
# 2. Model configurations
# ------------------------------------------------------------------------------

configs <- list(
  list(key = "nopi", model_type = "none"),
  list(key = "sym",  model_type = "symmetric"),
  list(key = "asym", model_type = "asymmetric")
)

stopifnot(length(configs) == N_MODELS, ncol(seed_matrix) == N_MODELS)

# Accumulate convergence log in memory to avoid re-reading .rds files in summary.
convergence_log <- vector("list", N_MODELS)

# ------------------------------------------------------------------------------
# 3. Bootstrap each model
# ------------------------------------------------------------------------------

for (k in seq_along(configs)) {
  cfg        <- configs[[k]]
  key        <- cfg$key
  model_type <- cfg$model_type
  out_path   <- file.path(boot_dir, sprintf("boot_%s_B%d.rds", key, B))
  seeds_k    <- seed_matrix[, k]

  # Skip if output already exists
  if (file.exists(out_path)) {
    cat(sprintf("  [SKIP] %s — boot file already exists: %s\n", key, out_path))
    next
  }

  # Load latest point-estimate fit
  fit_path <- tryCatch(
    .find_latest_fit(key, results_dir),
    error = function(e) {
      message(sprintf("  [SKIP] %s — could not find fit: %s", key,
                      conditionMessage(e)))
      NULL
    }
  )
  if (is.null(fit_path)) next

  fit_pt <- readRDS(fit_path)
  if (!is.list(fit_pt) || is.null(fit_pt$params) || is.null(fit_pt$loglik))
    stop(sprintf("Fit file '%s' is missing $params or $loglik.", fit_path))

  point_params  <- fit_pt$params
  point_loglik  <- fit_pt$loglik
  point_implied <- tryCatch(
    implied_ar2(point_params, model_type),
    error = function(e) NULL
  )

  cat(sprintf("  [START] %s (model_type=%s, LL=%.2f) ...\n",
              key, model_type, point_loglik))
  t0 <- proc.time()

  boot_results <- parallel::mclapply(
    seq_len(B),
    function(b) {
      bootstrap_one_ar2(
        df           = df_4w,
        seed         = seeds_k[[b]],
        model_type   = model_type,
        params_start = point_params,
        point_loglik = point_loglik,
        pi_cap       = PI_CAP
      )
    },
    mc.cores = N_CORES
  )

  flags <- vapply(boot_results, function(r) r$flag, character(1L))
  n_ok  <- sum(flags == "ok")

  summary_df <- summarise_bootstrap(
    boot_results  = boot_results,
    point_params  = point_params,
    point_implied = point_implied,
    point_loglik  = point_loglik
  )

  elapsed <- (proc.time() - t0)[["elapsed"]]
  cat(sprintf("  [DONE]  %s — n_ok=%d/%d | %.1fs\n", key, n_ok, B, elapsed))

  saveRDS(
    list(
      boot_results  = boot_results,
      summary       = summary_df,
      point_params  = point_params,
      point_implied = point_implied,
      point_loglik  = point_loglik,
      B             = B,
      n_obs         = n_obs,
      n_ok          = n_ok,
      label         = key,
      run_ts        = run_ts
    ),
    out_path
  )
  convergence_log[[k]] <- list(label = key, n_ok = n_ok, B = B)
  cat(sprintf("  Saved: %s\n\n", out_path))
}

# ------------------------------------------------------------------------------
# 4. Convergence summary (from in-memory log; no re-read of .rds files)
# ------------------------------------------------------------------------------

cat("\n========== CONVERGENCE SUMMARY ==========\n")
for (entry in convergence_log) {
  if (is.null(entry)) next
  cat(sprintf("  %-12s  n_ok=%d/%d  (%.1f%% ok)\n",
              entry$label, entry$n_ok, entry$B,
              100 * entry$n_ok / entry$B))
}

cat(sprintf("\nBootstrap pipeline complete. Results in: %s\n", boot_dir))
