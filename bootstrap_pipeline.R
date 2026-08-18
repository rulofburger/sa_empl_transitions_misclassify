# ==============================================================================
# Bootstrap inference pipeline for all EM models
# Created: 2026-05-06
#
# Computes bootstrap standard errors for all 22 EM model configurations:
#   - 6  baseline models  (EM-baseline/output/results/fit_*.rds)
#   - 16 extension models (EM-baseline-ext/output/results/fit_*.rds)
#
# For each model, runs B nonparametric bootstrap replicates, warm-started from
# the point estimate (no multi-start). Uses parallel::mclapply for within-model
# parallelism.  Results are saved per model so individual models can be
# re-run without repeating the full pipeline.
#
# Bootstrap design:
#   - Nonparametric: individuals are resampled with replacement.
#   - Weights: preserved as columns of the resampled data frame.
#   - B = 200 (labeled in filenames; increase B and re-run for journal submission).
#   - Seeds: deterministic per-rep seed vectors drawn from master seed.
#   - Quality flags: reps that fail to converge or have anomalously low loglik
#     (> 50 nats below point estimate) are stored but excluded from SE computation.
#
# Usage (from project root):
#   Rscript bootstrap_pipeline.R
#
# Prerequisites:
#   estimate_baseline_pipeline.R and estimate_extensions_pipeline.R must have
#   been run first to produce the point-estimate .rds files.
#
# Outputs:
#   EM-baseline/output/results/bootstrap/boot_{label}_B200.rds
#   EM-baseline-ext/output/results/bootstrap/boot_{label}_B200.rds
#
#   Each .rds file contains a list with elements:
#     $boot_results    — B-element list of replicate outputs (params, implied, flag)
#     $summary         — data frame from summarise_bootstrap()
#     $ame_summary     — data frame from summarise_bootstrap_ame() (cov models only)
#     $point_params    — point-estimate parameter list
#     $point_implied   — point-estimate implied quantities list
#     $point_loglik    — point-estimate log-likelihood
#     $B               — number of bootstrap reps attempted
#     $n_ok            — number of "ok" reps included in SEs
#     $label           — model label
#     $run_ts          — timestamp of this bootstrap run
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------

library(here)
library(dplyr)
library(parallel)

# Source all EM modules (baseline + extensions; includes bootstrap_utils and
# bootstrap_utils_ext via the updated source_all.R files)
source(here::here("EM-baseline",     "R", "source_all.R"))
source(here::here("EM-baseline-ext", "R", "source_all.R"))

# Bootstrap hyperparameters
# *** IMPORTANT: If you change B or MASTER_SEED, delete all existing bootstrap
# *** files first — the skip-if-exists logic will silently reuse cached results.
# *** Also update B in build_tables.R if you are not using auto-detection.
B           <- 200L             # number of bootstrap replicates (label in filenames)
MASTER_SEED <- 42L              # reproducibility seed
N_CORES     <- max(1L, parallel::detectCores() - 1L)

# Estimation hyperparameters (must match original pipelines)
THETA_CAP <- 0.999
PI_CAP    <- 0.49

# Output directories
boot_baseline_dir <- here::here("EM-baseline",     "output", "results", "bootstrap")
boot_ext_dir      <- here::here("EM-baseline-ext", "output", "results", "bootstrap")
dir.create(boot_baseline_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(boot_ext_dir,      recursive = TRUE, showWarnings = FALSE)

run_ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

cat(sprintf("Bootstrap pipeline — B=%d | cores=%d | seed=%d\n",
            B, N_CORES, MASTER_SEED))
cat(sprintf("Run timestamp: %s\n\n", run_ts))

# Draw a B × 22 matrix of deterministic per-rep seeds (rows = reps, cols = models)
# so each model gets its own independent seed vector regardless of parallelism order.
set.seed(MASTER_SEED)
N_MODELS   <- 22L
seed_matrix <- matrix(
  sample.int(.Machine$integer.max, B * N_MODELS, replace = FALSE),
  nrow = B, ncol = N_MODELS
)

# ------------------------------------------------------------------------------
# 1. Load and validate data
# ------------------------------------------------------------------------------

source(here::here("scripts", "ingest_data_3waves_SA.R"))  # loads df_qlfs

sector_source_path <- here::here("data", "raw", "QLFSmerged_mapped.rds")
if (!file.exists(sector_source_path))
  stop("Missing sector source: data/raw/QLFSmerged_mapped.rds")
sector_source <- readRDS(sector_source_path)
df_qlfs <- attach_transition_informal_sector(df_qlfs, sector_source)
rm(sector_source)

# Guard: factor check
for (y_col in c("y1", "y2", "y3")) {
  if (is.factor(df_qlfs[[y_col]]))
    stop(sprintf("Column '%s' is a factor — check ingest script.", y_col))
}

# Build df_baseline (same filter as estimate_baseline_pipeline.R)
df_baseline <- df_qlfs |>
  filter(!is.na(y1), !is.na(y2), !is.na(y3), !is.na(weight)) |>
  mutate(
    y1     = as.integer(y1),
    y2     = as.integer(y2),
    y3     = as.integer(y3),
    weight = as.numeric(weight)
  ) |>
  select(y1, y2, y3, weight)

# Build df_ext (same filter as estimate_extensions_pipeline.R)
keep <- !is.na(df_qlfs$y1) & !is.na(df_qlfs$y2) & !is.na(df_qlfs$y3) &
        !is.na(df_qlfs$weight) & df_qlfs$weight > 0 &
        !is.na(df_qlfs$age1)   & !is.na(df_qlfs$age2)   & !is.na(df_qlfs$age3) &
        !is.na(df_qlfs$educ1)  & !is.na(df_qlfs$educ2)  & !is.na(df_qlfs$educ3)
df_ext        <- df_qlfs[keep, , drop = FALSE]
df_ext$y1     <- as.integer(df_ext$y1)
df_ext$y2     <- as.integer(df_ext$y2)
df_ext$y3     <- as.integer(df_ext$y3)
df_ext$weight <- as.numeric(df_ext$weight)
for (nm in c("contracttype1", "contracttype2"))
  df_ext[[nm]] <- ifelse(is.na(df_ext[[nm]]), 0L, as.integer(df_ext[[nm]]))
df_ext <- as.data.frame(df_ext)

# Build covariate matrices and inconsistency matrix
cv_set1 <- prepare_covariate_matrix(df_ext, covariate_set = 1L)
cv_set2 <- prepare_covariate_matrix(df_ext, covariate_set = 2L)
cv_set3 <- prepare_covariate_matrix(df_ext, covariate_set = 3L)
X1 <- cv_set1$X
X2 <- cv_set2$X
X3 <- cv_set3$X_transition

df_incons  <- compute_inconsistencies(df_ext)
inc_mat    <- as.matrix(df_incons[, c("Y_age_1", "Y_age_2", "Y_age_3",
                                       "Y_edu_1", "Y_edu_2", "Y_edu_3")])
rm(df_qlfs, df_incons)

cat(sprintf("Data loaded: baseline N=%d | ext N=%d\n\n",
            nrow(df_baseline), nrow(df_ext)))

# ------------------------------------------------------------------------------
# 2. Helper: run bootstrap for one model and save
# ------------------------------------------------------------------------------

.run_and_save <- function(label, boot_dir, seeds,
                          fit_fn, implied_fn, save_ame = FALSE,
                          point_params, point_implied, point_loglik) {
  out_path <- file.path(boot_dir, sprintf("boot_%s_B%d.rds", label, B))

  if (file.exists(out_path)) {
    # NOTE: Skip is based on filename only (includes B count). If MASTER_SEED,
    # point estimates, or covariate matrices change without changing B, cached
    # bootstrap results will be silently reused. Delete bootstrap/ dir to force
    # a full rerun.
    cat(sprintf("  [SKIP] %s — file exists\n", label))
    return(invisible(NULL))
  }

  cat(sprintf("  [START] %s ...\n", label))
  t0 <- proc.time()

  boot_results <- parallel::mclapply(
    seq_len(B),
    function(b) fit_fn(seeds[b]),
    mc.cores = N_CORES
  )

  flags <- vapply(boot_results, function(r) r$flag, character(1L))
  n_ok  <- sum(flags == "ok")

  summary_df  <- summarise_bootstrap(boot_results, point_params,
                                      point_implied, point_loglik)
  ame_summary <- NULL
  if (save_ame && !is.null(point_implied$ame_entry)) {
    ame_summary <- summarise_bootstrap_ame(
      boot_results,
      point_ame_entry = point_implied$ame_entry,
      point_ame_exit  = point_implied$ame_exit
    )
  }

  elapsed <- (proc.time() - t0)[["elapsed"]]
  cat(sprintf("  [DONE]  %s — n_ok=%d/%d | %.1fs\n", label, n_ok, B, elapsed))

  saveRDS(list(
    boot_results  = boot_results,
    summary       = summary_df,
    ame_summary   = ame_summary,
    point_params  = point_params,
    point_implied = point_implied,
    point_loglik  = point_loglik,
    B             = B,
    n_ok          = n_ok,
    label         = label,
    run_ts        = run_ts
  ), out_path)
}

# ------------------------------------------------------------------------------
# 3. Bootstrap baseline models (6 models)
# ------------------------------------------------------------------------------

cat("========== BASELINE MODELS ==========\n")

baseline_configs <- list(
  list(label = "none_stat",  model_type = "none",       stationary = TRUE),
  list(label = "none_free",  model_type = "none",       stationary = FALSE),
  list(label = "sym_stat",   model_type = "symmetric",  stationary = TRUE),
  list(label = "sym_free",   model_type = "symmetric",  stationary = FALSE),
  list(label = "asym_stat",  model_type = "asymmetric", stationary = TRUE),
  list(label = "asym_free",  model_type = "asymmetric", stationary = FALSE)
)

for (k in seq_along(baseline_configs)) {
  cfg        <- baseline_configs[[k]]
  rds_path   <- here::here("EM-baseline", "output", "results",
                            sprintf("fit_%s.rds", cfg$label))
  if (!file.exists(rds_path)) {
    warning(sprintf("Point estimate missing: %s — skipping bootstrap.", rds_path))
    next
  }
  fit_pt       <- readRDS(rds_path)
  stopifnot(
    !is.null(fit_pt$params),
    is.list(fit_pt$params),
    is.numeric(fit_pt$loglik), length(fit_pt$loglik) == 1L,
    is.finite(fit_pt$loglik)
  )
  point_params  <- fit_pt$params
  point_loglik  <- fit_pt$loglik
  point_implied <- implied_baseline(point_params, cfg$model_type)

  model_type  <- cfg$model_type
  stationary  <- cfg$stationary
  seeds_k     <- seed_matrix[, k]

  .run_and_save(
    label         = cfg$label,
    boot_dir      = boot_baseline_dir,
    seeds         = seeds_k,
    fit_fn        = function(seed) {
      bootstrap_one_baseline(df_baseline, seed, model_type, stationary,
                             point_params, point_loglik, THETA_CAP, PI_CAP)
    },
    implied_fn    = NULL,
    save_ame      = FALSE,
    point_params  = point_params,
    point_implied = point_implied,
    point_loglik  = point_loglik
  )
}

# ------------------------------------------------------------------------------
# 4. Bootstrap extension models (16 models)
# ------------------------------------------------------------------------------

cat("\n========== EXTENSION MODELS ==========\n")

# Helper: build covariate model bootstrap runner
.make_cov_fn <- function(df, X, model_type, stationary,
                         point_params, point_loglik, pi_cap) {
  function(seed) {
    bootstrap_one_covariates(df, X, seed, model_type, stationary,
                             point_params, point_loglik, pi_cap)
  }
}

# Helper: build FMM model bootstrap runner
.make_fmm_fn <- function(df, model_type, stationary,
                         point_params, point_loglik, theta_cap, pi_cap) {
  function(seed) {
    bootstrap_one_fmm(df, seed, model_type, stationary,
                      point_params, point_loglik, theta_cap, pi_cap)
  }
}

# Helper: build inconsistency model bootstrap runner
.make_incons_fn <- function(df, incons_m, stationary,
                            point_params, point_loglik, theta_cap) {
  function(seed) {
    bootstrap_one_inconsistency(df, incons_m, seed, stationary,
                                point_params, point_loglik, theta_cap)
  }
}

ext_configs <- list(
  # Covariate Set 1
  list(label="cov_s1_sym_stat", family="cov", X=X1, model_type="symmetric", stationary=TRUE),
  list(label="cov_s1_sym_free", family="cov", X=X1, model_type="symmetric", stationary=FALSE),
  list(label="cov_s1_non_stat", family="cov", X=X1, model_type="none",      stationary=TRUE),
  list(label="cov_s1_non_free", family="cov", X=X1, model_type="none",      stationary=FALSE),
  # Covariate Set 2
  list(label="cov_s2_sym_stat", family="cov", X=X2, model_type="symmetric", stationary=TRUE),
  list(label="cov_s2_sym_free", family="cov", X=X2, model_type="symmetric", stationary=FALSE),
  list(label="cov_s2_non_stat", family="cov", X=X2, model_type="none",      stationary=TRUE),
  list(label="cov_s2_non_free", family="cov", X=X2, model_type="none",      stationary=FALSE),
  # Covariate Set 3
  list(label="cov_s3_sym_free", family="cov", X=X3, model_type="symmetric", stationary=FALSE),
  list(label="cov_s3_non_free", family="cov", X=X3, model_type="none",      stationary=FALSE),
  # FMM
  list(label="fmm_sym_stat", family="fmm", model_type="symmetric", stationary=TRUE),
  list(label="fmm_sym_free", family="fmm", model_type="symmetric", stationary=FALSE),
  list(label="fmm_non_stat", family="fmm", model_type="none",      stationary=TRUE),
  list(label="fmm_non_free", family="fmm", model_type="none",      stationary=FALSE),
  # Inconsistency
  list(label="incons_sym_stat", family="incons", stationary=TRUE),
  list(label="incons_sym_free", family="incons", stationary=FALSE)
)

# Release covariate matrices before FMM/inconsistency models reduce per-fork
# memory pressure (macOS fork() dirty-page copies).
rm(X1, X2, X3)
gc(invisible = TRUE)

for (k in seq_along(ext_configs)) {
  cfg      <- ext_configs[[k]]
  rds_path <- here::here("EM-baseline-ext", "output", "results",
                          sprintf("fit_%s.rds", cfg$label))
  if (!file.exists(rds_path)) {
    warning(sprintf("Point estimate missing: %s — skipping.", rds_path))
    next
  }

  fit_pt      <- readRDS(rds_path)
  stopifnot(
    !is.null(fit_pt$params),
    is.list(fit_pt$params),
    is.numeric(fit_pt$loglik), length(fit_pt$loglik) == 1L,
    is.finite(fit_pt$loglik)
  )
  point_params <- fit_pt$params
  point_loglik <- fit_pt$loglik
  model_col    <- length(baseline_configs) + k  # column index in seed_matrix
  seeds_k      <- seed_matrix[, model_col]

  if (cfg$family == "cov") {
    coef_names <- colnames(.as_transition_design(cfg$X)$X12)
    names(point_params$beta0) <- coef_names
    names(point_params$beta1) <- coef_names
    point_implied <- implied_covariates(point_params, cfg$X, cfg$model_type,
                                        df = df_ext, gamma = fit_pt$gamma)
    fit_fn <- .make_cov_fn(df_ext, cfg$X, cfg$model_type, cfg$stationary,
                            point_params, point_loglik, PI_CAP)
    save_ame <- TRUE

  } else if (cfg$family == "fmm") {
    point_implied <- implied_fmm(point_params, cfg$model_type)
    fit_fn <- .make_fmm_fn(df_ext, cfg$model_type, cfg$stationary,
                            point_params, point_loglik, THETA_CAP, PI_CAP)
    save_ame <- FALSE

  } else {  # inconsistency
    point_implied <- implied_inconsistency(point_params, inc_mat)
    fit_fn <- .make_incons_fn(df_ext, inc_mat, cfg$stationary,
                               point_params, point_loglik, THETA_CAP)
    save_ame <- FALSE
  }

  .run_and_save(
    label         = cfg$label,
    boot_dir      = boot_ext_dir,
    seeds         = seeds_k,
    fit_fn        = fit_fn,
    implied_fn    = NULL,
    save_ame      = save_ame,
    point_params  = point_params,
    point_implied = point_implied,
    point_loglik  = point_loglik
  )
}

# ------------------------------------------------------------------------------
# 5. Convergence summary
# ------------------------------------------------------------------------------

cat("\n========== CONVERGENCE SUMMARY ==========\n")

all_boot_files <- c(
  list.files(boot_baseline_dir, pattern = sprintf("_B%d\\.rds$", B), full.names = TRUE),
  list.files(boot_ext_dir,      pattern = sprintf("_B%d\\.rds$", B), full.names = TRUE)
)

if (length(all_boot_files) > 0L) {
  cat(sprintf("%-28s  %6s  %6s  %8s\n", "Model", "n_ok", "B", "pct_ok"))
  cat(strrep("-", 56), "\n")
  for (f in sort(all_boot_files)) {
    obj <- readRDS(f)
    cat(sprintf("%-28s  %6d  %6d  %7.1f%%\n",
                obj$label, obj$n_ok, obj$B, 100 * obj$n_ok / obj$B))
  }
}

cat(sprintf("\nBootstrap pipeline complete. Results in:\n  %s\n  %s\n",
            boot_baseline_dir, boot_ext_dir))
