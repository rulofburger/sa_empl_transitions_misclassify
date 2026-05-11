# ==============================================================================
# Bootstrap inference pipeline for the EM tenure contamination (epsilon) model
# Created: 2026-05-07
#
# Computes bootstrap standard errors for all 4 eps model variants:
#   - eps_free              (stationary=FALSE, linked=FALSE) — 7 free parameters
#   - eps_stationary        (stationary=TRUE,  linked=FALSE) — 6 free parameters
#   - eps_linked            (stationary=FALSE, linked=TRUE)  — 5 free parameters
#   - eps_stationary_linked (stationary=TRUE,  linked=TRUE)  — 4 free parameters
#
# Design (mirrors bootstrap_pipeline.R at project root):
#   - Nonparametric bootstrap: individuals resampled with replacement.
#   - Weights preserved as columns of the resampled data frame.
#   - B = 200 replicates (labeled in filenames; increase for journal submission).
#   - Warm-start from point estimate; no multi-start per rep.
#   - parallel::mclapply for within-model parallelism.
#   - Per-model .rds files: selective re-running without repeating the full pipeline.
#   - Quality flags: reps failing convergence or with loglik > 50 nats below
#     point estimate are stored but excluded from SE computation.
#
# Prerequisites:
#   EM-tenure/estimate_pipeline.R must have been run first to produce the
#   fit_eps_*.rds point-estimate files in output/results/.
#
# Usage (from project root):
#   Rscript EM-tenure/bootstrap_pipeline_tenure_contamination.R
#
# Outputs (in EM-tenure/output/results/bootstrap/):
#   boot_eps_free_B200.rds
#   boot_eps_stationary_B200.rds
#   boot_eps_linked_B200.rds
#   boot_eps_stationary_linked_B200.rds
#
#   Each .rds file contains a list with:
#     $boot_results   — B-element list of replicate outputs (params, implied, flag)
#     $summary        — data frame from summarise_bootstrap()
#     $point_params   — point-estimate parameter list
#     $point_implied  — point-estimate implied quantities list
#     $point_loglik   — point-estimate log-likelihood
#     $B              — number of bootstrap reps attempted
#     $n_ok           — number of "ok" reps included in SEs
#     $label          — model label string
#     $run_ts         — timestamp of this bootstrap run
# ==============================================================================

library(here)
library(parallel)

# --- Source all required modules ---
# EM-baseline must load FIRST: provides bootstrap_resample() and summarise_bootstrap().
# Its runtime guard calls stop() if 'e_step' already exists in .GlobalEnv, so it
# must precede EM-tenure which defines its own e_step.
source(here::here("EM-baseline", "R", "source_all.R"))
# EM-tenure: eps model estimator + utils (overwrites shared symbols with tenure versions)
source(here::here("EM-tenure", "R", "source_all.R"))

# Guard: .LL_THRESHOLD_EPS must match baseline threshold so quality flags are
# comparable across model families. Both are 50; assert equality at runtime.
stopifnot(".LL_THRESHOLD_EPS must equal baseline .LL_THRESHOLD" =
            .LL_THRESHOLD_EPS == .LL_THRESHOLD)

# --- Bootstrap hyperparameters ---
# IMPORTANT: If B or MASTER_SEED change, delete existing bootstrap files first.
# B and MASTER_SEED are encoded in output filenames; changing either without
# deleting cached files causes silent reuse of stale results.
B           <- 200L
MASTER_SEED <- 42L
N_CORES     <- max(1L, parallel::detectCores() - 1L)

# Estimation caps (must match estimate_pipeline.R)
THETA_CAP <- 0.999
PI_CAP    <- 0.49
EPS_CAP   <- 0.95
EPS_FLOOR <- 1e-4

# Output directory
boot_dir <- here::here("EM-tenure", "output", "results", "bootstrap")
dir.create(boot_dir, recursive = TRUE, showWarnings = FALSE)

run_ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
cat(sprintf("Bootstrap pipeline (eps model) — B=%d | cores=%d | seed=%d\n",
            B, N_CORES, MASTER_SEED))
cat(sprintf("Run timestamp: %s\n\n", run_ts))

# Draw a B × 4 matrix of deterministic per-rep seeds.
# Rows = replicates, columns = model variants (independent per-model seeds).
N_MODELS    <- 4L
set.seed(MASTER_SEED)
seed_matrix <- matrix(
  sample.int(.Machine$integer.max, B * N_MODELS, replace = FALSE),
  nrow = B, ncol = N_MODELS
)

# ==============================================================================
# 1. Load and prepare data
# ==============================================================================

if (!file.exists(here::here("data", "raw", "df_qlfs_A.rds"))) {
  stop("Missing input data: data/raw/df_qlfs_A.rds. ",
       "Obtain the QLFS data and place it at that path before running the pipeline.")
}

source(here::here("scripts", "ingest_data_3waves_SA.R"))

# Convert labelled columns to numeric (same as estimate_pipeline.R).
# Uses base R to avoid loading tidyverse in a reproducibility-critical pipeline.
lbl_cols <- names(df_qlfs)[vapply(df_qlfs, haven::is.labelled, logical(1L))]
if (length(lbl_cols) > 0L)
  df_qlfs[lbl_cols] <- lapply(df_qlfs[lbl_cols], as.numeric)
rm(lbl_cols)

# Never-worked imputation: impute timegap_cat from age for never-worked individuals
.nw_impute_timegap_cat <- function(age_years) {
  dur_months <- pmax((age_years - 16) * 12, 0)
  cut_vals   <- c(0, 3, 6, 9, 12, 36, 60, Inf)
  as.integer(cut(dur_months, breaks = cut_vals, right = FALSE,
                 include.lowest = TRUE))
}
for (.wave in 1:3) {
  nw_col  <- paste0("neverworked", .wave)
  tc_col  <- paste0("timegap_cat", .wave)
  age_col <- paste0("age", .wave)
  if (nw_col %in% names(df_qlfs) && age_col %in% names(df_qlfs)) {
    is_nw <- !is.na(df_qlfs[[nw_col]]) & as.numeric(df_qlfs[[nw_col]]) == 1
    if (any(is_nw))
      df_qlfs[[tc_col]][is_nw] <- .nw_impute_timegap_cat(df_qlfs[[age_col]][is_nw])
  }
}
rm(.wave, .nw_impute_timegap_cat)

# Weight column
if (!"weight" %in% names(df_qlfs))
  df_qlfs$weight <- df_qlfs$weight1

# Validate required columns exist
required_cols <- c("y1", "y2", "y3",
                   "tenure1", "tenure2", "tenure3",
                   "timegap_cat1", "timegap_cat2", "timegap_cat3",
                   "weight")
missing_cols <- setdiff(required_cols, names(df_qlfs))
if (length(missing_cols) > 0L)
  stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))

# Keep only rows used by the eps model (same filter as estimate_pipeline.R).
# complete.cases() on required_cols does a single C-level pass for NA detection,
# avoiding 10 separate is.na() traversals.
keep   <- complete.cases(df_qlfs[, required_cols]) & df_qlfs$weight > 0
df_eps <- as.data.frame(df_qlfs[keep, , drop = FALSE])
rm(df_qlfs, keep)

cat(sprintf("Data loaded: N = %d rows\n\n", nrow(df_eps)))

# ==============================================================================
# 2. Helper: load the most recent fit file for a given prefix
# ==============================================================================

#' Load the most recent .rds file matching a filename prefix pattern.
#'
#' @param results_dir Directory to search.
#' @param prefix      Filename prefix (e.g., "fit_eps_").
#' @return Loaded .rds object.
#' @keywords internal
.load_latest_fit <- function(results_dir, prefix) {
  pattern <- paste0("^", prefix, "\\d{8}_\\d{6}\\.rds$")
  files   <- sort(list.files(results_dir, pattern = pattern, full.names = TRUE))
  if (length(files) == 0L)
    stop(sprintf(
      "No fit file found for prefix '%s' in '%s'.\n",
      prefix, results_dir
    ))
  latest <- tail(files, 1L)
  cat(sprintf("  Loading point estimate: %s\n", basename(latest)))
  readRDS(latest)
}

# ==============================================================================
# 3. Helper: run bootstrap for one model and save
# ==============================================================================

#' Run bootstrap for one eps model variant and save results
#'
#' @param label         Character; model variant label (e.g., "eps_free").
#' @param variant_seeds Integer vector of length B; per-replicate random seeds.
#' @param point_fit     List; point-estimate fit object ($params, $loglik, $gamma).
#' @param stationary    Logical; passed to \code{bootstrap_one_eps()}.
#' @param linked        Logical; passed to \code{bootstrap_one_eps()}.
#'
#' @return Invisible NULL. Side effect: saves
#'   \code{boot_<label>_B<B>_s<MASTER_SEED>.rds} to \code{boot_dir}.
#' @keywords internal
.run_and_save_eps <- function(label, variant_seeds,
                              point_fit, stationary, linked) {
  out_path <- file.path(boot_dir,
    sprintf("boot_%s_B%d_s%d.rds", label, B, MASTER_SEED))

  if (file.exists(out_path)) {
    cat(sprintf("  [SKIP]  %s — already exists: %s\n", label, basename(out_path)))
    return(invisible(NULL))
  }

  point_params  <- point_fit$params
  point_loglik  <- point_fit$loglik
  point_implied <- implied_tenure_contamination(point_params)

  cat(sprintf("  [START] %s ...\n", label))
  t0 <- proc.time()

  boot_results <- parallel::mclapply(
    seq_len(B),
    function(b) {
      bootstrap_one_eps(
        df            = df_eps,
        seed          = variant_seeds[b],
        stationary    = stationary,
        linked        = linked,
        params_start  = point_params,
        point_loglik  = point_loglik,
        theta_cap     = THETA_CAP,
        pi_cap        = PI_CAP,
        eps_cap       = EPS_CAP,
        eps_floor     = EPS_FLOOR,
        max_iter      = 500L    # matches em_fit_tenure_eps() default in estimate_pipeline.R
      )
    },
    mc.cores = N_CORES
  )

  flags <- vapply(boot_results, function(r) r$flag, character(1L))
  n_ok  <- sum(flags == "ok")

  summary_df <- summarise_bootstrap(
    boot_results   = boot_results,
    point_params   = point_params,
    point_implied  = point_implied,
    point_loglik   = point_loglik
  )

  elapsed <- (proc.time() - t0)[["elapsed"]]
  cat(sprintf("  [DONE]  %s — n_ok=%d/%d | %.1fs\n",
              label, n_ok, B, elapsed))

  saveRDS(
    list(
      boot_results  = boot_results,
      summary       = summary_df,
      point_params  = point_params,
      point_implied = point_implied,
      point_loglik  = point_loglik,
      B             = B,
      master_seed   = MASTER_SEED,
      n_ok          = n_ok,
      label         = label,
      run_ts        = run_ts
    ),
    out_path
  )
  invisible(NULL)
}

# ==============================================================================
# 4. Run bootstrap for all 4 eps variants
# ==============================================================================

results_dir <- here::here("output", "results")

eps_configs <- list(
  list(label = "eps_free",              prefix = "fit_eps_",
       stationary = FALSE, linked = FALSE, col_idx = 1L),
  list(label = "eps_stationary",        prefix = "fit_eps_stationary_",
       stationary = TRUE,  linked = FALSE, col_idx = 2L),
  list(label = "eps_linked",            prefix = "fit_eps_linked_",
       stationary = FALSE, linked = TRUE,  col_idx = 3L),
  list(label = "eps_stationary_linked", prefix = "fit_eps_stationary_linked_",
       stationary = TRUE,  linked = TRUE,  col_idx = 4L)
)

# Note: "fit_eps_stationary_" must be listed before "fit_eps_" to avoid
# matching fit_eps_stationary* with the shorter prefix pattern. The
# pattern anchors to ^prefix + timestamp, so there is no ambiguity.

cat("========== EPS MODEL BOOTSTRAP ==========\n")
for (cfg in eps_configs) {
  fit_pt <- .load_latest_fit(results_dir, cfg$prefix)
  .run_and_save_eps(
    label          = cfg$label,
    variant_seeds  = seed_matrix[, cfg$col_idx],
    point_fit      = fit_pt,
    stationary     = cfg$stationary,
    linked         = cfg$linked
  )
}

# ==============================================================================
# 5. Convergence summary
# ==============================================================================

cat("\n========== CONVERGENCE SUMMARY ==========\n")

all_boot_files <- list.files(
  boot_dir,
  pattern    = sprintf("^boot_eps_.*_B%d\\.rds$", B),
  full.names = TRUE
)

if (length(all_boot_files) > 0L) {
  for (f in sort(all_boot_files)) {
    obj <- readRDS(f)
    flags_tbl <- table(vapply(obj$boot_results, `[[`, character(1L), "flag"))
    cat(sprintf("  %-35s  n_ok=%d/%d  flags: %s\n",
                obj$label, obj$n_ok, obj$B,
                paste(names(flags_tbl), flags_tbl, sep = "=", collapse = " ")))
  }
}

cat(sprintf(
  "\nBootstrap complete. Results in:\n  %s\n",
  boot_dir
))
