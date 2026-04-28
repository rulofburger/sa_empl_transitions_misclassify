# ==============================================================================
# EM-tenure: End-to-end estimation pipeline
# ==============================================================================
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

# --- Load EM-tenure module ---
source("EM-tenure/R/source_all.R")

# --- Ingest data ---
# Always re-ingest to avoid stale data from a previous session.
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

custom_init <- init_params(df_qlfs, discrete_timegap = TRUE, linked = FALSE)

# Warm-start: load converged params from the most recent fit_miscl run, if available.
# Falls back to init_params() defaults when no previous results exist.
.prev_rds_files <- list.files(
  "output/results",
  pattern    = "^fit_miscl_\\d{8}_\\d{6}\\.rds$",
  full.names = TRUE
)
if (length(.prev_rds_files) > 0) {
  .prev_rds <- sort(.prev_rds_files) |> tail(1)
  .prev_fit <- readRDS(.prev_rds)
  # Only overwrite keys present in custom_init to stay compatible with
  # stationary and non-stationary models downstream.
  for (.k in intersect(names(.prev_fit$params), names(custom_init))) {
    custom_init[[.k]] <- .prev_fit$params[[.k]]
  }
  use_custom <- TRUE
  message(sprintf(
    "Warm-start: loaded params from %s  [loglik = %.4f, iter = %d]",
    basename(.prev_rds), .prev_fit$loglik, .prev_fit$iterations
  ))
  rm(.prev_rds_files, .prev_rds, .prev_fit, .k)
} else {
  use_custom <- FALSE
  message("No previous results found — using default initial parameters.")
  rm(.prev_rds_files)
}
if (use_custom) {
  message("Using custom initial parameters.")
  custom_init$theta1 <- 0.95
  custom_init$theta0 <- 0.05
  custom_init$lambda_g <- 2
  custom_init$lambda_d <- 2
  custom_init$sigma2_g <- 0.5
  custom_init$pi <- 0.03
} else {
  message("Using default initial parameters.")
}



# --- Results storage setup ---
# run_id ties together the .rds files and the summary CSV row for this run.
set.seed(1234)  # reproducibility — placed here so both fits are seeded
run_id     <- format(Sys.time(), "%Y%m%d_%H%M%S")
results_dir <- "output/results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

message("=== Estimating non-stationary model ===")
fit_em <- em_fit_tenure(
  df                = df_qlfs,
  params0           = custom_init,
  stationary        = FALSE,
  linked            = FALSE,
  discrete_timegap  = TRUE,
  verbose           = 2L
)
message("=== Estimating stationary model ===")
fit_stationary_em <- em_fit_tenure(
  df                = df_qlfs,
  params0           = custom_init,
  stationary        = TRUE,
  linked            = FALSE,
  discrete_timegap  = TRUE,
  verbose           = 2L
)

# --- Save results ---
# Each fit is saved as a self-contained .rds (includes gamma, history, params).
# A flat summary row is appended to run_summary.csv for quick comparison
# across runs without loading the full .rds files.
saveRDS(fit_em,
        file.path(results_dir, sprintf("fit_miscl_%s.rds", run_id)))
saveRDS(fit_stationary_em,
        file.path(results_dir, sprintf("fit_miscl_stationary_%s.rds", run_id)))
message(sprintf("Fits saved [run_id = %s]", run_id))

message(sprintf("\nLog-likelihood (non-stationary): %.4f", fit_em$loglik))
message(sprintf("Log-likelihood (stationary):     %.4f", fit_stationary_em$loglik))

fits <- list(
  miscl             = fit_em,
  miscl_stationary  = fit_stationary_em
)

# Flat summary — one row per model per run
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
    sigma2_g    = p$sigma2_g,
    lambda_g    = p$lambda_g,
    lambda_d    = p$lambda_d,
    sigma2_d    = p$sigma2_d %||% NA_real_
  )
}

summary_rows <- data.table::rbindlist(
  lapply(names(fits), function(m) .make_summary_row(run_id, m, fits[[m]]))
)

summary_csv <- file.path(results_dir, "run_summary.csv")
data.table::fwrite(
  summary_rows,
  summary_csv,
  append    = file.exists(summary_csv),
  col.names = !file.exists(summary_csv)
)
message(sprintf("Summary appended: %s  [run_id = %s]", summary_csv, run_id))
