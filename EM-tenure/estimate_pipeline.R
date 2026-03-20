# ==============================================================================
# EM-tenure: End-to-end estimation pipeline
# ==============================================================================
# This script:
#   1. Sources the EM-tenure module
#   2. Sources the data ingestion script
#   3. Prepares the data (rename weight column if needed)
#   4. Runs the EM with and without misclassification
#   5. Prints a comparison table
#   6. Saves each fit as a timestamped .rds in output/results/
#   7. Appends a one-row summary per model to output/results/run_summary.csv
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
df_em <- df_qlfs |>
  filter(!(!is.na(neverworked1) & as.numeric(neverworked1) == 1))
df_em <- df_em |>
  filter(!(!is.na(neverworked2) & as.numeric(neverworked2) == 1))
df_em <- df_em |>
  filter(!(!is.na(neverworked3) & as.numeric(neverworked3) == 1))
df_qlfs <- df_em
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

# --- Run EM: with misclassification ---
message("=== Estimating model WITH misclassification ===")

# NOTE: The ad-hoc zero-duration filter has been removed. The ingest script
# now performs nearest-non-zero imputation (Issue 1 fix), so all tenure
# values for employed waves and all timegap_cat values are valid before
# reaching this pipeline. The never-worked filter above removes individuals
# who were never employed, eliminating the source of structural zeros.

custom_init <- init_params(df_qlfs, misclassification = TRUE, discrete_timegap = TRUE)
use_custom <- FALSE
if (use_custom) {
  message("Using custom initial parameters for misclassification model.")
  custom_init$theta1 <- 0.95
  custom_init$theta0 <- 0.05
  custom_init$lambda_g <- 2
  custom_init$lambda_d <- 2
  custom_init$sigma2_g <- 0.5
  custom_init$pi <- 0.03
} else {
  message("Using default initial parameters for misclassification model.")
}



# --- Results storage setup ---
# run_id ties together the .rds files and the summary CSV row for this run.
set.seed(2026L)  # reproducibility — placed here so all four fits are seeded
run_id     <- format(Sys.time(), "%Y%m%d_%H%M%S")
results_dir <- "output/results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

fit_miscl_em <- em_fit_tenure(
  df                = df_qlfs,
  params0           = custom_init,
  misclassification = TRUE,
  stationary        = FALSE,
  discrete_timegap  = TRUE,
  verbose           = 2L
)
fit_miscl_stationary_em <- em_fit_tenure(
  df                = df_qlfs,
  params0           = custom_init,
  misclassification = TRUE,
  stationary        = TRUE,
  discrete_timegap  = TRUE,
  verbose           = 2L
)


# --- Run EM: without misclassification ---
message("\n=== Estimating model WITHOUT misclassification ===")
fit_no_miscl_em <- em_fit_tenure(
  df                = df_qlfs,
  params0           = custom_init,
  misclassification = FALSE,
  stationary        = FALSE,
  discrete_timegap  = TRUE,
  verbose           = 2L
)
fit_no_miscl_stationary_em <- em_fit_tenure(
  df                = df_qlfs,
  params0           = custom_init,
  misclassification = FALSE,
  stationary        = TRUE,
  discrete_timegap  = TRUE,
  verbose           = 2L
)


# --- Comparison table ---
compare_params <- function(fit1, fit2, label1 = "With miscl.", label2 = "No miscl.") {
  p1 <- unlist(fit1$params)
  p2 <- unlist(fit2$params)
  data.frame(
    parameter = names(p1),
    !!label1 := p1,
    !!label2 := p2,
    row.names = NULL
  )
}

comparison <- compare_params(fit_miscl_em, fit_no_miscl_em)
message("\n=== Parameter comparison ===")
message(paste(capture.output(comparison), collapse = "\n"))
message(sprintf("\nLog-likelihood WITH miscl.:    %.4f", fit_miscl_em$loglik))
message(sprintf("Log-likelihood WITHOUT miscl.: %.4f", fit_no_miscl_em$loglik))

# --- Save results ---
# Each fit is saved as a self-contained .rds (includes gamma, history, params).
# A flat summary row is appended to run_summary.csv for quick comparison
# across runs without loading the full .rds files.

fits <- list(
  miscl             = fit_miscl_em,
  miscl_stationary  = fit_miscl_stationary_em,
  no_miscl          = fit_no_miscl_em,
  no_miscl_stationary = fit_no_miscl_stationary_em
)

for (model_name in names(fits)) {
  rds_path <- file.path(results_dir, sprintf("fit_%s_%s.rds", model_name, run_id))
  saveRDS(fits[[model_name]], rds_path)
  message(sprintf("Saved: %s", rds_path))
}

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
