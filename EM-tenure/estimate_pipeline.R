# ==============================================================================
# EM-tenure: End-to-end estimation pipeline
# ==============================================================================
# This script:
#   1. Sources the EM-tenure module
#   2. Sources the data ingestion script
#   3. Prepares the data (rename weight column if needed)
#   4. Runs the EM with and without misclassification
#   5. Prints a comparison table
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

# --- Prepare data ---
# Ensure weight column exists (the ingest script uses weight1/weight2/weight3;
# we pick weight1 as the baseline weight).
if (!"weight" %in% names(df_qlfs)) {
  df_qlfs$weight <- df_qlfs$weight1
}

# Durations should already be in years from the ingest script.
# Verify columns exist
required_cols <- c("y1", "y2", "y3",
                   "tenure1", "tenure2", "tenure3",
                   "timegap1", "timegap2", "timegap3",
                   "weight")
stopifnot(all(required_cols %in% names(df_qlfs)))

# --- Run EM: with misclassification ---
message("=== Estimating model WITH misclassification ===")

adhoc_filter <- TRUE
if (adhoc_filter == TRUE) {
  n_before <- nrow(df_qlfs)

for (t in 1:3) {
  y <- df_qlfs[[paste0("y", t)]]
  g <- df_qlfs[[paste0("tenure",  t)]]
  d <- df_qlfs[[paste0("timegap", t)]]
  keep <- !((y == 1 & g <= 0) | (y == 0 & d <= 0))
  df_qlfs <- df_qlfs[keep, ]
}

n_after <- nrow(df_qlfs)
message(sprintf("Dropped %d rows with zero/negative duration (%.2f%%)",
                n_before - n_after,
                100 * (n_before - n_after) / n_before))
}
custom_init <- init_params(df_qlfs, misclassification = TRUE)
use_custom <- TRUE
if (use_custom) {
  message("Using custom initial parameters for misclassification model.")
  custom_init$theta1 <- 0.95
  custom_init$theta0 <- 0.05
  custom_init$lambda_g <- 2
  custom_init$lambda_d <- 2
  custom_init$sigma2_g <- 0.5
  custom_init$sigma2_d <- 0.5
  custom_init$pi <- 0.03
} else {
  message("Using default initial parameters for misclassification model.")
}



fit_miscl <- em_fit_tenure(
  df                = df_qlfs,
  params0           = custom_init,
  misclassification = TRUE,
  stationary        = FALSE,
  verbose           = 2L
)
fit_miscl_stationary <- em_fit_tenure(
  df                = df_qlfs,
  params0           = custom_init,
  misclassification = TRUE,
  stationary        = TRUE,
  verbose           = 2L
)


# --- Run EM: without misclassification ---
message("\n=== Estimating model WITHOUT misclassification ===")
fit_no_miscl <- em_fit_tenure(
  df                = df_qlfs,
  params0           = custom_init,
  misclassification = FALSE,
  stationary        = FALSE,
  verbose           = 2L
)
fit_no_miscl_stationary <- em_fit_tenure(
  df                = df_qlfs,
  params0           = custom_init,
  misclassification = FALSE,
  stationary        = TRUE,
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

comparison <- compare_params(fit_miscl, fit_no_miscl)
message("\n=== Parameter comparison ===")
message(paste(capture.output(comparison), collapse = "\n"))
message(sprintf("\nLog-likelihood WITH miscl.:    %.4f", fit_miscl$loglik))
message(sprintf("Log-likelihood WITHOUT miscl.: %.4f", fit_no_miscl$loglik))
