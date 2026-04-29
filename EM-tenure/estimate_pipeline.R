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

# Warm-start: interactively select a previous fit_miscl run to warm-start from.
# Shows all available results with their log-likelihood and iteration count,
# then prompts the user to choose one (or use default initial parameters).
.prev_rds_files <- sort(list.files(
  "output/results",
  pattern    = "^fit_miscl_\\d{8}_\\d{6}\\.rds$",
  full.names = TRUE
))
if (length(.prev_rds_files) > 0) {
  # Build a display table: filename | loglik | iterations | converged
  .prev_meta <- lapply(.prev_rds_files, function(f) {
    tryCatch({
      obj <- readRDS(f)
      list(
        path      = f,
        label     = sprintf(
          "%-40s  loglik = %10.4f  iter = %3d  converged = %s",
          basename(f), obj$loglik, obj$iterations,
          if (isTRUE(obj$converged)) "YES" else "NO "
        )
      )
    }, error = function(e) {
      list(path = f, label = sprintf("%-40s  [unreadable]", basename(f)))
    })
  })

  message("\nAvailable warm-start candidates:")
  for (.i in seq_along(.prev_meta)) {
    message(sprintf("  [%d] %s", .i, .prev_meta[[.i]]$label))
  }

  .choice <- menu(
    choices = c(
      vapply(.prev_meta, `[[`, character(1), "label"),
      "--- Use default initial parameters (no warm-start) ---"
    ),
    title = "\nSelect a warm-start file (or choose the last option to skip):"
  )

  if (.choice > 0 && .choice <= length(.prev_meta)) {
    .prev_rds <- .prev_meta[[.choice]]$path
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
    rm(.prev_fit, .k)
  } else {
    use_custom <- FALSE
    message("Using default initial parameters (no warm-start selected).")
  }
  rm(.prev_rds_files, .prev_meta, .choice)
} else {
  use_custom <- FALSE
  message("No previous results found — using default initial parameters.")
}
if (use_custom) {
  message("Warm-starting from selected previous parameters.")
} else {
  message("Using default initial parameters.")
}



# --- Results storage setup ---
# run_id ties together the .rds files and the summary CSV row for this run.
set.seed(1234)  # reproducibility — placed here so both fits are seeded
run_id     <- format(Sys.time(), "%Y%m%d_%H%M%S")
results_dir <- "output/results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
summary_csv <- file.path(results_dir, "run_summary.csv")

message("=== Estimating non-stationary model (free) ===")
fit_em <- em_fit_tenure(
  df                = df_qlfs,
  params0           = custom_init,
  stationary        = FALSE,
  linked            = FALSE,
  discrete_timegap  = TRUE,
  verbose           = 2L
)
message("=== Estimating stationary model (free) ===")
fit_stationary_em <- em_fit_tenure(
  df                = df_qlfs,
  params0           = custom_init,
  stationary        = TRUE,
  linked            = FALSE,
  discrete_timegap  = TRUE,
  verbose           = 2L
)

# --- CTMC-linked estimation ---
# The linked specification constrains lambda_g = -log(theta1)/Delta and
# lambda_d = -log(1-theta0)/Delta via the CTMC link, so durations and
# transitions are jointly identified from a single pair (theta1, theta0).
# Initialise from the free-specification result to warm-start.
linked_init <- custom_init
linked_init$theta1  <- fit_em$params$theta1
linked_init$theta0  <- fit_em$params$theta0
linked_init$alpha   <- fit_em$params$alpha
linked_init$pi      <- fit_em$params$pi
linked_init$sigma2_g <- fit_em$params$sigma2_g
# Lambda will be set from theta inside em_fit_tenure when linked=TRUE
linked_init2 <- init_params(df_qlfs, discrete_timegap = TRUE, linked = TRUE)

message("=== Estimating non-stationary model (CTMC-linked) ===")
fit_linked <- em_fit_tenure(
  df                = df_qlfs,
  params0           = linked_init,
  stationary        = FALSE,
  linked            = TRUE,
  discrete_timegap  = TRUE,
  verbose           = 2L
)
fit_linked_data_driven_params <- em_fit_tenure(
  df                = df_qlfs,
  params0           = linked_init2,
  stationary        = FALSE,
  linked            = TRUE,
  discrete_timegap  = TRUE,
  verbose           = 2L
)
message("=== Estimating stationary model (CTMC-linked) ===")
fit_stationary_linked <- em_fit_tenure(
  df                = df_qlfs,
  params0           = linked_init,
  stationary        = TRUE,
  linked            = TRUE,
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
    rho         = NA_real_,   # base model has no rho; keeps CSV columns aligned
    sigma2_g    = p$sigma2_g,
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

# Warm-start from the free non-stationary fit
rho_init <- fit_em$params
rho_init$rho <- 0.15  # initial contamination rate

fit_rho <- em_fit_tenure_rho(
  df         = df_qlfs,
  params0    = rho_init,
  stationary = FALSE,
  linked     = FALSE,
  verbose    = 2L
)

message("\n=== Estimating rho model (free, stationary) ===")
fit_rho_stationary <- em_fit_tenure_rho(
  df         = df_qlfs,
  params0    = rho_init,
  stationary = TRUE,
  linked     = FALSE,
  verbose    = 2L
)

message("\n=== Estimating rho model (CTMC-linked, non-stationary) ===")
rho_linked_init <- rho_init
rho_linked_init$theta1  <- fit_rho$params$theta1
rho_linked_init$theta0  <- fit_rho$params$theta0

fit_rho_linked <- em_fit_tenure_rho(
  df         = df_qlfs,
  params0    = rho_linked_init,
  stationary = FALSE,
  linked     = TRUE,
  verbose    = 2L
)

message("\n=== Estimating rho model (CTMC-linked, stationary) ===")
fit_rho_stationary_linked <- em_fit_tenure_rho(
  df         = df_qlfs,
  params0    = rho_linked_init,
  stationary = TRUE,
  linked     = TRUE,
  verbose    = 2L
)

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

# Summary rows for rho models
.make_summary_row_rho <- function(run_id, model, fit) {
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
    rho         = p$rho,
    sigma2_g    = p$sigma2_g,
    lambda_g    = p$lambda_g,
    lambda_d    = p$lambda_d,
    sigma2_d    = NA_real_   # rho model uses discrete timegap; sigma2_d not estimated
  )
}

rho_fits <- list(
  rho_free              = fit_rho,
  rho_stationary        = fit_rho_stationary,
  rho_linked            = fit_rho_linked,
  rho_stationary_linked = fit_rho_stationary_linked
)

rho_summary <- data.table::rbindlist(
  lapply(names(rho_fits), function(m) .make_summary_row_rho(run_id, m, rho_fits[[m]]))
)

data.table::fwrite(
  rho_summary,
  summary_csv,
  append    = TRUE,
  col.names = FALSE
)
message(sprintf("Rho summary appended: %s  [run_id = %s]", summary_csv, run_id))
