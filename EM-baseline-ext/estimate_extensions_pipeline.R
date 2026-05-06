# ==============================================================================
# EM-baseline-ext: Estimation pipeline for baseline model extensions
# Created: 2026-05-06
#
# Estimates all 18 extension model configurations:
#   Extension I  (covariates):    3 sets × {symmetric, none} × {stat, free} = 12
#   Extension III (FMM):          {symmetric, none} × {stat, free} = 4
#   Extension IV  (inconsistency): {symmetric only} × {stat, free} = 2
#
# Outputs (per model):
#   EM-baseline-ext/output/results/fit_<label>.rds
# Outputs (aggregate):
#   EM-baseline-ext/output/results/run_summary.csv
#
# Usage (from project root):
#   Rscript EM-baseline-ext/estimate_extensions_pipeline.R
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------

library(here)
# Note: library(dplyr) is required by scripts/ingest_data_3waves_SA.R (uses |> with dplyr verbs).
library(dplyr)
library(tidyverse)

# Source baseline first (provides em_fit_baseline + shared utils)
source(here::here("EM-baseline", "R", "source_all.R"))

# Source extension modules (idempotent guard skips already-loaded baseline utils)
source(here::here("EM-baseline-ext", "R", "source_all.R"))

# Estimation hyperparameters
THETA_CAP   <- 0.999
PI_CAP      <- 0.49
N_STARTS    <- 5L
PERTURB_SD  <- 0.30
RANDOM_SEED <- 1234L
set.seed(RANDOM_SEED)

# Output directories
results_dir <- here::here("EM-baseline-ext", "output", "results")
tables_dir  <- here::here("EM-baseline-ext", "output", "tables")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir,  recursive = TRUE, showWarnings = FALSE)

run_ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

# ---- Parameter perturbation helpers ----------------------------------------
.perturb_p <- function(p, sd = PERTURB_SD) inv_logit(logit(p) + rnorm(1, 0, sd))
.clamp_pi  <- function(p, cap = PI_CAP) min(cap - 0.01, max(1e-6, .perturb_p(p)))
.perturb_beta <- function(beta, sd = PERTURB_SD) {
  beta + rnorm(length(beta), 0, sd)
}

# ---- Multi-start runner: covariate extension --------------------------------
.run_cov <- function(df, X, model_type, stationary, label) {
  best <- NULL
  for (s in seq_len(N_STARTS)) {
    p0 <- init_params_covariates(ncol(X), model_type)
    p0$beta0 <- .perturb_beta(p0$beta0)
    p0$beta1 <- .perturb_beta(p0$beta1)
    if (!is.null(p0$pi)) p0$pi <- .clamp_pi(p0$pi)
    fit <- tryCatch(
      em_fit_covariates(df, X, model_type = model_type, stationary = stationary,
                        params0 = p0, max_iter = 1000L, tol = 1e-8,
                        pi_cap = PI_CAP, verbose = 0L),
      error = function(e) { warning(sprintf("[%s] start %d: %s", label, s, e$message)); NULL }
    )
    if (!is.null(fit) && (is.null(best) || fit$loglik > best$loglik)) best <- fit
  }
  if (is.null(best)) stop(sprintf("All %d starts failed for %s", N_STARTS, label))
  saveRDS(best, file.path(results_dir, sprintf("fit_%s.rds", label)))
  cat(sprintf("  [%s] loglik=%.4f | converged=%s | iters=%d\n",
              label, best$loglik, best$converged, best$iterations))
  cat(sprintf("    beta0_intercept=%.4f  beta1_intercept=%.4f",
              best$params$beta0[1L], best$params$beta1[1L]))
  if (!is.null(best$params$pi)) cat(sprintf("  pi=%.4f", best$params$pi))
  cat("\n")
  best
}

# ---- Multi-start runner: FMM extension -------------------------------------
.run_fmm <- function(df, model_type, stationary, label) {
  best <- NULL
  for (s in seq_len(N_STARTS)) {
    p0 <- init_params_fmm(model_type, stationary)
    p0$theta0_A <- .perturb_p(p0$theta0_A)
    p0$theta1_A <- .perturb_p(p0$theta1_A)
    p0$theta0_B <- .perturb_p(p0$theta0_B)
    p0$theta1_B <- .perturb_p(p0$theta1_B)
    p0$phi      <- .perturb_p(p0$phi)
    if (!is.null(p0$pi)) p0$pi <- .clamp_pi(p0$pi)
    fit <- tryCatch(
      em_fit_fmm(df, model_type = model_type, stationary = stationary,
                 params0 = p0, max_iter = 1000L, tol = 1e-8,
                 theta_cap = THETA_CAP, pi_cap = PI_CAP, verbose = 0L),
      error = function(e) { warning(sprintf("[%s] start %d: %s", label, s, e$message)); NULL }
    )
    if (!is.null(fit) && (is.null(best) || fit$loglik > best$loglik)) best <- fit
  }
  if (is.null(best)) stop(sprintf("All %d starts failed for %s", N_STARTS, label))
  saveRDS(best, file.path(results_dir, sprintf("fit_%s.rds", label)))
  p <- best$params
  cat(sprintf("  [%s] loglik=%.4f | converged=%s | iters=%d\n",
              label, best$loglik, best$converged, best$iterations))
  cat(sprintf("    TypeA: theta0=%.4f theta1=%.4f alpha=%.4f\n",
              p$theta0_A, p$theta1_A, p$alpha_A))
  cat(sprintf("    TypeB: theta0=%.4f theta1=%.4f alpha=%.4f\n",
              p$theta0_B, p$theta1_B, p$alpha_B))
  cat(sprintf("    phi=%.4f", p$phi))
  if (!is.null(p$pi)) cat(sprintf("  pi=%.4f", p$pi))
  cat("\n")
  best
}

# ---- Multi-start runner: inconsistency extension ---------------------------
.run_incons <- function(df, inc_mat, stationary, label) {
  best <- NULL
  for (s in seq_len(N_STARTS)) {
    p0 <- init_params_inconsistency(stationary)
    p0$theta0    <- .perturb_p(p0$theta0)
    p0$theta1    <- .perturb_p(p0$theta1)
    p0$delta[1L] <- p0$delta[1L] + rnorm(1, 0, PERTURB_SD)
    # slopes start near 0 with small noise
    p0$delta[2L] <- rnorm(1, 0, 0.10)
    p0$delta[3L] <- rnorm(1, 0, 0.10)
    fit <- tryCatch(
      em_fit_inconsistency(df, inc_mat, stationary = stationary,
                           params0 = p0, max_iter = 1000L, tol = 1e-8,
                           theta_cap = THETA_CAP, verbose = 0L),
      error = function(e) { warning(sprintf("[%s] start %d: %s", label, s, e$message)); NULL }
    )
    if (!is.null(fit) && (is.null(best) || fit$loglik > best$loglik)) best <- fit
  }
  if (is.null(best)) stop(sprintf("All %d starts failed for %s", N_STARTS, label))
  saveRDS(best, file.path(results_dir, sprintf("fit_%s.rds", label)))
  p <- best$params
  cat(sprintf("  [%s] loglik=%.4f | converged=%s | iters=%d\n",
              label, best$loglik, best$converged, best$iterations))
  cat(sprintf("    theta0=%.4f theta1=%.4f alpha=%.4f\n", p$theta0, p$theta1, p$alpha))
  cat(sprintf("    delta=(%s)  pi_base=%.4f\n",
              paste(round(p$delta, 4L), collapse=", "),
              0.5 * plogis(p$delta[1L])))
  best
}

# ------------------------------------------------------------------------------
# 1. Ingest and prepare data
# ------------------------------------------------------------------------------

source(here::here("scripts", "ingest_data_3waves_SA.R"))  # loads df_qlfs

cat(sprintf("Raw data loaded: N = %d\n", nrow(df_qlfs)))

# Validate employment columns
for (y_col in c("y1", "y2", "y3")) {
  if (is.factor(df_qlfs[[y_col]]))
    stop(sprintf("Column '%s' is a factor — convert to integer 0/1 in ingest script.", y_col))
  non_na <- df_qlfs[[y_col]][!is.na(df_qlfs[[y_col]])]
  if (!all(non_na %in% c(0, 1)))
    stop(sprintf("Column '%s' has non-binary values.", y_col))
}

# Validate weight column type before coercion
if (!is.numeric(df_qlfs$weight) && !is.integer(df_qlfs$weight))
  stop("Column 'weight' is not numeric/integer. Check the ingest script.")

# Keep rows needed for extension models (require age + educ for inconsistency/covariates)
# Also require strictly positive weight.
keep <- !is.na(df_qlfs$y1) & !is.na(df_qlfs$y2) & !is.na(df_qlfs$y3) &
        !is.na(df_qlfs$weight) & df_qlfs$weight > 0 &
        !is.na(df_qlfs$age1)  & !is.na(df_qlfs$age2)  & !is.na(df_qlfs$age3) &
        !is.na(df_qlfs$educ1) & !is.na(df_qlfs$educ2) & !is.na(df_qlfs$educ3)
df_ext        <- df_qlfs[keep, , drop = FALSE]
df_ext$y1     <- as.integer(df_ext$y1)
df_ext$y2     <- as.integer(df_ext$y2)
df_ext$y3     <- as.integer(df_ext$y3)
df_ext$weight <- as.numeric(df_ext$weight)
# contracttype1 is NA for non-employed — recode to 0 (no contract)
df_ext$contracttype1 <- ifelse(is.na(df_ext$contracttype1), 0L,
                               as.integer(df_ext$contracttype1))
df_ext        <- as.data.frame(df_ext)

cat(sprintf("Analysis sample (complete age/educ): N = %d\n", nrow(df_ext)))

stopifnot(
  all(df_ext$y1 %in% c(0L, 1L)),
  all(df_ext$y2 %in% c(0L, 1L)),
  all(df_ext$y3 %in% c(0L, 1L)),
  all(df_ext$weight > 0)
)

# ---- Prepare inconsistency matrix ------------------------------------------
df_incons <- compute_inconsistencies(df_ext)
inc_mat   <- as.matrix(df_incons[, c("Y_age_1", "Y_age_2", "Y_age_3",
                                      "Y_edu_1", "Y_edu_2", "Y_edu_3")])

cat(sprintf("Inconsistency indicators: %d observations with at least one age/edu flag\n",
            sum(rowSums(inc_mat) > 0)))

# ---- Prepare covariate matrices (3 sets) -----------------------------------
cv_set1 <- prepare_covariate_matrix(df_ext, covariate_set = 1L)
cv_set2 <- prepare_covariate_matrix(df_ext, covariate_set = 2L)
cv_set3 <- prepare_covariate_matrix(df_ext, covariate_set = 3L)

X1 <- cv_set1$X
X2 <- cv_set2$X
X3 <- cv_set3$X

cat(sprintf("Covariate matrices: Set1 p=%d | Set2 p=%d | Set3 p=%d\n",
            ncol(X1), ncol(X2), ncol(X3)))

# Free large raw object
rm(df_qlfs, df_incons)

# Container for all fits and run summary
fits     <- list()
run_rows <- list()

# Helper: build a run summary row (returns a data.frame; caller accumulates)
.make_row <- function(label, family, model_type, stationary, fit) {
  data.frame(
    timestamp   = run_ts,
    family      = family,
    model_type  = model_type,
    stationary  = stationary,
    label       = label,
    converged   = fit$converged,
    iterations  = fit$iterations,
    loglik      = fit$loglik,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# SECTION 1: Extension I — Observable heterogeneity (covariate models)
#   12 models: Set {1,2,3} × model_type {symmetric, none} × {stat, free}
# ==============================================================================

cat("\n\n========== EXTENSION I: COVARIATE MODELS ==========\n")

# --- Set 1: intercept + age + age^2 + educ ------------------------------------

cat("\n--- Set 1 (age + educ) ---\n")

fits[["cov_s1_sym_stat"]] <- .run_cov(df_ext, X1, "symmetric", TRUE,  "cov_s1_sym_stat")
run_rows[["cov_s1_sym_stat"]] <- .make_row("cov_s1_sym_stat", "covariate", "symmetric", TRUE,  fits[["cov_s1_sym_stat"]])

fits[["cov_s1_sym_free"]] <- .run_cov(df_ext, X1, "symmetric", FALSE, "cov_s1_sym_free")
run_rows[["cov_s1_sym_free"]] <- .make_row("cov_s1_sym_free", "covariate", "symmetric", FALSE, fits[["cov_s1_sym_free"]])

fits[["cov_s1_non_stat"]] <- .run_cov(df_ext, X1, "none",      TRUE,  "cov_s1_non_stat")
run_rows[["cov_s1_non_stat"]] <- .make_row("cov_s1_non_stat", "covariate", "none",      TRUE,  fits[["cov_s1_non_stat"]])

fits[["cov_s1_non_free"]] <- .run_cov(df_ext, X1, "none",      FALSE, "cov_s1_non_free")
run_rows[["cov_s1_non_free"]] <- .make_row("cov_s1_non_free", "covariate", "none",      FALSE, fits[["cov_s1_non_free"]])

# --- Set 2: Set 1 + race + female ---------------------------------------------

cat("\n--- Set 2 (age + educ + race + female) ---\n")

fits[["cov_s2_sym_stat"]] <- .run_cov(df_ext, X2, "symmetric", TRUE,  "cov_s2_sym_stat")
run_rows[["cov_s2_sym_stat"]] <- .make_row("cov_s2_sym_stat", "covariate", "symmetric", TRUE,  fits[["cov_s2_sym_stat"]])

fits[["cov_s2_sym_free"]] <- .run_cov(df_ext, X2, "symmetric", FALSE, "cov_s2_sym_free")
run_rows[["cov_s2_sym_free"]] <- .make_row("cov_s2_sym_free", "covariate", "symmetric", FALSE, fits[["cov_s2_sym_free"]])

fits[["cov_s2_non_stat"]] <- .run_cov(df_ext, X2, "none",      TRUE,  "cov_s2_non_stat")
run_rows[["cov_s2_non_stat"]] <- .make_row("cov_s2_non_stat", "covariate", "none",      TRUE,  fits[["cov_s2_non_stat"]])

fits[["cov_s2_non_free"]] <- .run_cov(df_ext, X2, "none",      FALSE, "cov_s2_non_free")
run_rows[["cov_s2_non_free"]] <- .make_row("cov_s2_non_free", "covariate", "none",      FALSE, fits[["cov_s2_non_free"]])

# --- Set 3: Set 2 + contract type ---------------------------------------------

cat("\n--- Set 3 (age + educ + race + female + contracttype) ---\n")

fits[["cov_s3_sym_stat"]] <- .run_cov(df_ext, X3, "symmetric", TRUE,  "cov_s3_sym_stat")
run_rows[["cov_s3_sym_stat"]] <- .make_row("cov_s3_sym_stat", "covariate", "symmetric", TRUE,  fits[["cov_s3_sym_stat"]])

fits[["cov_s3_sym_free"]] <- .run_cov(df_ext, X3, "symmetric", FALSE, "cov_s3_sym_free")
run_rows[["cov_s3_sym_free"]] <- .make_row("cov_s3_sym_free", "covariate", "symmetric", FALSE, fits[["cov_s3_sym_free"]])

fits[["cov_s3_non_stat"]] <- .run_cov(df_ext, X3, "none",      TRUE,  "cov_s3_non_stat")
run_rows[["cov_s3_non_stat"]] <- .make_row("cov_s3_non_stat", "covariate", "none",      TRUE,  fits[["cov_s3_non_stat"]])

fits[["cov_s3_non_free"]] <- .run_cov(df_ext, X3, "none",      FALSE, "cov_s3_non_free")
run_rows[["cov_s3_non_free"]] <- .make_row("cov_s3_non_free", "covariate", "none",      FALSE, fits[["cov_s3_non_free"]])

# ==============================================================================
# SECTION 2: Extension III — 2-type FMM
#   4 models: model_type {symmetric, none} × {stat, free}
# ==============================================================================

cat("\n\n========== EXTENSION III: FMM MODELS ==========\n")

fits[["fmm_sym_stat"]] <- .run_fmm(df_ext, "symmetric", TRUE,  "fmm_sym_stat")
run_rows[["fmm_sym_stat"]] <- .make_row("fmm_sym_stat", "fmm", "symmetric", TRUE,  fits[["fmm_sym_stat"]])

fits[["fmm_sym_free"]] <- .run_fmm(df_ext, "symmetric", FALSE, "fmm_sym_free")
run_rows[["fmm_sym_free"]] <- .make_row("fmm_sym_free", "fmm", "symmetric", FALSE, fits[["fmm_sym_free"]])

fits[["fmm_non_stat"]] <- .run_fmm(df_ext, "none",      TRUE,  "fmm_non_stat")
run_rows[["fmm_non_stat"]] <- .make_row("fmm_non_stat", "fmm", "none",      TRUE,  fits[["fmm_non_stat"]])

fits[["fmm_non_free"]] <- .run_fmm(df_ext, "none",      FALSE, "fmm_non_free")
run_rows[["fmm_non_free"]] <- .make_row("fmm_non_free", "fmm", "none",      FALSE, fits[["fmm_non_free"]])

# ==============================================================================
# SECTION 3: Extension IV — Inconsistency-augmented misclassification
#   2 models: {stat, free} (symmetric only)
# ==============================================================================

cat("\n\n========== EXTENSION IV: INCONSISTENCY MODELS ==========\n")

fits[["incons_sym_stat"]] <- .run_incons(df_ext, inc_mat, TRUE,  "incons_sym_stat")
run_rows[["incons_sym_stat"]] <- .make_row("incons_sym_stat", "inconsistency", "symmetric", TRUE,  fits[["incons_sym_stat"]])

fits[["incons_sym_free"]] <- .run_incons(df_ext, inc_mat, FALSE, "incons_sym_free")
run_rows[["incons_sym_free"]] <- .make_row("incons_sym_free", "inconsistency", "symmetric", FALSE, fits[["incons_sym_free"]])

# ==============================================================================
# SECTION 4: Cross-extension comparison and summary
# ==============================================================================

cat("\n\n========== SUMMARY ==========\n")

run_summary <- do.call(rbind, run_rows)
write.csv(run_summary,
          file.path(results_dir, "run_summary.csv"),
          row.names = FALSE)

# Print compact table
cat(sprintf("\n%-24s %10s %10s %10s\n", "Label", "LogLik", "Converged", "Iters"))
cat(strrep("-", 58), "\n")
for (lbl in names(fits)) {
  fit <- fits[[lbl]]
  cat(sprintf("%-24s %10.4f %10s %10d\n",
              lbl, fit$loglik, fit$converged, fit$iterations))
}

# LR tests within covariate family (Set 1 vs Set 2 vs Set 3, symmetric+stationary)
cat("\n--- LR tests: covariate model complexity (symmetric, stationary) ---\n")
.lr_test <- function(fit_null, fit_alt, df_diff, label) {
  if (!is.numeric(df_diff) || df_diff <= 0L)
    stop(sprintf(".lr_test: df_diff=%s must be > 0 for '%s'. Check X matrix dimensions.",
                 df_diff, label))
  lr_stat <- 2 * (fit_alt$loglik - fit_null$loglik)
  p_val   <- pchisq(lr_stat, df_diff, lower.tail = FALSE)
  cat(sprintf("  %-30s LR=%.3f df=%d p=%.4f\n", label, lr_stat, df_diff, p_val))
}

.lr_test(fits[["cov_s1_sym_stat"]], fits[["cov_s2_sym_stat"]],
         ncol(X2) - ncol(X1), "Set1 vs Set2 (sym, stat)")
.lr_test(fits[["cov_s2_sym_stat"]], fits[["cov_s3_sym_stat"]],
         ncol(X3) - ncol(X2), "Set2 vs Set3 (sym, stat)")

# Stationarity LR test within each extension family (symmetric only)
cat("\n--- LR tests: free vs stationary alpha (symmetric) ---\n")

.lr_test(fits[["cov_s1_sym_stat"]], fits[["cov_s1_sym_free"]], 1L, "Cov Set1 stat vs free")
.lr_test(fits[["cov_s2_sym_stat"]], fits[["cov_s2_sym_free"]], 1L, "Cov Set2 stat vs free")
.lr_test(fits[["cov_s3_sym_stat"]], fits[["cov_s3_sym_free"]], 1L, "Cov Set3 stat vs free")
.lr_test(fits[["fmm_sym_stat"]],    fits[["fmm_sym_free"]],    1L, "FMM stat vs free")
.lr_test(fits[["incons_sym_stat"]], fits[["incons_sym_free"]], 1L, "Inconsistency stat vs free")

# Error vs no-error LR tests (stationary only, Set 1 for covariates)
cat("\n--- LR tests: symmetric vs no-error (stationary) ---\n")

.lr_test(fits[["cov_s1_non_stat"]], fits[["cov_s1_sym_stat"]], 1L, "Cov Set1: none vs sym")
.lr_test(fits[["cov_s2_non_stat"]], fits[["cov_s2_sym_stat"]], 1L, "Cov Set2: none vs sym")
.lr_test(fits[["cov_s3_non_stat"]], fits[["cov_s3_sym_stat"]], 1L, "Cov Set3: none vs sym")
.lr_test(fits[["fmm_non_stat"]],    fits[["fmm_sym_stat"]],    1L, "FMM: none vs sym")

cat(sprintf("\nRun complete. Results in: %s\n", results_dir))
