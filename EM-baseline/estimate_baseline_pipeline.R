# ==============================================================================
# EM-baseline: End-to-end estimation pipeline
# Created: 2026-05-05
#
# Estimates all six baseline EM model configurations on the QLFS 3-wave panel:
#   - Model variants: symmetric, asymmetric, none (no misclassification)
#   - Stationarity: stationary (alpha derived) and free (alpha estimated)
#
# Outputs:
#   - EM-baseline/output/results/fit_{model_type}_{stat_label}.rds  (per model)
#   - EM-baseline/output/results/run_summary.csv                    (all runs)
#   - EM-baseline/output/tables/table_baseline_results.tex          (LaTeX)
#
# Usage: Rscript EM-baseline/estimate_baseline_pipeline.R
#        (from project root)
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------

library(here)
library(dplyr)

source(here::here("EM-baseline", "R", "source_all.R"))

# Estimation hyperparameters (named constants — edit here to change globally)
THETA_CAP  <- 0.999
PI_CAP     <- 0.49
N_STARTS   <- 5L
PERTURB_SD <- 0.30

# Random seed for reproducibility (multi-start initialisation and stargazer fake-lm)
RANDOM_SEED <- 1234L
set.seed(RANDOM_SEED)

# Output directories
results_dir <- here::here("EM-baseline", "output", "results")
tables_dir  <- here::here("EM-baseline", "output", "tables")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir,  recursive = TRUE, showWarnings = FALSE)

# Timestamp for this run (used in run_summary.csv)
run_ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

# ------------------------------------------------------------------------------
# 1. Ingest and prepare data
# ------------------------------------------------------------------------------

source(here::here("scripts", "ingest_data_3waves_SA.R"))  # loads df_qlfs

# Guard against factor-encoded employment status (as.integer(factor) = 1/2, not 0/1)
for (y_col in c("y1", "y2", "y3")) {
  if (is.factor(df_qlfs[[y_col]]))
    stop(sprintf(
      "Column '%s' in df_qlfs is a factor. Convert to integer 0/1 in the ingest script.",
      y_col
    ))
}
# Binary check on raw values before coercion — catches non-integer numerics (e.g. 0.5)
# that as.integer() would silently truncate to valid-looking 0L or 1L.
for (y_col in c("y1", "y2", "y3")) {
  non_na <- df_qlfs[[y_col]][!is.na(df_qlfs[[y_col]])]
  if (!all(non_na %in% c(0, 1)))
    stop(sprintf(
      "Column '%s' in df_qlfs contains non-binary values: %s",
      y_col, paste(setdiff(unique(non_na), c(0, 1)), collapse = ", ")
    ))
}

# Keep only the columns needed for the baseline model (no tenure/timegap)
df_baseline <- df_qlfs |>
  filter(!is.na(y1), !is.na(y2), !is.na(y3), !is.na(weight)) |>
  mutate(
    y1     = as.integer(y1),
    y2     = as.integer(y2),
    y3     = as.integer(y3),
    weight = as.numeric(weight)
  ) |>
  select(y1, y2, y3, weight)

stopifnot(
  all(df_baseline$y1 %in% c(0L, 1L)),
  all(df_baseline$y2 %in% c(0L, 1L)),
  all(df_baseline$y3 %in% c(0L, 1L)),
  all(df_baseline$weight > 0)
)

# Free the larger raw frame once df_baseline is validated
rm(df_qlfs)

cat(sprintf("Data loaded: N = %d observations\n", nrow(df_baseline)))

# ------------------------------------------------------------------------------
# 2. Define model configurations
# ------------------------------------------------------------------------------

configs <- list(
  list(model_type = "none",       stationary = TRUE,  label = "none_stat"),
  list(model_type = "none",       stationary = FALSE, label = "none_free"),
  list(model_type = "symmetric",  stationary = TRUE,  label = "sym_stat"),
  list(model_type = "symmetric",  stationary = FALSE, label = "sym_free"),
  list(model_type = "asymmetric", stationary = TRUE,  label = "asym_stat"),
  list(model_type = "asymmetric", stationary = FALSE, label = "asym_free")
)

# Number of random starts per model
# (N_STARTS defined in constants block above)

# ------------------------------------------------------------------------------
# 3. Estimate all models
# ------------------------------------------------------------------------------

fits       <- list()
run_rows   <- list()
start_time <- proc.time()

# Helper: perturb a probability on the logit scale and back-transform
.perturb_param <- function(p, sd = PERTURB_SD) inv_logit(logit(p) + rnorm(1, 0, sd))

# Helper: perturb and clamp a misclassification probability into (eps, PI_CAP - margin).
# inv_logit() always returns (0,1), but explicit eps lower-bound prevents a
# one-in-a-billion IEEE underflow case (inv_logit returning exactly 0.0) from
# producing a params0$pi <= 0 error downstream.
.clamp_pi <- function(p, cap = PI_CAP, margin = 0.01, eps = 1e-6) {
  min(cap - margin, max(eps, .perturb_param(p)))
}

for (cfg in configs) {
  cat(sprintf("\n--- Fitting: %s ---\n", cfg$label))
  t0 <- proc.time()

  # Multi-start: run EM from multiple random initialisations and keep the
  # solution with the highest log-likelihood
  best_fit <- NULL

  for (s in seq_len(N_STARTS)) {
    # Perturb default starting values
    p0 <- init_params(model_type = cfg$model_type, stationary = cfg$stationary)

    # Add small random perturbation (on logit scale, then back-transform)
    p0$theta0 <- .perturb_param(p0$theta0)
    p0$theta1 <- .perturb_param(p0$theta1)
    if (!cfg$stationary) p0$alpha <- .perturb_param(p0$alpha)
    if (!is.null(p0$pi))  p0$pi  <- .clamp_pi(p0$pi)
    if (!is.null(p0$pi0)) p0$pi0 <- .clamp_pi(p0$pi0)
    if (!is.null(p0$pi1)) p0$pi1 <- .clamp_pi(p0$pi1)

    fit <- tryCatch(
      em_fit_baseline(
        df         = df_baseline,
        model_type = cfg$model_type,
        stationary = cfg$stationary,
        params0    = p0,
        max_iter   = 1000L,
        tol        = 1e-8,
        theta_cap  = THETA_CAP,
        pi_cap     = PI_CAP,
        verbose    = 0L
      ),
      error = function(e) {
        warning(sprintf("Start %d for %s failed: %s", s, cfg$label, e$message))
        NULL
      }
    )

    if (!is.null(fit) && (is.null(best_fit) || fit$loglik > best_fit$loglik)) {
      best_fit <- fit
    }
  }

  if (is.null(best_fit)) {
    warning(sprintf("All starts failed for model: %s", cfg$label))
    next
  }

  elapsed <- (proc.time() - t0)[["elapsed"]]
  cat(sprintf("  loglik = %.4f | converged = %s | iters = %d | time = %.1fs\n",
              best_fit$loglik, best_fit$converged, best_fit$iterations, elapsed))

  # Store fit
  fits[[cfg$label]] <- best_fit

  # Save .rds checkpoint
  rds_path <- file.path(results_dir, sprintf("fit_%s.rds", cfg$label))
  saveRDS(best_fit, rds_path)

  # Build run summary row
  p <- best_fit$params
  run_rows[[cfg$label]] <- data.frame(
    timestamp  = run_ts,
    model_type = cfg$model_type,
    stationary = cfg$stationary,
    label      = cfg$label,
    converged  = best_fit$converged,
    iterations = best_fit$iterations,
    loglik     = best_fit$loglik,
    alpha      = p$alpha,
    theta0     = p$theta0,
    theta1     = p$theta1,
    pi         = p$pi  %||% NA_real_,
    pi0        = p$pi0 %||% NA_real_,
    pi1        = p$pi1 %||% NA_real_,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# 4. Append to run_summary.csv
# ------------------------------------------------------------------------------

run_summary <- do.call(rbind, run_rows)
csv_path    <- file.path(results_dir, "run_summary.csv")

if (file.exists(csv_path)) {
  existing <- read.csv(csv_path, stringsAsFactors = FALSE)
  # Validate columns match before appending to avoid silent corruption
  if (!identical(names(existing), names(run_summary))) {
    backup_path <- sub("\\.csv$",
                       sprintf("_backup_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")),
                       csv_path)
    file.copy(csv_path, backup_path)
    stop(sprintf(paste0(
      "run_summary.csv column schema changed \u2014 cannot append safely.\n",
      "Existing file backed up to: %s\n",
      "Reconcile column names before re-running."
    ), backup_path))
  }
  run_summary <- rbind(existing, run_summary)
  # Deduplicate: keep the most recent run per model label
  run_summary <- run_summary[!duplicated(run_summary[["label"]], fromLast = TRUE), ]
}
write.csv(run_summary, csv_path, row.names = FALSE)
cat(sprintf("\nRun summary written to %s\n", csv_path))

# ------------------------------------------------------------------------------
# 5. Likelihood ratio tests
# ------------------------------------------------------------------------------

cat("\n--- Likelihood Ratio Tests ---\n")

lr_test <- function(label_restricted, label_unrestricted, df_test,
                    description) {
  fit_r <- fits[[label_restricted]]
  fit_u <- fits[[label_unrestricted]]
  if (is.null(fit_r) || is.null(fit_u)) {
    cat(sprintf("  %-45s  [skipped: model not available]\n", description))
    return(invisible(NULL))
  }
  lr  <- 2 * (fit_u$loglik - fit_r$loglik)
  pv  <- pchisq(lr, df = df_test, lower.tail = FALSE)
  cat(sprintf("  %-45s  LR = %7.3f  df = %d  p = %.4f\n",
              description, lr, df_test, pv))
}

# Misclassification tests (stationary models)
lr_test("none_stat",  "sym_stat",  df_test = 1, "H0: pi=0  (stationary)")
lr_test("sym_stat",   "asym_stat", df_test = 1, "H0: pi0=pi1  (stationary)")
lr_test("none_free",  "sym_free",  df_test = 1, "H0: pi=0  (free alpha)")
lr_test("sym_free",   "asym_free", df_test = 1, "H0: pi0=pi1  (free alpha)")

# Stationarity tests (within each misclassification variant)
lr_test("none_stat",  "none_free",  df_test = 1, "H0: stationary  (no error)")
lr_test("sym_stat",   "sym_free",   df_test = 1, "H0: stationary  (symmetric)")
lr_test("asym_stat",  "asym_free",  df_test = 1, "H0: stationary  (asymmetric)")

# ------------------------------------------------------------------------------
# 6. Overidentification tests (TeX Eq 22)
# ------------------------------------------------------------------------------

cat("\n--- Overidentification Tests (LR vs saturated model) ---\n")

# Saturated log-likelihood: sum_abc n_abc * log(p_hat_abc)
W          <- sum(df_baseline$weight)
cell_counts <- df_baseline |>
  group_by(y1, y2, y3) |>
  summarise(n = sum(weight), .groups = "drop")
ll_sat <- sum(cell_counts$n * log(cell_counts$n / W))

for (cfg in configs) {
  fit <- fits[[cfg$label]]
  if (is.null(fit)) next
  # Degrees of freedom: free moments (7) minus free parameters
  n_params <- length(fit$params)
  df_overid <- 7L - n_params
  if (df_overid <= 0) {
    cat(sprintf("  %-20s  [saturated or over-parameterised]\n", cfg$label))
    next
  }
  lr <- 2 * (ll_sat - fit$loglik)
  pv <- pchisq(lr, df = df_overid, lower.tail = FALSE)
  cat(sprintf("  %-20s  LR = %7.3f  df = %d  p = %.4f  (n_params = %d)\n",
              cfg$label, lr, df_overid, pv, n_params))
}

# ------------------------------------------------------------------------------
# 7. LaTeX output table (stargazer)
# ------------------------------------------------------------------------------

cat("\n--- Creating LaTeX output tables ---\n")

# Build a results summary data frame for display
param_table <- run_summary[run_summary$timestamp == run_ts, ]
param_table$alpha  <- round(param_table$alpha,  4)
param_table$theta0 <- round(param_table$theta0, 4)
param_table$theta1 <- round(param_table$theta1, 4)
param_table$pi     <- round(param_table$pi,     4)
param_table$pi0    <- round(param_table$pi0,    4)
param_table$pi1    <- round(param_table$pi1,    4)
param_table$loglik <- round(param_table$loglik, 2)

# Stargazer table: treat each model as a fake lm() object populated with EM params.
# Re-seed immediately before rnorm() calls so residuals are reproducible.
make_fake_lm <- function(params_named, n_obs, seed) {
  set.seed(seed)
  shell <- lm(rnorm(n_obs) ~ 1)
  shell$coefficients <- unlist(params_named)
  names(shell$coefficients) <- names(params_named)
  set.seed(seed + 1L)
  shell$residuals <- rnorm(n_obs)
  shell
}

n_obs <- nrow(df_baseline)

# Each model gets a non-overlapping seed pair: model i uses seeds
# (RANDOM_SEED + 2L*i) and (RANDOM_SEED + 2L*i + 1L) — make_fake_lm()
# calls set.seed(seed) and set.seed(seed + 1L) internally.
model_list <- lapply(seq_along(configs), function(i) {
  cfg <- configs[[i]]
  fit <- fits[[cfg$label]]
  if (is.null(fit)) return(NULL)
  make_fake_lm(fit$params, n_obs, seed = RANDOM_SEED + 2L * i)
})
names(model_list) <- sapply(configs, `[[`, "label")
model_list        <- Filter(Negate(is.null), model_list)

# Column headers and loglik rows for stargazer
col_labels <- gsub("_", " ", names(model_list))
ll_vals    <- sapply(names(model_list), function(lbl) round(fits[[lbl]]$loglik, 2))
n_iters    <- sapply(names(model_list), function(lbl) fits[[lbl]]$iterations)
conv       <- sapply(names(model_list), function(lbl) as.integer(fits[[lbl]]$converged))

table_baseline <- stargazer::stargazer(
  model_list,
  type        = "latex",
  title       = "EM baseline model estimates (SA QLFS, 3 waves)",
  label       = "tab:em_baseline",
  column.labels = col_labels,
  covariate.labels = c(
    "alpha (initial employment)",
    "theta0 (job-finding)",
    "theta1 (employment persistence)",
    "pi (symmetric misclassification)",
    "pi0 (false positive rate)",
    "pi1 (false negative rate)"
  ),
  omit.stat   = c("adj.rsq", "ser", "f", "rsq"),
  keep.stat   = c("n"),
  add.lines   = list(
    c("Log-likelihood", as.character(ll_vals)),
    c("EM iterations",  as.character(n_iters)),
    c("Converged",      as.character(conv))
  ),
  notes       = "EM algorithm estimates. Standard errors not reported (use profile likelihood or bootstrap for inference).",
  header      = FALSE
)

tex_path <- file.path(tables_dir, "table_baseline_results.tex")
cat(table_baseline, file = tex_path, sep = "\n")
cat(sprintf("LaTeX table written to %s\n", tex_path))

# ------------------------------------------------------------------------------
# 8. Runtime summary
# ------------------------------------------------------------------------------

total_elapsed <- (proc.time() - start_time)[["elapsed"]]
cat(sprintf("\nTotal pipeline runtime: %.1f seconds\n", total_elapsed))
cat("Done.\n")
