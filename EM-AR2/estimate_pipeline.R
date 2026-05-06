# ==============================================================================
# EM-AR2: End-to-end estimation pipeline
# ==============================================================================
# Created: 2026-05-05
# This script:
#   1. Sources the EM-AR2 module
#   2. Sources the 4-wave data ingestion script
#   3. Filters and prepares the data (period1 >= 30 & period1 <= 32; normalised
#      weights)
#   4. Runs three EM variants:
#        fit_sym   — symmetric misclassification (estimate_pi = TRUE)
#        fit_nopi  — no misclassification (estimate_pi = FALSE, fixed_pi = 0)
#        fit_asym  — asymmetric misclassification (asymmetric = TRUE)
#   5. Saves each fit as a timestamped .rds in EM-AR2/output/results/
#   6. Appends one summary row per model to EM-AR2/output/results/run_summary.csv
#   7. Computes goodness-of-fit and implied transitions for each fit
#   8. Writes three LaTeX tables to EM-AR2/output/tables/
#
# Usage (from project root):
#   source("EM-AR2/estimate_pipeline.R")
# ==============================================================================

# Verify working directory -------------------------------------------------------
if (!file.exists("EM-AR2/R/source_all.R")) {
  stop(
    "estimate_pipeline.R must be sourced from the project root. ",
    "Expected to find 'EM-AR2/R/source_all.R' relative to cwd."
  )
}

# Load EM-AR2 module -------------------------------------------------------------
source("EM-AR2/R/source_all.R")

# Ingest 4-wave data -------------------------------------------------------------
if (!file.exists("data/raw/df_qlfs_A.rds")) {
  stop("Missing input data: 'data/raw/df_qlfs_A.rds'. ",
       "Obtain the QLFS data and place it at that path before running the pipeline.")
}
if (!file.exists("scripts/ingest_data_4waves_SA.R")) {
  stop("Missing ingestion script: 'scripts/ingest_data_4waves_SA.R'.")
}

library(tidyverse)
source("scripts/ingest_data_4waves_SA.R")

# Filter to comparable cohort (QLFS period 2010 Q2 – 2011 Q1) -------------------
df_raw <- df_qlfs %>%
  filter(period1 >= 30 & period1 <= 32)

# Normalise weights so that sum(weight) = N  ------------------------------------
n_obs        <- nrow(df_raw)
weight_total <- sum(df_raw$weight)
df_4w <- df_raw %>%
  mutate(weight = n_obs * weight / weight_total)

message(sprintf("EM-AR2: N = %d observations after filtering.", n_obs))

# Verify required columns -------------------------------------------------------
required_cols <- c("y1", "y2", "y3", "y4", "weight")
missing <- setdiff(required_cols, names(df_4w))
if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

# Output directories ------------------------------------------------------------
results_dir <- "EM-AR2/output/results"
tables_dir  <- "EM-AR2/output/tables"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir,  recursive = TRUE, showWarnings = FALSE)

run_id      <- format(Sys.time(), "%Y%m%d_%H%M%S")
summary_csv <- file.path(results_dir, "run_summary.csv")

# ==============================================================================
# Estimation
# ==============================================================================
message("\n=== EM-AR2: Starting estimation ===")

# Shared EM settings
EM_MAX_ITER <- 500L
EM_TOL      <- 1e-8
EM_VERBOSE  <- 1L

# 1. Symmetric misclassification -----------------------------------------------
message("\n--- Model 1: Symmetric misclassification ---")
fit_sym <- em_fit_ar2(
  df         = df_4w,
  estimate_pi = TRUE,
  max_iter   = EM_MAX_ITER,
  tol        = EM_TOL,
  verbose    = EM_VERBOSE
)
message(sprintf("  Converged: %s  |  LL = %.4f  |  Iterations: %d",
                fit_sym$converged, fit_sym$loglik, fit_sym$iterations))

# 2. No misclassification -------------------------------------------------------
message("\n--- Model 2: No misclassification ---")
fit_nopi <- em_fit_ar2(
  df          = df_4w,
  estimate_pi = FALSE,
  fixed_pi    = 0,
  max_iter    = EM_MAX_ITER,
  tol         = EM_TOL,
  verbose     = EM_VERBOSE
)
message(sprintf("  Converged: %s  |  LL = %.4f  |  Iterations: %d",
                fit_nopi$converged, fit_nopi$loglik, fit_nopi$iterations))

# 3. Asymmetric misclassification -----------------------------------------------
message("\n--- Model 3: Asymmetric misclassification ---")
fit_asym <- em_fit_ar2(
  df         = df_4w,
  asymmetric = TRUE,
  max_iter   = EM_MAX_ITER,
  tol        = EM_TOL,
  verbose    = EM_VERBOSE
)
message(sprintf("  Converged: %s  |  LL = %.4f  |  Iterations: %d",
                fit_asym$converged, fit_asym$loglik, fit_asym$iterations))

# ==============================================================================
# Save raw fits
# ==============================================================================
fits <- list(sym = fit_sym, nopi = fit_nopi, asym = fit_asym)
for (nm in names(fits)) {
  rds_path <- file.path(results_dir, sprintf("em_ar2_%s_%s.rds", nm, run_id))
  saveRDS(fits[[nm]], rds_path)
  message(sprintf("Saved: %s", rds_path))
}

# ==============================================================================
# Goodness-of-fit & implied transitions
# ==============================================================================
gof_sym  <- goodness_of_fit_ar2(df_4w, fit_sym$params)
gof_nopi <- goodness_of_fit_ar2(df_4w, fit_nopi$params)
gof_asym <- goodness_of_fit_ar2(df_4w, fit_asym$params)

trans_sym  <- implied_transitions_ar2(fit_sym$params)
trans_nopi <- implied_transitions_ar2(fit_nopi$params)
trans_asym <- implied_transitions_ar2(fit_asym$params)

# ==============================================================================
# Append to run_summary.csv
# ==============================================================================
.summary_row <- function(nm, fit, gof) {
  p <- fit$params
  data.frame(
    run_id       = run_id,
    model        = nm,
    converged    = fit$converged,
    iterations   = fit$iterations,
    loglik       = fit$loglik,
    theta0       = p$theta0,
    theta01      = p$theta01,
    theta1       = p$theta1,
    theta10      = p$theta10,
    pi           = if (!is.null(p$pi))  p$pi  else NA_real_,
    pi0          = if (!is.null(p$pi0)) p$pi0 else NA_real_,
    pi1          = if (!is.null(p$pi1)) p$pi1 else NA_real_,
    ssr          = gof$ssr,
    n_obs        = n_obs,
    stringsAsFactors = FALSE
  )
}

summary_rows <- rbind(
  .summary_row("sym",  fit_sym,  gof_sym),
  .summary_row("nopi", fit_nopi, gof_nopi),
  .summary_row("asym", fit_asym, gof_asym)
)

if (file.exists(summary_csv)) {
  write.table(summary_rows, summary_csv,
              append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
} else {
  write.csv(summary_rows, summary_csv, row.names = FALSE)
}
message(sprintf("\nAppended 3 rows to %s", summary_csv))

# ==============================================================================
# LaTeX tables
# ==============================================================================

# Helper: build a dummy lm model for stargazer.
# The residuals and intercept are random but seeded to guarantee reproducible
# LaTeX output across pipeline runs (stargazer uses residuals for some stats).
.make_lm_for_stargazer <- function(coef_vec, n_obs_val) {
  set.seed(20260505L)  # fixed seed: reproducible table output
  dummy_lm <- lm(rnorm(2) ~ 1)
  dummy_lm$coefficients  <- coef_vec
  dummy_lm$residuals     <- rnorm(n_obs_val)
  dummy_lm
}

# --- Table 1: AR(2) EM parameter estimates ------------------------------------
# Parameters: theta0, theta01, theta1, theta10, pi (or pi0/pi1)
.extract_params_vec <- function(fit) {
  p <- fit$params
  c(
    theta0  = p$theta0,
    theta01 = p$theta01,
    theta1  = p$theta1,
    theta10 = p$theta10,
    pi      = if (!is.null(p$pi))  p$pi  else NA_real_,
    pi0     = if (!is.null(p$pi0)) p$pi0 else NA_real_,
    pi1     = if (!is.null(p$pi1)) p$pi1 else NA_real_
  )
}

params_sym  <- .extract_params_vec(fit_sym)
params_nopi <- .extract_params_vec(fit_nopi)
params_asym <- .extract_params_vec(fit_asym)

lm_sym  <- .make_lm_for_stargazer(params_sym,  n_obs)
lm_nopi <- .make_lm_for_stargazer(params_nopi, n_obs)
lm_asym <- .make_lm_for_stargazer(params_asym, n_obs)

params_table_path <- file.path(tables_dir, "table_ar2_em_params.tex")
stargazer::stargazer(
  lm_sym, lm_nopi, lm_asym,
  type       = "latex",
  out        = params_table_path,
  se         = list(rep(NA, length(params_sym)),
                    rep(NA, length(params_nopi)),
                    rep(NA, length(params_asym))),
  keep.stat  = c("n"),
  column.labels = c("AR(2) Sym ME", "AR(2) No ME", "AR(2) Asym ME"),
  covariate.labels = c(
    "Entry rate ($\\theta_0$)",
    "Surplus entry rate ($\\theta_{01}$)",
    "Exit rate ($\\theta_1$)",
    "Surplus exit rate ($\\theta_{10}$)",
    "Misclassification ($\\pi$)",
    "Misclassification from 0 ($\\pi_0$)",
    "Misclassification from 1 ($\\pi_1$)"
  ),
  dep.var.labels = "Employment",
  add.lines = list(
    c("Log-likelihood",
      round(fit_sym$loglik, 1),
      round(fit_nopi$loglik, 1),
      round(fit_asym$loglik, 1)),
    c("Converged",
      as.character(fit_sym$converged),
      as.character(fit_nopi$converged),
      as.character(fit_asym$converged))
  ),
  title  = "AR(2) EM Estimates",
  label  = "tab:ar2_em_params",
  header = FALSE
)
message(sprintf("Saved: %s", params_table_path))

# --- Table 2: Goodness-of-fit --------------------------------------------------
gof_combined <- rbind(
  cbind(model = "Sym ME",  gof_sym$table),
  cbind(model = "No ME",   gof_nopi$table),
  cbind(model = "Asym ME", gof_asym$table)
)

gof_table_path <- file.path(tables_dir, "table_ar2_em_gof.tex")

# Write a simple LaTeX table (no stargazer wrapper needed for cell-level data)
.write_gof_latex <- function(gof_tbl, ssr_sym, ssr_nopi, ssr_asym, path) {
  lines <- c(
    "\\begin{table}[ht]",
    "\\centering",
    "\\caption{AR(2) EM: Goodness-of-Fit by Cell}",
    "\\label{tab:ar2_em_gof}",
    "\\begin{tabular}{lrrrr}",
    "\\hline",
    "Cell & Observed & Sym ME & No ME & Asym ME \\\\",
    "\\hline"
  )

  cells <- unique(interaction(gof_sym$table$y1, gof_sym$table$y2,
                               gof_sym$table$y3, gof_sym$table$y4, sep=""))
  sym_tbl  <- gof_sym$table
  nopi_tbl <- gof_nopi$table
  asym_tbl <- gof_asym$table

  for (i in seq_len(nrow(sym_tbl))) {
    cell_label <- paste0(sym_tbl$y1[i], sym_tbl$y2[i],
                          sym_tbl$y3[i], sym_tbl$y4[i])
    obs_p   <- round(sym_tbl$empirical[i], 4)
    mod_sym  <- round(sym_tbl$model_prob[i], 4)
    mod_nopi <- round(nopi_tbl$model_prob[i], 4)
    mod_asym <- round(asym_tbl$model_prob[i], 4)
    lines <- c(lines,
               sprintf("%s & %.4f & %.4f & %.4f & %.4f \\\\",
                       cell_label, obs_p, mod_sym, mod_nopi, mod_asym))
  }
  lines <- c(lines,
             "\\hline",
             sprintf("SSR & & %.6f & %.6f & %.6f \\\\",
                     ssr_sym, ssr_nopi, ssr_asym),
             "\\hline",
             "\\end{tabular}",
             "\\end{table}")
  writeLines(lines, path)
}

.write_gof_latex(gof_combined, gof_sym$ssr, gof_nopi$ssr, gof_asym$ssr,
                 gof_table_path)
message(sprintf("Saved: %s", gof_table_path))

# --- Table 3: Implied transitions ----------------------------------------------
trans_table_path <- file.path(tables_dir, "table_ar2_em_transitions.tex")

.write_trans_latex <- function(t_sym, t_nopi, t_asym, path) {
  row_labels <- c(
    "Prob employed in w4, given not in w1",
    "Ever employed, given not in w1",
    "Prob not employed in w4, given employed in w1",
    "Ever not employed, given employed in w1"
  )
  fields <- c("p_emp_w4_from0", "p_ever_emp_from0",
              "p_nonemp_w4_from1", "p_ever_nonemp_from1")

  lines <- c(
    "\\begin{table}[ht]",
    "\\centering",
    "\\caption{AR(2) EM: Implied Transition Probabilities (Stationary Distribution)}",
    "\\label{tab:ar2_em_transitions}",
    "\\begin{tabular}{lrrr}",
    "\\hline",
    "Quantity & Sym ME & No ME & Asym ME \\\\",
    "\\hline"
  )
  for (i in seq_along(fields)) {
    f <- fields[i]
    lines <- c(lines,
               sprintf("%s & %.4f & %.4f & %.4f \\\\",
                       row_labels[i], t_sym[f], t_nopi[f], t_asym[f]))
  }
  lines <- c(lines, "\\hline", "\\end{tabular}", "\\end{table}")
  writeLines(lines, path)
}

.write_trans_latex(trans_sym, trans_nopi, trans_asym, trans_table_path)
message(sprintf("Saved: %s", trans_table_path))

# ==============================================================================
# Summary print
# ==============================================================================
message("\n=== EM-AR2: Pipeline complete ===")
message(sprintf("Run ID: %s", run_id))
message(sprintf("N obs:  %d", n_obs))
message("\nParameter estimates:")
print(summary_rows[, c("model","theta0","theta01","theta1","theta10","pi","pi0","pi1","loglik","converged")])
message("\nImplied transitions (Sym ME):")
print(round(trans_sym, 4))
message("\nImplied transitions (No ME):")
print(round(trans_nopi, 4))
message("\nImplied transitions (Asym ME):")
print(round(trans_asym, 4))
