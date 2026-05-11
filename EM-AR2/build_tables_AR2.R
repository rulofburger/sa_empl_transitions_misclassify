# ==============================================================================
# EM-AR2: Publication-ready LaTeX tables
# Created: 2026-05-07
#
# Reads AR(2) EM point estimates and bootstrap results, then writes two
# publication-ready LaTeX tables:
#
#   Table 1: AR(2) parameter estimates (raw scale, 4 d.p.) with bootstrap SEs
#            and significance stars testing H0: θ = 0.
#            Output: EM-AR2/output/tables/table_AR2_params.tex
#
#   Table 2: Implied transition probabilities p_{jk} and employment rate,
#            expressed as percentages (2 d.p.) with bootstrap SEs.
#            Output: EM-AR2/output/tables/table_AR2_implied.tex
#
# Columns: (1) No misclassification, (2) Symmetric π, (3) Asymmetric (π₀, π₁)
#
# The four AR(2) transition probabilities follow paper eq. 4:
#   p_{00} = θ₀          (entry after two periods not employed)
#   p_{10} = θ₀ + θ_{01} (entry after employed then not employed)
#   p_{01} = 1−θ₁−θ_{10} (retention after not employed then employed)
#   p_{11} = 1−θ₁         (retention after two employed periods)
#
# Bootstrap SEs: 200 nonparametric resamples via bootstrap_pipeline_AR2.R
#
# Usage (from project root):
#   Rscript EM-AR2/build_tables_AR2.R
#   # or from R console:
#   source("EM-AR2/build_tables_AR2.R")
#
# Prerequisites:
#   bootstrap_pipeline_AR2.R must have been run first.
# ==============================================================================

# Verify working directory -------------------------------------------------------
if (!file.exists("EM-AR2/R/source_all.R")) {
  stop(
    "build_tables_AR2.R must be sourced from the project root. ",
    "Expected to find 'EM-AR2/R/source_all.R' relative to cwd."
  )
}

# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------

library(tidyverse)   # needed by ingest_data_4waves_SA.R

source("EM-baseline/R/source_all.R")  # bootstrap_resample, summarise_bootstrap
source("EM-AR2/R/source_all.R")       # em_fit_ar2, implied_ar2, .find_latest_fit

boot_dir   <- "EM-AR2/output/results/bootstrap"
tables_dir <- "EM-AR2/output/tables"
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Detect B from existing bootstrap files
.detect_B_ar2 <- function(boot_dir) {
  fls <- list.files(boot_dir, pattern = "boot_.*_B[0-9]+\\.rds$", full.names = FALSE)
  if (length(fls) == 0L) {
    warning(".detect_B_ar2: no bootstrap files found; defaulting to B=200.")
    return(200L)
  }
  B_str <- sub(".*_B([0-9]+)\\.rds$", "\\1", fls[[1L]])
  B_val <- suppressWarnings(as.integer(B_str))
  if (is.na(B_val))
    stop(".detect_B_ar2: failed to parse B from '", fls[[1L]], "'.")
  B_val
}
B <- .detect_B_ar2(boot_dir)
cat(sprintf("B detected from bootstrap files: %d\n", B))

# Load bootstrap objects now (needed for N_AR2 and later for SEs) ---------
# This avoids re-ingesting the full QLFS dataset just to count rows.
# N_AR2 is stored in each boot object by bootstrap_pipeline_AR2.R ($n_obs).
# If no boot objects exist yet, fall back to ingesting the data.
ar2_labels     <- c("nopi", "sym", "asym")
ar2_model_types <- c(nopi = "none", sym = "symmetric", asym = "asymmetric")
results_dir    <- "EM-AR2/output/results"

.load_boot_AR2_early <- function(key) {
  path <- file.path(boot_dir, sprintf("boot_%s_B%d.rds", key, B))
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}
ar2_boots_early <- lapply(ar2_labels, .load_boot_AR2_early)
not_null_idx <- which(!vapply(ar2_boots_early, is.null, logical(1L)))

if (length(not_null_idx) > 0L && !is.null(ar2_boots_early[[not_null_idx[[1L]]]]$n_obs)) {
  N_AR2 <- ar2_boots_early[[not_null_idx[[1L]]]]$n_obs
  cat(sprintf("N_AR2 = %s (from boot object)\n",
              formatC(N_AR2, format = "d", big.mark = ",")))
} else {
  # Fallback: ingest data to count rows
  if (!file.exists("scripts/ingest_data_4waves_SA.R"))
    stop("Missing ingestion script: 'scripts/ingest_data_4waves_SA.R'.")
  source("scripts/ingest_data_4waves_SA.R")
  if (!exists("df_qlfs") || !is.data.frame(df_qlfs))
    stop("ingest_data_4waves_SA.R did not produce a 'df_qlfs' data frame.")
  N_AR2 <- sum(df_qlfs$period1 >= 30 & df_qlfs$period1 <= 32)
  rm(df_qlfs)
  cat(sprintf("N_AR2 = %s (from ingestion fallback)\n",
              formatC(N_AR2, format = "d", big.mark = ",")))
}

# ------------------------------------------------------------------------------
# 1. Table formatting helpers
# (Mirror of build_tables.R helpers, standalone to avoid coupling)
# ------------------------------------------------------------------------------

`%+%` <- function(a, b) paste0(a, b)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Significance threshold constants (standard normal quantiles, two-sided)
.CRIT_p01 <- qnorm(0.995)  # ≈ 2.576  (p < 0.01, two-sided)
.CRIT_p05 <- qnorm(0.975)  # ≈ 1.960  (p < 0.05, two-sided)
.CRIT_p10 <- qnorm(0.95)   # ≈ 1.645  (p < 0.10, two-sided)

#' Significance stars based on z = |est / se|
.stars <- function(est, se) {
  if (is.na(est) || is.na(se) || se <= 0) return("")
  z <- abs(est / se)
  if (z > .CRIT_p01) return("$^{**}$")
  if (z > .CRIT_p05) return("$^{*}$")
  if (z > .CRIT_p10) return("$^{.}$")
  ""
}

#' Scale, format, and attach significance stars.
.fmt_estimate <- function(est, se, factor = 1, digits = 4L) {
  if (is.na(est)) return(c("---", "(---)"))
  star    <- .stars(est, se)
  est_str <- paste0(formatC(est * factor, format = "f", digits = digits), star)
  se_str  <- if (!is.na(se)) {
    sprintf("(%s)", formatC(se * factor, format = "f", digits = digits))
  } else {
    "(---)"
  }
  c(est_str, se_str)
}

#' Format a raw parameter (4 d.p.)
.fmt_param <- function(est, se, digits = 4L) {
  .fmt_estimate(est, se, factor = 1, digits = digits)
}

#' Format as percentage (×100, 2 d.p.)
.fmt_pct <- function(est, se, digits = 2L) {
  .fmt_estimate(est, se, factor = 100, digits = digits)
}

#' Format log-likelihood in millions
.fmt_ll <- function(ll) {
  if (is.null(ll) || is.na(ll)) return("---")
  sprintf("%.1fM", ll / 1e6)
}

#' Format sample size with comma separator
.fmt_n <- function(n) {
  formatC(n, format = "d", big.mark = ",")
}

#' Write a LaTeX booktabs table to a file.
.write_latex_table <- function(col_headers, row_data, caption, label, path,
                               col_align = NULL, note = NULL,
                               sub_headers = NULL) {
  n_cols <- length(col_headers)
  if (is.null(col_align)) {
    col_align <- paste0("l", paste(rep("c", n_cols - 1), collapse = ""))
  }

  buf <- vector("list", 300L)
  i   <- 0L
  push <- function(x) { i <<- i + 1L; buf[[i]] <<- x }

  push("\\begin{table}[htbp]")
  push("\\centering")
  push(sprintf("\\caption{%s}", caption))
  push(sprintf("\\label{%s}", label))
  push(sprintf("\\begin{tabular}{%s}", col_align))
  push("\\toprule")
  push(paste(col_headers, collapse = " & ") %+% " \\\\")

  if (!is.null(sub_headers)) for (sh in sub_headers) push(sh)

  push("\\midrule")

  for (block in row_data) {
    if (!is.null(block$header))
      push(sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\", n_cols, block$header))
    for (row in block$rows) {
      push(paste(row[["est"]], collapse = " & ") %+% " \\\\")
      se_row <- row[["se"]]
      if (length(se_row) > 0 && !all(se_row == ""))
        push(paste(se_row, collapse = " & ") %+% " \\\\")
    }
  }

  push("\\bottomrule")

  if (!is.null(note))
    push(sprintf("\\multicolumn{%d}{l}{\\footnotesize \\textit{Note:} %s} \\\\",
                 n_cols, note))

  push("\\end{tabular}")
  push("\\end{table}")

  cat(paste(unlist(buf[seq_len(i)]), collapse = "\n"), "\n", file = path)
  cat(sprintf("Written: %s\n", path))
}

# ------------------------------------------------------------------------------
# 2. Load AR(2) point estimates and bootstrap summaries
# ------------------------------------------------------------------------------

#' Load and validate a bootstrap summary .rds for a single AR(2) model.
.load_boot_AR2 <- function(key) {
  path <- file.path(boot_dir, sprintf("boot_%s_B%d.rds", key, B))
  if (!file.exists(path)) {
    message(sprintf(".load_boot_AR2: '%s' not found — SEs will be NA.", basename(path)))
    return(NULL)
  }
  obj <- readRDS(path)
  if (!is.list(obj) || !is.data.frame(obj$summary))
    stop(sprintf(".load_boot_AR2: '%s' missing required $summary data frame.",
                 basename(path)))
  obj
}

#' Pre-index bootstrap summary as a named SE vector for O(1) lookup.
.se_map_ar2 <- function(boot_obj) {
  if (is.null(boot_obj) || is.null(boot_obj$summary)) return(NULL)
  setNames(boot_obj$summary$se, boot_obj$summary$quantity)
}

#' Retrieve SE for a quantity name from a pre-indexed se_map.
.get_se_ar2 <- function(se_map, quantity) {
  if (is.null(se_map)) return(NA_real_)
  v <- se_map[[quantity]]
  if (is.null(v)) NA_real_ else v
}

# ar2_labels, ar2_model_types, results_dir declared above (near N_AR2 block).

# Load point-estimate fits (latest timestamped .rds per model)
ar2_fits <- lapply(ar2_labels, function(lbl) {
  path <- tryCatch(.find_latest_fit(lbl, results_dir), error = function(e) NULL)
  if (is.null(path)) {
    message(sprintf("No point estimate found for '%s'. Run estimate_pipeline.R.", lbl))
    return(NULL)
  }
  obj <- readRDS(path)
  if (!is.list(obj) || is.null(obj$params) || is.null(obj$loglik))
    stop(sprintf("Fit file for '%s' missing $params or $loglik.", lbl))
  obj
})
names(ar2_fits) <- ar2_labels

# Compute point-estimate implied quantities
ar2_implied <- lapply(ar2_labels, function(lbl) {
  fit <- ar2_fits[[lbl]]
  if (is.null(fit)) return(NULL)
  tryCatch(
    implied_ar2(fit$params, ar2_model_types[[lbl]]),
    error = function(e) NULL
  )
})
names(ar2_implied) <- ar2_labels

# Load bootstrap summaries
ar2_boots <- lapply(ar2_labels, .load_boot_AR2)
names(ar2_boots) <- ar2_labels
ar2_se    <- lapply(ar2_boots, .se_map_ar2)
names(ar2_se) <- ar2_labels

# Helpers to extract estimate and SE for a given quantity
.ar2_est <- function(lbl, quantity) {
  imp <- ar2_implied[[lbl]]
  if (is.null(imp)) return(NA_real_)
  v <- imp[[quantity]]
  if (is.null(v)) NA_real_ else v
}

.ar2_param_est <- function(lbl, quantity) {
  fit <- ar2_fits[[lbl]]
  if (is.null(fit)) return(NA_real_)
  v <- fit$params[[quantity]]
  if (is.null(v)) NA_real_ else v
}

# Column headers (1-based numbering; first col is row labels)
ar2_col_headers  <- c("", "(1)", "(2)", "(3)")
ar2_sub_headers  <- c(
  paste(c("", "No miscl.", "Symmetric $\\pi$", "Asymmetric $(\\pi_0, \\pi_1)$"),
        collapse = " & ") %+% " \\\\"
)

# Table note
ar2_note <- paste(
  "Bootstrap standard errors (B~=~", B, ") in parentheses.",
  " Stars: $^{.}$ p$<$0.10, $^{*}$ p$<$0.05, $^{**}$ p$<$0.01",
  " (two-sided z-test, H$_0$: estimate $=$ 0).",
  sep = ""
)

# ------------------------------------------------------------------------------
# 3. Table 1 — Parameter estimates (raw scale, 4 d.p.)
# ------------------------------------------------------------------------------

cat("\n--- Building Table 1: AR(2) parameter estimates ---\n")

.ar2_param_row <- function(param_name, row_label) {
  cells <- lapply(ar2_labels, function(lbl) {
    est <- .ar2_param_est(lbl, param_name)
    se  <- .get_se_ar2(ar2_se[[lbl]], param_name)
    .fmt_param(est, se)
  })
  list(
    est = c(row_label,    vapply(cells, `[[`, character(1L), 1L)),
    se  = c("",           vapply(cells, `[[`, character(1L), 2L))
  )
}

.ar2_pi_row <- function(pi_name, row_label) {
  # Only the relevant model shows a value; others show dashes
  cells <- lapply(ar2_labels, function(lbl) {
    est <- .ar2_param_est(lbl, pi_name)
    se  <- .get_se_ar2(ar2_se[[lbl]], pi_name)
    if (is.na(est)) c("---", "(---)") else .fmt_param(est, se)
  })
  list(
    est = c(row_label,    vapply(cells, `[[`, character(1L), 1L)),
    se  = c("",           vapply(cells, `[[`, character(1L), 2L))
  )
}

params_rows <- list(
  list(
    header = NULL,
    rows   = list(
      .ar2_param_row("theta0",  "$\\theta_0$"),
      .ar2_param_row("theta01", "$\\theta_{01}$"),
      .ar2_param_row("theta1",  "$\\theta_1$"),
      .ar2_param_row("theta10", "$\\theta_{10}$")
    )
  ),
  list(
    header = "Misclassification",
    rows   = list(
      .ar2_pi_row("pi",  "$\\pi$"),
      .ar2_pi_row("pi0", "$\\pi_0$"),
      .ar2_pi_row("pi1", "$\\pi_1$")
    )
  ),
  list(
    header = NULL,
    rows   = list(
      # Log-likelihood
      list(
        est = c("Log-likelihood",
                vapply(ar2_labels, function(lbl) {
                  fit <- ar2_fits[[lbl]]
                  if (is.null(fit)) "---" else .fmt_ll(fit$loglik)
                }, character(1L))),
        se  = c("", "", "", "")
      ),
      # Sample size
      list(
        est = c("$N$",
                .fmt_n(N_AR2), .fmt_n(N_AR2), .fmt_n(N_AR2)),
        se  = c("", "", "", "")
      )
    )
  )
)

.write_latex_table(
  col_headers = ar2_col_headers,
  row_data    = params_rows,
  caption     = "AR(2) EM model: parameter estimates",
  label       = "tab:ar2_params",
  path        = file.path(tables_dir, "table_AR2_params.tex"),
  col_align   = "lccc",
  note        = ar2_note,
  sub_headers = ar2_sub_headers
)

# ------------------------------------------------------------------------------
# 4. Table 2 — Implied transition probabilities (as %, 2 d.p.)
# ------------------------------------------------------------------------------

cat("\n--- Building Table 2: AR(2) implied transition probabilities ---\n")

.ar2_implied_row <- function(quantity, row_label) {
  cells <- lapply(ar2_labels, function(lbl) {
    est <- .ar2_est(lbl, quantity)
    se  <- .get_se_ar2(ar2_se[[lbl]], quantity)
    .fmt_pct(est, se)
  })
  list(
    est = c(row_label, vapply(cells, `[[`, character(1L), 1L)),
    se  = c("",        vapply(cells, `[[`, character(1L), 2L))
  )
}

.ar2_implied_pi_row <- function(quantity, row_label) {
  # Dash for models where the quantity is not estimated
  cells <- lapply(ar2_labels, function(lbl) {
    est <- .ar2_est(lbl, quantity)
    if (is.na(est)) {
      return(c("---", "(---)"))
    }
    se <- .get_se_ar2(ar2_se[[lbl]], quantity)
    .fmt_pct(est, se)
  })
  list(
    est = c(row_label, vapply(cells, `[[`, character(1L), 1L)),
    se  = c("",        vapply(cells, `[[`, character(1L), 2L))
  )
}

implied_rows <- list(
  list(
    header = "Transition probabilities",
    rows   = list(
      .ar2_implied_row("p_00", "$p_{00}$ (\\%)"),
      .ar2_implied_row("p_10", "$p_{10}$ (\\%)"),
      .ar2_implied_row("p_01", "$p_{01}$ (\\%)"),
      .ar2_implied_row("p_11", "$p_{11}$ (\\%)")
    )
  ),
  list(
    header = "Steady state",
    rows   = list(
      .ar2_implied_row("employment_rate", "Employment rate $\\alpha^*$ (\\%)")
    )
  ),
  list(
    header = "Misclassification",
    rows   = list(
      .ar2_implied_pi_row("pi",  "$\\pi$ (\\%)"),
      .ar2_implied_pi_row("pi0", "$\\pi_0$ (\\%)"),
      .ar2_implied_pi_row("pi1", "$\\pi_1$ (\\%)")
    )
  ),
  list(
    header = NULL,
    rows   = list(
      list(
        est = c("Log-likelihood",
                vapply(ar2_labels, function(lbl) {
                  fit <- ar2_fits[[lbl]]
                  if (is.null(fit)) "---" else .fmt_ll(fit$loglik)
                }, character(1L))),
        se  = c("", "", "", "")
      ),
      list(
        est = c("$N$",
                .fmt_n(N_AR2), .fmt_n(N_AR2), .fmt_n(N_AR2)),
        se  = c("", "", "", "")
      )
    )
  )
)

.write_latex_table(
  col_headers = ar2_col_headers,
  row_data    = implied_rows,
  caption     = paste0(
    "AR(2) EM model: implied transition probabilities (\\%)",
    "\\newline",
    "\\small $p_{jk}$ = Pr(employed in period $t$ $|$ status in $t-2$ was $j$, ",
    "status in $t-1$ was $k$)"
  ),
  label       = "tab:ar2_implied",
  path        = file.path(tables_dir, "table_AR2_implied.tex"),
  col_align   = "lccc",
  note        = ar2_note,
  sub_headers = ar2_sub_headers
)

cat(sprintf("\nDone. Tables written to: %s\n", tables_dir))
