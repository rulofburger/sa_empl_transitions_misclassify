# ==============================================================================
# Table builder for the EM tenure contamination (epsilon) model
# Created: 2026-05-07
#
# Reads point estimates (fit_eps_*.rds) and bootstrap results
# (boot_eps_*_B200.rds), then produces three publication-ready LaTeX tables:
#
#   Table 1: Raw parameter estimates          (table_eps_params.tex)
#   Table 2: Implied probabilities (%)        (table_eps_implied.tex)
#   Table 3: Tenure interpretation            (table_eps_tenure.tex)
#
# Four model columns:
#   (1) Free (non-stationary)     — 7 free parameters
#   (2) Stationary                — 6 free parameters
#   (3) CTMC-linked               — 5 free parameters
#   (4) Stationary + linked       — 4 free parameters (preferred)
#
# Formatting conventions (matching build_tables.R):
#   - Significance stars: . (10%), * (5%), ** (1%), z = |est/SE|
#   - Percentages: ×100, 2 d.p.
#   - Duration quantities: 2 d.p. (years) or 1 d.p. (months)
#   - Log-likelihood: in millions (e.g., -1745.4M)
#   - Constrained parameters: shown as "---"
#   - Bootstrap SE in parentheses below each estimate
#
# Prerequisites:
#   bootstrap_pipeline_tenure_contamination.R must have been run first.
#
# Usage (from project root):
#   Rscript EM-tenure/build_tables_tenure_contamination.R
# ==============================================================================

library(here)

source(here::here("EM-tenure", "R", "source_all.R"))
source(here::here("EM-tenure", "R", "implied_quantities_tenure_contamination.R"))

# Output directory
tables_dir <- here::here("EM-tenure", "output", "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Point-estimate results directory
results_dir <- here::here("output", "results")
boot_dir    <- here::here("EM-tenure", "output", "results", "bootstrap")

# ==============================================================================
# 0. Data: sample size N
# ==============================================================================

# N is read from the first available fit file (all variants share the same data)
.detect_N <- function() {
  # Load first available fit to get N from gamma dimensions
  for (prefix in c("fit_eps_", "fit_eps_stationary_",
                   "fit_eps_linked_", "fit_eps_stationary_linked_")) {
    pattern <- paste0("^", prefix, "\\d{8}_\\d{6}\\.rds$")
    f <- sort(list.files(results_dir, pattern = pattern, full.names = TRUE))
    if (length(f) > 0L) {
      obj <- readRDS(tail(f, 1L))
      if (!is.null(obj$gamma)) return(nrow(obj$gamma))
    }
  }
  NA_integer_
}
N <- .detect_N()
cat(sprintf("N = %s\n", if (is.na(N)) "unknown" else formatC(N, big.mark = ",")))

# ==============================================================================
# 1. Table formatting helpers (mirrors build_tables.R)
# ==============================================================================

`%+%` <- function(a, b) paste0(a, b)

.CRIT_p01 <- qnorm(0.995)
.CRIT_p05 <- qnorm(0.975)
.CRIT_p10 <- qnorm(0.95)

#' Apply significance stars based on z-statistic.
#' @keywords internal
.stars <- function(est, se) {
  if (is.na(est) || is.na(se) || se <= 0) return("")
  z <- abs(est / se)
  if (z > .CRIT_p01) return("$^{**}$")
  if (z > .CRIT_p05) return("$^{*}$")
  if (z > .CRIT_p10) return("$^{.}$")
  ""
}

#' Format a point estimate + SE pair for a LaTeX table cell.
#' @keywords internal
.fmt_estimate <- function(est, se, factor = 1, digits = 4L) {
  if (is.null(est) || is.na(est)) return(c("---", "(---)"))
  star    <- .stars(est, se)
  est_str <- paste0(formatC(est * factor, format = "f", digits = digits), star)
  se_str  <- if (!is.null(se) && !is.na(se)) {
    sprintf("(%s)", formatC(se * factor, format = "f", digits = digits))
  } else {
    "(---)"
  }
  c(est_str, se_str)
}

#' Format a raw parameter estimate (no scaling). @keywords internal
.fmt_param <- function(est, se, digits = 4L) .fmt_estimate(est, se, 1,   digits)

#' Format a probability as a percentage (×100). @keywords internal
.fmt_pct   <- function(est, se, digits = 2L) .fmt_estimate(est, se, 100, digits)

#' Format a duration quantity (custom digits). @keywords internal
.fmt_dur <- function(est, se, digits = 2L) {
  if (is.null(est) || is.infinite(est) || is.na(est)) return(c("Inf", "(---)"))
  .fmt_estimate(est, se, 1, digits)
}

#' Format a duration in months (1 d.p.). @keywords internal
.fmt_months <- function(est, se) .fmt_dur(est, se, digits = 1L)

#' Format a log-likelihood in millions. @keywords internal
.fmt_ll <- function(ll) {
  if (is.null(ll) || is.na(ll)) return("---")
  sprintf("%.1fM", ll / 1e6)
}

#' Format a sample size with comma separator. @keywords internal
.fmt_n <- function(n) formatC(n, format = "d", big.mark = ",")

#' Write a LaTeX booktabs table to file (identical logic to build_tables.R)
.write_latex_table <- function(col_headers, row_data, caption, label, path,
                               col_align = NULL, note = NULL,
                               sub_headers = NULL) {
  n_cols <- length(col_headers)
  if (is.null(col_align))
    col_align <- paste0("l", paste(rep("c", n_cols - 1L), collapse = ""))

  buf  <- vector("list", 300L)
  i    <- 0L
  push <- function(x) { i <<- i + 1L; buf[[i]] <<- x }

  push("\\begin{table}[htbp]")
  push("\\centering")
  push(sprintf("\\caption{%s}", caption))
  push(sprintf("\\label{%s}", label))
  push(sprintf("\\begin{tabular}{%s}", col_align))
  push("\\toprule")
  push(paste(col_headers, collapse = " & "))
  push("\\\\")
  if (!is.null(sub_headers)) for (sh in sub_headers) push(sh)
  push("\\midrule")

  for (block in row_data) {
    if (!is.null(block$header))
      push(sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\",
                   n_cols, block$header))
    for (row in block$rows) {
      push(paste(row[["est"]], collapse = " & "))
      se_row <- row[["se"]]
      if (length(se_row) > 0L && !all(se_row == ""))
        push(paste(se_row, collapse = " & "))
      push("\\\\")
    }
  }

  push("\\bottomrule")
  if (!is.null(note))
    push(sprintf(
      "\\multicolumn{%d}{l}{\\footnotesize \\textit{Note:} %s} \\\\",
      n_cols, note
    ))
  push("\\end{tabular}")
  push("\\end{table}")

  cat(paste(unlist(buf[seq_len(i)]), collapse = "\n"), "\n", file = path)
  cat(sprintf("Written: %s\n", path))
}

# ==============================================================================
# 2. Load fits and bootstraps
# ==============================================================================

eps_labels <- c("eps_free", "eps_stationary", "eps_linked", "eps_stationary_linked")

# Prefixes for auto-discovery of timestamped fit files
.fit_prefix <- c(
  eps_free              = "fit_eps_",
  eps_stationary        = "fit_eps_stationary_",
  eps_linked            = "fit_eps_linked_",
  eps_stationary_linked = "fit_eps_stationary_linked_"
)

# Which parameters are constrained in each variant
.is_stationary <- c(
  eps_free              = FALSE,
  eps_stationary        = TRUE,
  eps_linked            = FALSE,
  eps_stationary_linked = TRUE
)
.is_linked <- c(
  eps_free              = FALSE,
  eps_stationary        = FALSE,
  eps_linked            = TRUE,
  eps_stationary_linked = TRUE
)

#' Load the most recent fit file for a given prefix (timestamped naming).
.load_latest_eps_fit <- function(prefix) {
  pattern <- paste0("^", prefix, "\\d{8}_\\d{6}\\.rds$")
  files   <- sort(list.files(results_dir, pattern = pattern, full.names = TRUE))
  if (length(files) == 0L) {
    warning(sprintf(".load_latest_eps_fit: no file for prefix '%s'", prefix))
    return(NULL)
  }
  obj <- readRDS(tail(files, 1L))
  if (!is.list(obj) || is.null(obj$params) || is.null(obj$loglik))
    stop(sprintf(".load_latest_eps_fit: file missing $params or $loglik"))
  obj
}

#' Load bootstrap .rds and extract SE map.
#'
#' Discovers the file by glob pattern (B and seed are encoded in the filename).
#' Selects the most recently modified file to avoid lexicographic ordering errors
#' when multiple files with different B or seed values coexist.
#' Reads B and master_seed from the object so the table builder does not need
#' to hard-code them.
#' @keywords internal
.load_eps_boot <- function(label) {
  # Match files like: boot_eps_free_B200_s42.rds
  pattern <- sprintf("^boot_%s_B\\d+_s\\d+\\.rds$", label)
  files   <- list.files(boot_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) {
    warning(sprintf(".load_eps_boot: no bootstrap file for '%s' — SEs will be NA.",
                    label))
    return(NULL)
  }
  # Select by mtime (not name) to avoid lexicographic ordering bugs when B
  # or seed values differ in digit length (e.g., "B1000" < "B200" lexicographically).
  latest <- files[which.max(file.info(files)$mtime)]
  obj <- readRDS(latest)
  if (!is.data.frame(obj$summary))
    stop(sprintf(".load_eps_boot: '%s' missing $summary data frame", basename(latest)))
  obj
}

#' Extract named SE vector from bootstrap object (O(1) lookup).
.se_map <- function(boot_obj) {
  if (is.null(boot_obj) || is.null(boot_obj$summary)) return(NULL)
  setNames(boot_obj$summary$se, boot_obj$summary$quantity)
}

#' Look up SE for a named quantity.
.get_se <- function(se_vec, quantity) {
  if (is.null(se_vec)) return(NA_real_)
  v <- se_vec[[quantity]]
  if (is.null(v)) NA_real_ else v
}

cat("\n--- Loading fits and bootstraps ---\n")

eps_fits <- lapply(eps_labels, function(lbl)
  .load_latest_eps_fit(.fit_prefix[[lbl]]))
names(eps_fits) <- eps_labels

eps_boots <- lapply(eps_labels, .load_eps_boot)
names(eps_boots) <- eps_labels
eps_se    <- lapply(eps_boots, .se_map)

# Read B from the first available bootstrap object so this script does not need
# to hard-code it. $B and $master_seed are stored in each bootstrap .rds file.
.filtered_boots <- Filter(Negate(is.null), eps_boots)
.first_boot     <- if (length(.filtered_boots) > 0L) .filtered_boots[[1L]] else NULL
B <- if (!is.null(.first_boot)) .first_boot$B else {
  warning("No bootstrap files loaded — B unknown; footnote will show B=NA")
  NA_integer_
}
rm(.filtered_boots, .first_boot)

eps_implied <- lapply(eps_labels, function(lbl) {
  fit <- eps_fits[[lbl]]
  if (is.null(fit)) return(NULL)
  implied_tenure_contamination(fit$params)
})
names(eps_implied) <- eps_labels

# ==============================================================================
# 3. Column headers
# ==============================================================================

eps_col_headers <- c("", "(1)", "(2)", "(3)", "(4)")

eps_sub_headers <- c(
  paste0(
    "\\multicolumn{1}{l}{} & ",
    "\\multicolumn{1}{c}{Free} & ",
    "\\multicolumn{1}{c}{Stationary} & ",
    "\\multicolumn{1}{c}{CTMC-linked} & ",
    "\\multicolumn{1}{c}{Stat.+linked} \\\\"
  ),
  paste0(
    "\\multicolumn{1}{l}{} & ",
    "\\multicolumn{1}{c}{(7 param.)} & ",
    "\\multicolumn{1}{c}{(6 param.)} & ",
    "\\multicolumn{1}{c}{(5 param.)} & ",
    "\\multicolumn{1}{c}{(4 param.)} \\\\"
  ),
  "\\midrule"
)

star_note <- sprintf(
  "Bootstrap SE in parentheses ($B=%d$). Significance: $^{**}$ $p<0.01$, $^{*}$ $p<0.05$, $^{.}$ $p<0.10$.",
  B
)

# ==============================================================================
# 4. Table 1: Parameter Estimates
# ==============================================================================

cat("\n--- Building Table 1: Parameter Estimates ---\n")

# Helper: one parameter row across all 4 columns.
# Returns "---" (no SE) for constrained parameters.
.eps_param_row <- function(param_name, row_label,
                           constrained_when_stationary = FALSE,
                           constrained_when_linked     = FALSE) {
  cells <- lapply(eps_labels, function(lbl) {
    is_stat   <- .is_stationary[[lbl]]
    is_linked <- .is_linked[[lbl]]
    if ((constrained_when_stationary && is_stat) ||
        (constrained_when_linked     && is_linked)) {
      return(c("---", ""))
    }
    fit <- eps_fits[[lbl]]
    v   <- if (!is.null(fit)) fit$params[[param_name]] else NA_real_
    if (is.null(v)) v <- NA_real_
    se  <- .get_se(eps_se[[lbl]], param_name)
    .fmt_param(v, se)
  })
  list(
    est = c(row_label, vapply(cells, `[[`, character(1L), 1L)),
    se  = c("",        vapply(cells, `[[`, character(1L), 2L))
  )
}

param_rows <- list(
  list(header = "Transition parameters", rows = list(
    .eps_param_row("theta0", "$\\theta_0$ (entry)"),
    .eps_param_row("theta1", "$\\theta_1$ (persistence)"),
    # alpha is derived (not free) in stationary variants
    .eps_param_row("alpha",  "$\\alpha$ (initial empl.)",
                   constrained_when_stationary = TRUE)
  )),
  list(header = "Misclassification \\& contamination", rows = list(
    .eps_param_row("pi",  "$\\pi$ (misclassification)"),
    .eps_param_row("eps", "$\\varepsilon$ (contamination)")
  )),
  list(header = "Exponential spell rates", rows = list(
    # lambda_g and lambda_d are derived in linked variants
    .eps_param_row("lambda_g", "$\\lambda_g$ (employment rate)",
                   constrained_when_linked = TRUE),
    .eps_param_row("lambda_d", "$\\lambda_d$ (non-employment rate)",
                   constrained_when_linked = TRUE)
  )),
  list(header = NULL, rows = list(
    list(
      est = c("Log-likelihood",
              vapply(eps_labels, function(lbl) {
                fit <- eps_fits[[lbl]]
                if (is.null(fit)) "---" else .fmt_ll(fit$loglik)
              }, character(1L))),
      se  = rep("", 5L)
    ),
    list(
      est = c("$N$", rep(if (is.na(N)) "---" else .fmt_n(N), 4L)),
      se  = rep("", 5L)
    )
  ))
)

.write_latex_table(
  col_headers = eps_col_headers,
  row_data    = param_rows,
  caption     = paste0(
    "EM tenure contamination model: parameter estimates. ",
    "Preferred specification: column~(4) (stationary~+~CTMC-linked)."
  ),
  label       = "tab:eps_params",
  path        = file.path(tables_dir, "table_eps_params.tex"),
  sub_headers = eps_sub_headers,
  note        = star_note
)

# ==============================================================================
# 5. Table 2: Implied Probabilities
# ==============================================================================

cat("\n--- Building Table 2: Implied Probabilities ---\n")

# Helper: one implied quantity row across all 4 columns (as percentages)
.eps_implied_row_pct <- function(qty_name, row_label) {
  cells <- lapply(eps_labels, function(lbl) {
    imp <- eps_implied[[lbl]]
    v   <- if (!is.null(imp)) imp[[qty_name]] else NA_real_
    if (is.null(v)) v <- NA_real_
    se  <- .get_se(eps_se[[lbl]], qty_name)
    .fmt_pct(v, se)
  })
  list(
    est = c(row_label, vapply(cells, `[[`, character(1L), 1L)),
    se  = c("",        vapply(cells, `[[`, character(1L), 2L))
  )
}

implied_rows <- list(
  list(header = "Labour-market flow rates", rows = list(
    .eps_implied_row_pct("entry_rate",      "Entry rate ($\\theta_0$, \\%)"),
    .eps_implied_row_pct("exit_rate",       "Exit rate ($1-\\theta_1$, \\%)"),
    .eps_implied_row_pct("employment_rate", "Employment rate ($\\alpha^*$, \\%)")
  )),
  list(header = "Measurement error rates", rows = list(
    .eps_implied_row_pct("pi",  "$\\pi$: misclassification (\\%)"),
    .eps_implied_row_pct("eps", "$\\varepsilon$: contamination (\\%)")
  )),
  list(header = NULL, rows = list(
    list(
      est = c("$N$", rep(if (is.na(N)) "---" else .fmt_n(N), 4L)),
      se  = rep("", 5L)
    )
  ))
)

.write_latex_table(
  col_headers = eps_col_headers,
  row_data    = implied_rows,
  caption     = paste0(
    "EM tenure contamination model: implied probabilities (\\%). ",
    "Employment rate is steady-state $\\alpha^* = \\theta_0 / (\\theta_0 + 1 - \\theta_1)$."
  ),
  label       = "tab:eps_implied",
  path        = file.path(tables_dir, "table_eps_implied.tex"),
  sub_headers = eps_sub_headers,
  note        = paste0(
    "All quantities in \\%. ",
    star_note
  )
)

# ==============================================================================
# 6. Table 3: Tenure Interpretation
# ==============================================================================

cat("\n--- Building Table 3: Tenure Interpretation ---\n")

# Helper: duration row in years (2 d.p.)
.eps_dur_row_yrs <- function(qty_name, row_label) {
  cells <- lapply(eps_labels, function(lbl) {
    imp <- eps_implied[[lbl]]
    v   <- if (!is.null(imp)) imp[[qty_name]] else NA_real_
    if (is.null(v)) v <- NA_real_
    se  <- .get_se(eps_se[[lbl]], qty_name)
    .fmt_dur(v, se, digits = 2L)
  })
  list(
    est = c(row_label, vapply(cells, `[[`, character(1L), 1L)),
    se  = c("",        vapply(cells, `[[`, character(1L), 2L))
  )
}

# Helper: duration row in months (1 d.p.)
.eps_dur_row_mths <- function(qty_name, row_label) {
  cells <- lapply(eps_labels, function(lbl) {
    imp <- eps_implied[[lbl]]
    v   <- if (!is.null(imp)) imp[[qty_name]] else NA_real_
    if (is.null(v)) v <- NA_real_
    se  <- .get_se(eps_se[[lbl]], qty_name)
    .fmt_months(v, se)
  })
  list(
    est = c(row_label, vapply(cells, `[[`, character(1L), 1L)),
    se  = c("",        vapply(cells, `[[`, character(1L), 2L))
  )
}

# Helper: probability row as percentage (2 d.p.)
.eps_prob_row <- function(qty_name, row_label) {
  .eps_implied_row_pct(qty_name, row_label)
}

# Helper: contaminated per 1000 (1 d.p., no percentage scaling)
.eps_contam_row <- function(qty_name, row_label) {
  cells <- lapply(eps_labels, function(lbl) {
    imp <- eps_implied[[lbl]]
    v   <- if (!is.null(imp)) imp[[qty_name]] else NA_real_
    if (is.null(v)) v <- NA_real_
    se  <- .get_se(eps_se[[lbl]], qty_name)
    # Contaminated per 1000 is already on the ×1000 scale; show 1 d.p.
    .fmt_dur(v, se, digits = 1L)
  })
  list(
    est = c(row_label, vapply(cells, `[[`, character(1L), 1L)),
    se  = c("",        vapply(cells, `[[`, character(1L), 2L))
  )
}

tenure_rows <- list(
  list(header = "Employment spell duration ($\\hat{\\lambda}_g$)", rows = list(
    .eps_dur_row_yrs( "mean_spell_g_years",    "Mean spell (years)"),
    .eps_dur_row_mths("mean_spell_g_months",   "Mean spell (months)"),
    .eps_dur_row_mths("median_spell_g_months", "Median spell (months)")
  )),
  list(header = "Non-employment spell duration ($\\hat{\\lambda}_d$)", rows = list(
    .eps_dur_row_yrs( "mean_spell_d_years",  "Mean spell (years)"),
    .eps_dur_row_mths("mean_spell_d_months", "Mean spell (months)")
  )),
  list(header = "Tenure measurement quality ($\\hat{\\varepsilon}$)", rows = list(
    .eps_prob_row("p_clock_consistent",
                  "P(clock-consistent), $1-\\varepsilon$ (\\%)"),
    .eps_prob_row("p_pair_clean",
                  "P(spell-pair both clean), $(1-\\varepsilon)^2$ (\\%)"),
    .eps_prob_row("p_triple_clean",
                  "P(K=3 all clean), $(1-\\varepsilon)^3$ (\\%)"),
    .eps_contam_row("contaminated_per_1000",
                    "Contaminated per 1,000 observations")
  )),
  list(header = NULL, rows = list(
    list(
      est = c("$N$", rep(if (is.na(N)) "---" else .fmt_n(N), 4L)),
      se  = rep("", 5L)
    )
  ))
)

.write_latex_table(
  col_headers = eps_col_headers,
  row_data    = tenure_rows,
  caption     = paste0(
    "EM tenure contamination model: tenure and spell-length interpretation. ",
    "Durations derived from $\\hat{\\lambda}_g$ and $\\hat{\\lambda}_d$ assuming ",
    "Exponential distributions. Contamination probabilities derived from ",
    "$\\hat{\\varepsilon}$."
  ),
  label       = "tab:eps_tenure",
  path        = file.path(tables_dir, "table_eps_tenure.tex"),
  sub_headers = eps_sub_headers,
  note        = paste0(
    "Spell durations are $1/\\hat{\\lambda}$ (mean) and $\\log(2)/\\hat{\\lambda}$ ",
    "(median) under the Exp$(\\hat{\\lambda})$ model. ",
    star_note
  )
)

cat("\nDone. Tables written to:", tables_dir, "\n")
