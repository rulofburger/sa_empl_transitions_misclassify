# ==============================================================================
# Table builder: parameter estimates and implied probabilities
# Created: 2026-05-06 | Revised: 2026-05-07
#
# Reads point estimates and analytical/bootstrap inference results, then
# produces publication-ready LaTeX tables with:
#   - Significance stars (. 10%, * 5%, ** 1%) testing H0: θ=0
#   - Implied probabilities expressed as percentages (2 d.p.)
#   - AMEs expressed as percentage points (2 d.p.)
#   - Log-likelihoods in millions (e.g., -357.2M)
#   - Numbered columns: (1), (2), ...
#   - No-misclassification models first
#   - N (observations) row in all tables
#
# Tables produced:
#   Paper Table 2: Baseline implied probabilities         (table_baseline_implied.tex)
#   Paper Table 3: Baseline parameter estimates           (table_baseline_params.tex)
#
# Usage (from project root):
#   Rscript build_tables_v2.R
#
# Baseline prerequisite: estimate_baseline_pipeline.R must have been run first.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------

library(here)
library(dplyr)
library(ggplot2)  # required by ingest_data_3waves_SA.R diagnostic plots

# Null-coalescing (also defined in EM-baseline/R/utils.R via source_all.R;
# explicit here to avoid silent coupling to the source chain).
`%||%` <- function(a, b) if (!is.null(a)) a else b

source(here::here("EM-baseline",     "R", "source_all.R"))
source(here::here("EM-baseline-ext", "R", "source_all.R"))

boot_baseline_dir <- here::here("EM-baseline",     "output", "results", "bootstrap")
boot_ext_dir      <- here::here("EM-baseline-ext", "output", "results", "bootstrap")

.detect_B <- function(boot_dir) {
  fls <- list.files(boot_dir, pattern = "boot_.*_B[0-9]+\\.rds$")
  if (length(fls) == 0L) return(NULL)

  as.integer(sub(".*_B([0-9]+)\\.rds$", "\\1", fls[1]))
}
B_baseline <- .detect_B(boot_baseline_dir)
B_ext      <- .detect_B(boot_ext_dir)
B          <- B_baseline %||% B_ext %||% 200L
if (!is.null(B_baseline) && !is.null(B_ext) && B_baseline != B_ext) {
  warning(sprintf("Detected different B in baseline (%d) vs ext (%d) bootstrap dirs.",
                  B_baseline, B_ext))
}
cat(sprintf("B detected from bootstrap files: %d\n", B))

tables_baseline_dir <- here::here("EM-baseline",     "output", "tables")
tables_ext_dir      <- here::here("EM-baseline-ext", "output", "tables")
dir.create(tables_baseline_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_ext_dir,      recursive = TRUE, showWarnings = FALSE)

# Load data and compute sample sizes
source(here::here("scripts", "ingest_data_3waves_SA.R"))

sector_source_path <- here::here("data", "raw", "QLFSmerged_mapped.rds")
if (!file.exists(sector_source_path))
  stop("Missing sector source: data/raw/QLFSmerged_mapped.rds")
sector_source <- readRDS(sector_source_path)
df_qlfs <- attach_transition_informal_sector(df_qlfs, sector_source)
rm(sector_source)

# Baseline sample: y1, y2, y3, weight valid
keep_baseline <- !is.na(df_qlfs$y1) & !is.na(df_qlfs$y2) & !is.na(df_qlfs$y3) &
                 !is.na(df_qlfs$weight) & df_qlfs$weight > 0
N_baseline <- sum(keep_baseline)
df_baseline_check <- data.frame(
  y1 = as.integer(df_qlfs$y1[keep_baseline]),
  y2 = as.integer(df_qlfs$y2[keep_baseline]),
  y3 = as.integer(df_qlfs$y3[keep_baseline]),
  weight = as.numeric(df_qlfs$weight[keep_baseline])
)
baseline_signature <- baseline_sample_signature(collapse_baseline_cells(df_baseline_check))

# Extension sample: additionally requires covariates
keep <- keep_baseline &
        !is.na(df_qlfs$age1)   & !is.na(df_qlfs$age2)   & !is.na(df_qlfs$age3) &
        !is.na(df_qlfs$educ1)  & !is.na(df_qlfs$educ2)  & !is.na(df_qlfs$educ3)
df_ext        <- df_qlfs[keep, , drop = FALSE]
df_ext$y1     <- as.integer(df_ext$y1)
df_ext$y2     <- as.integer(df_ext$y2)
df_ext$y3     <- as.integer(df_ext$y3)
df_ext$weight <- as.numeric(df_ext$weight)
for (nm in c("contracttype1", "contracttype2"))
  df_ext[[nm]] <- ifelse(is.na(df_ext[[nm]]), 0L, as.integer(df_ext[[nm]]))
df_ext <- as.data.frame(df_ext)
N_ext  <- nrow(df_ext)

cv_set1 <- prepare_covariate_matrix(df_ext, covariate_set = 1L)
cv_set2 <- prepare_covariate_matrix(df_ext, covariate_set = 2L)
cv_set3 <- prepare_covariate_matrix(df_ext, covariate_set = 3L)
X_list  <- list(X1 = cv_set1$X, X2 = cv_set2$X,
                X3 = cv_set3$X_transition)
rm(df_qlfs)

cat(sprintf("N_baseline = %s, N_ext = %s\n",
            formatC(N_baseline, big.mark = ","), formatC(N_ext, big.mark = ",")))

# ------------------------------------------------------------------------------
# 1. Table formatting helpers
# ------------------------------------------------------------------------------

# String concatenation helper
`%+%` <- function(a, b) paste0(a, b)

# Significance threshold constants (standard normal quantiles)
.CRIT_p01 <- qnorm(0.995)  # ≈ 2.576  (p < 0.01, two-sided)
.CRIT_p05 <- qnorm(0.975)  # ≈ 1.960  (p < 0.05)
.CRIT_p10 <- qnorm(0.945)  # ≈ 1.645  (p < 0.10)

#' Significance stars based on z = |est/se|
.stars <- function(est, se) {
  if (is.na(est) || is.na(se) || se <= 0) return("")
  z <- abs(est / se)
  if (z > .CRIT_p01) return("$^{**}$")
  if (z > .CRIT_p05) return("$^{*}$")
  if (z > .CRIT_p10) return("$^{.}$")
  return("")
}

#' Common formatter: scale by factor, format to digits d.p., append stars
.fmt_estimate <- function(est, se, factor = 1, digits = 4L) {
  if (is.na(est)) return(c("---", ""))
  star    <- .stars(est, se)
  est_str <- paste0(formatC(est * factor, format = "f", digits = digits), star)
  se_str  <- if (!is.na(se)) {
    sprintf("(%s)", formatC(se * factor, format = "f", digits = digits))
  } else {
    ""
  }
  c(est_str, se_str)
}

#' Format a raw parameter estimate (4 d.p.) with stars + SE row
.fmt_param <- function(est, se, digits = 4L) .fmt_estimate(est, se, factor = 1,   digits = digits)
.fmt_param_plain <- function(est, se, digits = 4L) {
  if (is.na(est)) return(c("---", ""))
  c(formatC(est, format = "f", digits = digits),
    if (is.na(se)) "" else sprintf("(%s)", formatC(se, format = "f", digits = digits)))
}

#' Format a probability as percentage (×100, 2 d.p.) with stars + SE row
.fmt_pct   <- function(est, se, digits = 2L) .fmt_estimate(est, se, factor = 100, digits = digits)
.fmt_pct_plain <- function(est, se, digits = 2L) {
  if (is.na(est)) return(c("---", ""))
  c(formatC(est * 100, format = "f", digits = digits),
    if (is.na(se)) "" else sprintf("(%s)", formatC(se * 100, format = "f", digits = digits)))
}

#' Format a marginal effect as percentage points (×100, 2 d.p.) with stars
.fmt_pp    <- function(est, se, digits = 2L) .fmt_pct(est, se, digits)

#' Format log-likelihood in millions
.fmt_ll <- function(ll, digits = 1L) {
  if (is.na(ll) || is.null(ll)) return("---")
  sprintf(paste0("%.", digits, "fM"), ll / 1e6)
}

#' Format sample size with comma separator
.fmt_n <- function(n) {
  formatC(n, format = "d", big.mark = ",")
}

#' Write a LaTeX booktabs table to file
.write_latex_table <- function(col_headers, row_data, caption, label, path,
                               col_align = NULL, note = NULL,
                               sub_headers = NULL) {
  n_cols <- length(col_headers)
  if (is.null(col_align)) {
    col_align <- paste0("l", paste(rep("c", n_cols - 1), collapse = ""))
  }

  # Use a pre-allocated list buffer to avoid O(n²) vector copies from c(lines,...)
  buf <- vector("list", 200L)
  i   <- 0L
  push <- function(x) { i <<- i + 1L; buf[[i]] <<- x }

  push("\\begin{table}[htbp]")
  push("\\centering")
  push(sprintf("\\caption{%s}", caption))
  push(sprintf("\\label{%s}", label))
  push(sprintf("\\begin{tabular}{%s}", col_align))
  push("\\toprule")
  push(paste(col_headers, collapse = " & "))
  push("\\\\")

  # Optional sub-header rows (e.g., model descriptions under column numbers)
  if (!is.null(sub_headers)) for (sh in sub_headers) push(sh)

  push("\\midrule")

  for (block in row_data) {
    if (!is.null(block$header))
      push(sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\", n_cols, block$header))
    for (row in block$rows) {
      est_row <- row[["est"]]
      se_row  <- row[["se"]]
      push(paste(est_row, collapse = " & "))
      push("\\\\")
      # Only add SE row if it has content
      if (length(se_row) > 0 && !all(se_row == "")) {
        push(paste(c(se_row[1], se_row[-1]), collapse = " & "))
        push("\\\\")
      }
    }
  }

  push("\\bottomrule")
  push("\\end{tabular}")
  if (!is.null(note)) {
    push("\\begin{minipage}{0.98\\linewidth}")
    push(sprintf("\\footnotesize \\textit{Note:} %s", note))
    push("\\end{minipage}")
  }
  push("\\end{table}")

  cat(paste(unlist(buf[seq_len(i)]), collapse = "\n"), "\n", file = path)
  cat(sprintf("Written: %s\n", path))
}

#' Load a point-estimate fit .rds, with structural validation.
.load_fit <- function(label, subdir) {
  p <- here::here(subdir, "output", "results", sprintf("fit_%s.rds", label))
  if (!file.exists(p)) return(NULL)
  obj <- readRDS(p)
  if (!is.list(obj) || is.null(obj$params) || is.null(obj$loglik))
    stop(sprintf(".load_fit: fit_%s.rds is missing $params or $loglik.", label))
  obj
}

#' Load bootstrap summary .rds, with schema validation.
.load_boot <- function(label, boot_dir) {
  path <- file.path(boot_dir, sprintf("boot_%s_B%d.rds", label, B))
  if (!file.exists(path)) return(NULL)
  obj <- readRDS(path)
  if (!is.list(obj) || !is.data.frame(obj$summary))
    stop(sprintf(".load_boot: '%s' is missing required $summary data frame.",
                 basename(path)))
  obj
}

.load_analytic <- function(label, subdir = "EM-baseline-ext") {
  path <- here::here(subdir, "output", "results",
                     sprintf("analytical_se_%s.rds", label))
  if (!file.exists(path)) return(NULL)
  obj <- readRDS(path)
  if (!is.list(obj) || !is.data.frame(obj$summary))
    stop(sprintf(".load_analytic: '%s' has an invalid schema.", basename(path)))
  obj
}

#' Look up SE for a quantity from a bootstrap summary (data.frame or named vector).
.get_se <- function(summary_or_map, quantity) {
  if (is.null(summary_or_map)) return(NA_real_)
  if (is.data.frame(summary_or_map)) {
    idx <- which(summary_or_map$quantity == quantity)
    if (length(idx) == 0L) NA_real_ else summary_or_map$se[idx[1L]]
  } else {
    # Named numeric vector (pre-indexed se_map)
    if (is.null(names(summary_or_map)) || !quantity %in% names(summary_or_map))
      NA_real_
    else
      unname(summary_or_map[[quantity]])
  }
}

#' Pre-index bootstrap summary as a named SE vector for O(1) lookup.
.se_map <- function(boot_obj) {
  if (is.null(boot_obj) || is.null(boot_obj$summary)) return(NULL)
  setNames(boot_obj$summary$se, boot_obj$summary$quantity)
}

# AME covariate label map (raw name → publication label)
.cov_label_map <- c(
  intercept      = "Intercept",
  age            = "Age",
  age_sq         = "Age$^2$",
  educ           = "Education",

  race_2         = "Race: Coloured",
  race_3         = "Race: Indian",
  race_4         = "Race: White",
  female         = "Female",
  contracttype_1 = "Permanent contract",
  informal_sector = "Informal sector"
)

# ------------------------------------------------------------------------------
# 2. Table 1 & 2: Baseline models
# ------------------------------------------------------------------------------

cat("\n--- Building baseline tables ---\n")

# Column order: no-misclassification first, then symmetric, then asymmetric
baseline_labels <- c("none_stat", "none_free", "sym_stat", "sym_free",
                     "asym_stat", "asym_free")

baseline_col_headers <- c("",
                           "(1)", "(2)", "(3)", "(4)", "(5)", "(6)")
baseline_sub_headers <- c(
  paste(c("", "\\multicolumn{2}{c}{No miscl.}",
          "\\multicolumn{2}{c}{Symmetric}",
          "\\multicolumn{2}{c}{Asymmetric}"), collapse = " & ") %+% " \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}",
  paste(c("", "Stat.", "Free", "Stat.", "Free", "Stat.", "Free"),
        collapse = " & ") %+% " \\\\"
)

# Load fits and bootstraps
bl_fits  <- lapply(baseline_labels, function(lbl) .load_fit(lbl, "EM-baseline"))
names(bl_fits) <- baseline_labels

for (lbl in baseline_labels) {
  fit <- bl_fits[[lbl]]
  if (is.null(fit)) stop("Missing baseline fit: ", lbl)
  if (!identical(fit$estimator, "direct_eight_cell_mle"))
    stop("Baseline fit is not a validated direct MLE: ", lbl)
  if (!isTRUE(fit$converged)) stop("Baseline fit failed diagnostics: ", lbl)
  if (!identical(fit$sample$signature, baseline_signature))
    stop("Baseline fit was estimated on a different or stale sample: ", lbl)
}
check_baseline_nesting(bl_fits)

bl_analytic <- lapply(baseline_labels, function(lbl)
  .load_analytic(lbl, "EM-baseline"))
names(bl_analytic) <- baseline_labels
if (any(vapply(bl_analytic, is.null, logical(1L))))
  stop("Analytical baseline inference is missing; rerun estimate_baseline_pipeline.R")
bl_se <- lapply(bl_analytic, .se_map)

model_types_bl <- c(none_stat = "none", none_free = "none",
                    sym_stat = "symmetric", sym_free = "symmetric",
                    asym_stat = "asymmetric", asym_free = "asymmetric")

bl_implied <- lapply(baseline_labels, function(lbl) {
  fit <- bl_fits[[lbl]]
  if (is.null(fit)) return(NULL)
  implied_baseline(fit$params, model_types_bl[[lbl]])
})
names(bl_implied) <- baseline_labels

# Helper: parameter row with stars (raw scale, 4 d.p.)
.bl_param_row <- function(param_name, row_label) {
  cells <- lapply(baseline_labels, function(lbl) {
    fit <- bl_fits[[lbl]]
    v   <- if (!is.null(fit)) fit$params[[param_name]] else NA_real_
    if (is.null(v)) v <- NA_real_
    se  <- .get_se(bl_se[[lbl]], param_name)
    .fmt_param_plain(v, se)
  })
  list(
    est = c(row_label, sapply(cells, `[[`, 1L)),
    se  = c("",        sapply(cells, `[[`, 2L))
  )
}

# Paper Table 3: baseline parameters
param_rows_bl <- list(
  list(header = "Transition parameters", rows = list(
    .bl_param_row("theta0", "$\\theta_0$ (entry)"),
    .bl_param_row("theta1", "$\\theta_1$ (persistence)"),
    .bl_param_row("alpha",  "$\\alpha$ (initial empl.)")
  )),
  list(header = "Misclassification", rows = list(
    .bl_param_row("pi",  "$\\pi$ (symmetric)"),
    .bl_param_row("pi0", "$\\pi_0$ (false positive)"),
    .bl_param_row("pi1", "$\\pi_1$ (false negative)")
  )),
  list(header = NULL, rows = list(
    list(
      est = c("Log-likelihood", vapply(baseline_labels, function(lbl) {
        fit <- bl_fits[[lbl]]
        if (is.null(fit)) "---" else .fmt_ll(fit$loglik, digits = 3L)
      }, character(1L))),
      se = rep("", 7L)
    ),
    list(
      est = c("$N$", rep(.fmt_n(N_baseline), 6L)),
      se  = rep("", 7L)
    )
  ))
)

.write_latex_table(
  col_headers = baseline_col_headers,
  row_data    = param_rows_bl,
  caption     = "Baseline latent-state model: parameter estimates",
  label       = "tab:baseline_params",
  path        = file.path(tables_baseline_dir, "table_baseline_params.tex"),
  sub_headers = baseline_sub_headers,
  note        = paste0("Analytical SE in parentheses use an individual-level ",
                       "survey-weighted sandwich covariance matrix. $\\alpha$ is the ",
                       "initial employment probability; in stationary columns it is ",
                       "derived from the transition parameters. The calculation does ",
                       "not incorporate strata or PSU clustering.")
)

# Paper Table 2: baseline implied probabilities (as %)
.bl_implied_row <- function(qty_name, row_label) {
  cells <- lapply(baseline_labels, function(lbl) {
    imp <- bl_implied[[lbl]]
    v   <- if (!is.null(imp)) imp[[qty_name]] else NA_real_
    if (is.null(v)) v <- NA_real_
    se  <- .get_se(bl_se[[lbl]], qty_name)
    .fmt_pct_plain(v, se)
  })
  list(
    est = c(row_label, sapply(cells, `[[`, 1L)),
    se  = c("",        sapply(cells, `[[`, 2L))
  )
}

implied_rows_bl <- list(
  list(header = NULL, rows = list(
    .bl_implied_row("entry_rate",      "Entry rate ($\\theta_0$)"),
    .bl_implied_row("exit_rate",       "Exit rate ($1-\\theta_1$)"),
    .bl_implied_row("employment_rate", "Employment rate ($\\alpha^*$)"),
    .bl_implied_row("pi",              "$\\pi$ (misclassification)"),
    .bl_implied_row("pi0",             "$\\pi_0$ (false positive)"),
    .bl_implied_row("pi1",             "$\\pi_1$ (false negative)")
  )),
  list(header = NULL, rows = list(
    list(
      est = c("$N$", rep(.fmt_n(N_baseline), 6L)),
      se  = rep("", 7L)
    )
  ))
)

.write_latex_table(
  col_headers = baseline_col_headers,
  row_data    = implied_rows_bl,
  caption     = "Baseline latent-state model: implied probabilities (\\%)",
  label       = "tab:baseline_implied",
  path        = file.path(tables_baseline_dir, "table_baseline_implied.tex"),
  sub_headers = baseline_sub_headers,
  note        = paste0("Estimates in \\%. Analytical SE in parentheses use an ",
                       "individual-level survey-weighted sandwich covariance matrix ",
                       "and the delta method. Employment rate is steady-state ",
                       "$\\alpha^* = \\theta_0 / (\\theta_0 + 1 - \\theta_1)$, ",
                       "including in free-initial-probability columns. The calculation ",
                       "does not incorporate strata or PSU clustering.")
)

# ------------------------------------------------------------------------------
# 3. Covariate models: implied probabilities (split into stat/free panels)
# ------------------------------------------------------------------------------

cat("\n--- Building covariate tables ---\n")

# Set 3 has transition-varying contract and sector measures, so stationarity is
# not imposed on it. The main table consistently compares free-alpha models.
cov_labels_stat <- c("cov_s1_non_stat", "cov_s2_non_stat",
                     "cov_s1_sym_stat", "cov_s2_sym_stat")
cov_labels_free <- c("cov_s1_non_free", "cov_s2_non_free", "cov_s3_non_free",
                     "cov_s1_sym_free", "cov_s2_sym_free", "cov_s3_sym_free")
cov_labels_all  <- unique(c(cov_labels_stat, cov_labels_free))

cov_col_headers_panel <- c("", "(1)", "(2)", "(3)", "(4)", "(5)", "(6)")
cov_sub_headers_panel <- c(
  paste(c("", "\\multicolumn{3}{c}{No miscl.}",
          "\\multicolumn{3}{c}{Symmetric}"), collapse = " & ") %+% " \\\\",
  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}",
  paste(c("", "Set 1", "Set 2", "Set 3", "Set 1", "Set 2", "Set 3"),
        collapse = " & ") %+% " \\\\"
)

# Load all covariate fits and bootstraps
cov_fits <- lapply(cov_labels_all, function(lbl) .load_fit(lbl, "EM-baseline-ext"))
names(cov_fits) <- cov_labels_all

cov_boots <- lapply(cov_labels_all, function(lbl) .load_boot(lbl, boot_ext_dir))
names(cov_boots) <- cov_labels_all
cov_analytic <- lapply(cov_labels_all, .load_analytic)
names(cov_analytic) <- cov_labels_all
cov_se <- lapply(cov_labels_all, function(lbl) {
  # Prefer bootstrap inference when it exists; otherwise use the analytical
  # individual-level sandwich/delta approximation.
  .se_map(cov_boots[[lbl]]) %||% .se_map(cov_analytic[[lbl]])
})
names(cov_se) <- cov_labels_all

cov_mt <- c(
  cov_s1_sym_stat = "symmetric", cov_s1_sym_free = "symmetric",
  cov_s1_non_stat = "none",      cov_s1_non_free = "none",
  cov_s2_sym_stat = "symmetric", cov_s2_sym_free = "symmetric",
  cov_s2_non_stat = "none",      cov_s2_non_free = "none",
  cov_s3_sym_stat = "symmetric", cov_s3_sym_free = "symmetric",
  cov_s3_non_stat = "none",      cov_s3_non_free = "none"
)
cov_xmat <- c(
  cov_s1_sym_stat = "X1", cov_s1_sym_free = "X1",
  cov_s1_non_stat = "X1", cov_s1_non_free = "X1",
  cov_s2_sym_stat = "X2", cov_s2_sym_free = "X2",
  cov_s2_non_stat = "X2", cov_s2_non_free = "X2",
  cov_s3_sym_stat = "X3", cov_s3_sym_free = "X3",
  cov_s3_non_stat = "X3", cov_s3_non_free = "X3"
)

cov_implied <- lapply(cov_labels_all, function(lbl) {
  fit <- cov_fits[[lbl]]
  if (is.null(fit)) return(NULL)
  X <- X_list[[cov_xmat[[lbl]]]]
  tryCatch(implied_covariates(fit$params, X, cov_mt[[lbl]],
                              df = df_ext, gamma = fit$gamma),
           error = function(e) NULL)
})
names(cov_implied) <- cov_labels_all

# Helper: one implied probability row for a panel (as %)
.cov_implied_row_panel <- function(qty_name, row_label, panel_labels) {
  cells <- lapply(panel_labels, function(lbl) {
    imp <- cov_implied[[lbl]]
    v   <- if (!is.null(imp)) imp[[qty_name]] else NA_real_
    if (is.null(v)) v <- NA_real_
    se  <- .get_se(cov_se[[lbl]], qty_name)
    .fmt_pct(v, se)
  })
  list(
    est = c(row_label, sapply(cells, `[[`, 1L)),
    se  = c("",        sapply(cells, `[[`, 2L))
  )
}

# Covariate inclusion indicator rows
.cov_indicator_rows <- function(panel_labels) {
  # Determine which set each column belongs to
  set_of <- function(lbl) {
    if (grepl("_s1_", lbl)) return(1L)
    if (grepl("_s2_", lbl)) return(2L)
    return(3L)
  }
  sets <- vapply(panel_labels, set_of, integer(1L))

  # Indicator: Yes if set includes covariate group
  # Set 1: Age, Age^2, Educ
  # Set 2: Set 1 + Race, Gender
  # Set 3: Set 2 + contract type and informal sector in the exit equation
  make_row <- function(label, min_set) {
    cells <- ifelse(sets >= min_set, "Yes", "No")
    list(est = c(label, cells), se = rep("", length(panel_labels) + 1L))
  }
  list(
    make_row("Age, Age$^2$, Educ.", 1L),
    make_row("+ Race, Gender", 2L),
    make_row("+ Contract type, informal sector (exit only)", 3L)
  )
}

# Build and write Panel A (Stationary)
.build_cov_implied_panel <- function(panel_labels, panel_name, suffix, caption_extra) {
  row_data <- list(
    list(header = NULL, rows = list(
      .cov_implied_row_panel("mean_entry_rate",      "Mean entry rate (\\%)", panel_labels),
      .cov_implied_row_panel("mean_exit_rate",       "Mean exit rate (\\%)", panel_labels),
      .cov_implied_row_panel("mean_employment_rate", "Mean employment rate (\\%)", panel_labels),
      .cov_implied_row_panel("pi",                   "$\\pi$ misclassification (\\%)", panel_labels)
    )),
    list(header = NULL, rows = list(
      list(
        est = c("Log-likelihood", vapply(panel_labels, function(lbl) {
          fit <- cov_fits[[lbl]]
          if (is.null(fit)) "---" else .fmt_ll(fit$loglik)
        }, character(1L))),
        se = rep("", length(panel_labels) + 1L)
      ),
      list(
        est = c("$N$", rep(.fmt_n(N_ext), length(panel_labels))),
        se  = rep("", length(panel_labels) + 1L)
      )
    )),
    list(header = "Covariates included", rows = .cov_indicator_rows(panel_labels))
  )

  .write_latex_table(
    col_headers = cov_col_headers_panel,
    row_data    = row_data,
    caption     = sprintf("Covariate EM models: implied probabilities (%s)", caption_extra),
    label       = sprintf("tab:cov_implied_%s", suffix),
    path        = file.path(tables_ext_dir, sprintf("table_cov_implied_%s.tex", suffix)),
    sub_headers = cov_sub_headers_panel,
    note        = sprintf("Estimates in \\%%. Bootstrap SE in parentheses ($B=%d$). Significance: $^{**}$ $p<0.01$, $^{*}$ $p<0.05$, $^{.}$ $p<0.10$.", B)
  )
}

.build_cov_implied_panel(cov_labels_free, "Free alpha", "free",
                         "free $\\alpha$")

# ---- Main Table 4: risk-set and survey-weighted transition summaries --------

.cov_plain_row <- function(qty_name, row_label, panel_labels, formatter = .fmt_pct) {
  cells <- lapply(panel_labels, function(lbl) {
    imp <- cov_implied[[lbl]]
    v <- if (!is.null(imp)) imp[[qty_name]] else NA_real_
    if (is.null(v)) v <- NA_real_
    formatter(v, .get_se(cov_se[[lbl]], qty_name))
  })
  list(est = c(row_label, vapply(cells, `[[`, character(1L), 1L)),
       se  = c("", vapply(cells, `[[`, character(1L), 2L)))
}

main_table4_rows <- list(
  list(header = "Risk-set weighted transition probabilities", rows = list(
    .cov_plain_row("mean_entry_rate", "Entry rate", cov_labels_free, .fmt_pct_plain),
    .cov_plain_row("mean_exit_rate", "Exit rate", cov_labels_free, .fmt_pct_plain),
    .cov_plain_row("mean_employment_rate", "Posterior employed risk share", cov_labels_free, .fmt_pct_plain)
  )),
  list(header = "Employment turnover", rows = list(
    .cov_plain_row("total_churn_flow", "Total employment churn", cov_labels_free, .fmt_pct_plain)
  )),
  list(header = "Other model quantities", rows = list(
    .cov_plain_row("pi", "Misclassification probability $\\pi$", cov_labels_free, .fmt_pct_plain),
    .cov_plain_row("alpha", "Initial employment probability $\\alpha$", cov_labels_free, .fmt_pct_plain),
    list(est = c("Log-likelihood", vapply(cov_labels_free, function(lbl) {
      fit <- cov_fits[[lbl]]
      if (is.null(fit)) "---" else .fmt_ll(fit$loglik)
    }, character(1L))), se = rep("", 7L)),
    list(est = c("$N$", rep(.fmt_n(N_ext), 6L)), se = rep("", 7L))
  )),
  list(header = "Covariates included", rows = .cov_indicator_rows(cov_labels_free))
)

.write_latex_table(
  col_headers = cov_col_headers_panel,
  row_data = main_table4_rows,
  caption = "Risk-set and survey-weighted quarterly employment transitions",
  label = "tab:cov_risk_weighted",
  path = file.path(tables_ext_dir, "table_cov_risk_weighted.tex"),
  sub_headers = cov_sub_headers_panel,
  note = paste0("All specifications estimate the initial employment probability freely. ",
                "Rates average predicted hazards using survey weights and posterior ",
                "probabilities of belonging to the relevant origin-state risk set. ",
                "Total churn is the expected number of entry and exit transitions per ",
                "person-quarter. ",
                "SE are bootstrap estimates when available and otherwise use ",
                "an individual-level survey-weighted sandwich/delta approximation. ",
                "The approximation does not incorporate strata or PSU clustering.")
)

# ---- Appendix: raw probit coefficients and distributional summaries --------

.raw_parameter_row <- function(block, index, row_label) {
  cells <- lapply(cov_labels_free, function(lbl) {
    fit <- cov_fits[[lbl]]
    if (is.null(fit)) return(c("---", ""))
    xnames <- colnames(.as_transition_design(X_list[[cov_xmat[[lbl]]]])$X12)
    j <- match(index, xnames)
    if (is.na(j)) return(c("---", ""))
    active <- attr(.as_transition_design(X_list[[cov_xmat[[lbl]]]])$X12,
                   "entry_active")
    if (is.null(active)) active <- rep(TRUE, length(xnames))
    if (block == "beta0" && !active[j]) return(c("0.0000", "[fixed]"))
    value <- fit$params[[block]][j]
    .fmt_param(value, .get_se(cov_se[[lbl]], paste0(block, "_", index)))
  })
  list(est = c(row_label, vapply(cells, `[[`, character(1L), 1L)),
       se = c("", vapply(cells, `[[`, character(1L), 2L)))
}

all_coef_names <- unique(unlist(lapply(cov_labels_free, function(lbl) {
  colnames(.as_transition_design(X_list[[cov_xmat[[lbl]]]])$X12)
})))
coef_labels <- function(nm) if (nm %in% names(.cov_label_map)) .cov_label_map[[nm]] else nm

appendix_rows <- list(
  list(header = "Entry equation: $\\Phi(X_{it}\\beta_0)$", rows =
         lapply(all_coef_names, function(nm) .raw_parameter_row("beta0", nm, coef_labels(nm)))),
  list(header = "Employment-persistence equation: $\\Phi(X_{it}\\beta_1)$", rows =
         lapply(all_coef_names, function(nm) .raw_parameter_row("beta1", nm, coef_labels(nm)))),
  list(header = "Distribution of predicted entry hazards", rows = list(
    .cov_plain_row("entry_p10", "10th percentile", cov_labels_free),
    .cov_plain_row("entry_median", "Median", cov_labels_free),
    .cov_plain_row("entry_p90", "90th percentile", cov_labels_free)
  )),
  list(header = "Distribution of predicted exit hazards", rows = list(
    .cov_plain_row("exit_p10", "10th percentile", cov_labels_free),
    .cov_plain_row("exit_median", "Median", cov_labels_free),
    .cov_plain_row("exit_p90", "90th percentile", cov_labels_free)
  )),
  list(header = "Set 3 discrete effects on the exit probability", rows = list(
    .cov_plain_row("contract_exit_effect", "Permanent contract", cov_labels_free),
    .cov_plain_row("informal_exit_effect", "Informal sector", cov_labels_free)
  ))
)

.write_latex_table(
  col_headers = cov_col_headers_panel,
  row_data = appendix_rows,
  caption = "Covariate-model coefficients and predicted-hazard distributions",
  label = "tab:cov_coefficients_appendix",
  path = file.path(tables_ext_dir, "table_cov_coefficients_appendix.tex"),
  sub_headers = cov_sub_headers_panel,
  note = paste0("Raw probit coefficients are reported. Contract type and sector are ",
                "origin-wave characteristics and may change between transitions; their ",
                "entry coefficients are fixed at zero because these attributes are not ",
                "defined for non-employed origin states. Hazard percentiles use survey ",
                "and posterior risk-set weights. SE are bootstrap estimates when ",
                "available and otherwise use an individual-level survey-weighted ",
                "sandwich/delta approximation; quantile SE are omitted.")
)

cat("\n--- Main Table 4 (point estimates, percent) ---\n")
screen_table4 <- data.frame(
  model = cov_labels_free,
  entry = 100 * vapply(cov_implied[cov_labels_free], `[[`, numeric(1L), "mean_entry_rate"),
  exit = 100 * vapply(cov_implied[cov_labels_free], `[[`, numeric(1L), "mean_exit_rate"),
  employed_risk_share = 100 * vapply(cov_implied[cov_labels_free], `[[`, numeric(1L), "mean_employment_rate"),
  total_churn = 100 * vapply(cov_implied[cov_labels_free], `[[`, numeric(1L), "total_churn_flow"),
  pi = 100 * vapply(cov_implied[cov_labels_free], function(x) x$pi, numeric(1L))
)
print(screen_table4, row.names = FALSE, digits = 4)

cat("\n--- Appendix hazard distribution (percent) ---\n")
screen_hazards <- data.frame(
  model = cov_labels_free,
  entry_p10 = 100 * vapply(cov_implied[cov_labels_free], `[[`, numeric(1L), "entry_p10"),
  entry_p50 = 100 * vapply(cov_implied[cov_labels_free], `[[`, numeric(1L), "entry_median"),
  entry_p90 = 100 * vapply(cov_implied[cov_labels_free], `[[`, numeric(1L), "entry_p90"),
  exit_p10 = 100 * vapply(cov_implied[cov_labels_free], `[[`, numeric(1L), "exit_p10"),
  exit_p50 = 100 * vapply(cov_implied[cov_labels_free], `[[`, numeric(1L), "exit_median"),
  exit_p90 = 100 * vapply(cov_implied[cov_labels_free], `[[`, numeric(1L), "exit_p90")
)
print(screen_hazards, row.names = FALSE, digits = 4)

# ------------------------------------------------------------------------------
# 4. Table: AMEs — percentage points, proper labels, stars
# ------------------------------------------------------------------------------

cat("\n--- Building AME table ---\n")

# Focus on symmetric + stationary (main specification)
ame_labels <- c("cov_s1_sym_stat", "cov_s2_sym_stat", "cov_s3_sym_stat")

ame_col_headers <- c("", "(1)", "(2)", "(3)", "(4)", "(5)", "(6)")
ame_sub_headers <- c(
  paste(c("", "\\multicolumn{3}{c}{Entry rate}",
          "\\multicolumn{3}{c}{Exit rate}"), collapse = " & ") %+% " \\\\",
  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}",
  paste(c("", "Set 1", "Set 2", "Set 3", "Set 1", "Set 2", "Set 3"),
        collapse = " & ") %+% " \\\\"
)

# Build row-by-row
all_cov_nms <- if (!is.null(cov_implied[["cov_s3_sym_stat"]])) {
  names(cov_implied[["cov_s3_sym_stat"]]$ame_entry)
} else {
  character(0L)
}

ame_rows <- lapply(all_cov_nms, function(cn) {
  .get_ame_cell_pp <- function(lbl, outcome) {
    imp <- cov_implied[[lbl]]
    if (is.null(imp)) return(c("---", "(---)"))
    ame_vec <- imp[[paste0("ame_", outcome)]]
    ame_se  <- NA_real_
    boot    <- cov_boots[[lbl]]
    if (!is.null(boot) && !is.null(boot$ame_summary)) {
      sub <- boot$ame_summary[boot$ame_summary$outcome == outcome &
                                boot$ame_summary$covariate == cn, ]
      if (nrow(sub) > 0L) ame_se <- sub$se[1L]
    }
    v <- if (!is.null(ame_vec) && cn %in% names(ame_vec)) ame_vec[cn] else NA_real_
    .fmt_pp(v, ame_se)
  }

  entry_cells <- lapply(ame_labels, .get_ame_cell_pp, outcome = "entry")
  exit_cells  <- lapply(ame_labels, .get_ame_cell_pp, outcome = "exit")

  # Proper label
  lbl <- if (cn %in% names(.cov_label_map)) .cov_label_map[[cn]] else cn

  list(
    est = c(lbl,
            sapply(entry_cells, `[[`, 1L),
            sapply(exit_cells,  `[[`, 1L)),
    se  = c("",
            sapply(entry_cells, `[[`, 2L),
            sapply(exit_cells,  `[[`, 2L))
  )
})

# Add N row
ame_rows_with_footer <- list(
  list(header = NULL, rows = ame_rows),
  list(header = NULL, rows = list(
    list(
      est = c("$N$", rep(.fmt_n(N_ext), 6L)),
      se  = rep("", 7L)
    )
  ))
)

.write_latex_table(
  col_headers = ame_col_headers,
  row_data    = ame_rows_with_footer,
  caption     = "Average Marginal Effects on entry and exit rates (p.p.): symmetric, stationary models",
  label       = "tab:cov_ame",
  path        = file.path(tables_ext_dir, "table_cov_ame.tex"),
  col_align   = "lcccccc",
  sub_headers = ame_sub_headers,
  note        = sprintf("Estimates in percentage points. AME$_k = N^{-1}\\sum_i \\phi(x_i'\\hat\\beta)\\hat\\beta_k$. Bootstrap SE in parentheses ($B=%d$). Significance: $^{**}$ $p<0.01$, $^{*}$ $p<0.05$, $^{.}$ $p<0.10$.", B)
)

# ------------------------------------------------------------------------------
# 5. Table: FMM models (reordered: no-miscl first)
# ------------------------------------------------------------------------------

cat("\n--- Building FMM table ---\n")

# Reorder: None first, then Symmetric
fmm_labels <- c("fmm_non_stat", "fmm_non_free", "fmm_sym_stat", "fmm_sym_free")
fmm_col_headers <- c("", "(1)", "(2)", "(3)", "(4)")
fmm_sub_headers <- c(
  paste(c("", "\\multicolumn{2}{c}{No miscl.}",
          "\\multicolumn{2}{c}{Symmetric}"), collapse = " & ") %+% " \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  paste(c("", "Stat.", "Free", "Stat.", "Free"), collapse = " & ") %+% " \\\\"
)
fmm_mt <- c(fmm_non_stat = "none",      fmm_non_free = "none",
            fmm_sym_stat = "symmetric", fmm_sym_free = "symmetric")

fmm_fits <- lapply(fmm_labels, function(lbl) .load_fit(lbl, "EM-baseline-ext"))
names(fmm_fits) <- fmm_labels

fmm_boots <- lapply(fmm_labels, function(lbl) .load_boot(lbl, boot_ext_dir))
names(fmm_boots) <- fmm_labels
fmm_se    <- lapply(fmm_boots, .se_map)  # named SE vectors for O(1) lookup

fmm_implied <- lapply(fmm_labels, function(lbl) {
  fit <- fmm_fits[[lbl]]
  if (is.null(fit)) return(NULL)
  tryCatch(implied_fmm(fit$params, fmm_mt[[lbl]]), error = function(e) NULL)
})
names(fmm_implied) <- fmm_labels

.fmm_pct_row <- function(qty_name, row_label) {
  cells <- lapply(fmm_labels, function(lbl) {
    imp <- fmm_implied[[lbl]]
    v   <- if (!is.null(imp)) imp[[qty_name]] else NA_real_
    if (is.null(v)) v <- NA_real_
    se  <- .get_se(fmm_se[[lbl]], qty_name)
    .fmt_pct(v, se)
  })
  list(
    est = c(row_label, sapply(cells, `[[`, 1L)),
    se  = c("",        sapply(cells, `[[`, 2L))
  )
}

.write_latex_table(
  col_headers = fmm_col_headers,
  row_data    = list(
    list(header = "Type A (stable workers)", rows = list(
      .fmm_pct_row("entry_rate_A",      "Entry rate $\\theta_0^A$ (\\%)"),
      .fmm_pct_row("exit_rate_A",       "Exit rate $1-\\theta_1^A$ (\\%)"),
      .fmm_pct_row("employment_rate_A", "Employment rate $\\alpha^A$ (\\%)")
    )),
    list(header = "Type B (unstable workers)", rows = list(
      .fmm_pct_row("entry_rate_B",      "Entry rate $\\theta_0^B$ (\\%)"),
      .fmm_pct_row("exit_rate_B",       "Exit rate $1-\\theta_1^B$ (\\%)"),
      .fmm_pct_row("employment_rate_B", "Employment rate $\\alpha^B$ (\\%)")
    )),
    list(header = "Mixture and misclassification", rows = list(
      .fmm_pct_row("phi",                  "$\\phi$ prob.~type A (\\%)"),
      .fmm_pct_row("weighted_entry_rate",  "Wtd. avg. entry rate (\\%)"),
      .fmm_pct_row("weighted_exit_rate",   "Wtd. avg. exit rate (\\%)"),
      .fmm_pct_row("pi",                   "$\\pi$ misclassification (\\%)")
    )),
    list(header = NULL, rows = list(
      list(
        est = c("Log-likelihood", vapply(fmm_labels, function(lbl) {
          fit <- fmm_fits[[lbl]]
          if (is.null(fit)) "---" else .fmt_ll(fit$loglik)
        }, character(1L))),
        se = rep("", 5L)
      ),
      list(
        est = c("$N$", rep(.fmt_n(N_ext), 4L)),
        se  = rep("", 5L)
      )
    ))
  ),
  caption     = "Finite mixture model: implied probabilities",
  label       = "tab:fmm",
  path        = file.path(tables_ext_dir, "table_fmm.tex"),
  sub_headers = fmm_sub_headers,
  note        = sprintf("Estimates in \\%%. Type A has higher $\\theta_1$ (more stable). Bootstrap SE in parentheses ($B=%d$). Significance: $^{**}$ $p<0.01$, $^{*}$ $p<0.05$, $^{.}$ $p<0.10$.", B)
)

# ------------------------------------------------------------------------------
# 6. Table: Inconsistency models
# ------------------------------------------------------------------------------

cat("\n--- Building inconsistency table ---\n")

incons_labels      <- c("incons_sym_stat", "incons_sym_free")
incons_col_headers <- c("", "(1)", "(2)")
incons_sub_headers <- c(
  paste(c("", "Stationary", "Free $\\alpha$"), collapse = " & ") %+% " \\\\"
)

incons_fits <- lapply(incons_labels, function(lbl) .load_fit(lbl, "EM-baseline-ext"))
names(incons_fits) <- incons_labels

incons_boots <- lapply(incons_labels, function(lbl) .load_boot(lbl, boot_ext_dir))
names(incons_boots) <- incons_labels
incons_se    <- lapply(incons_boots, .se_map)  # named SE vectors for O(1) lookup

df_incons_tbl <- compute_inconsistencies(df_ext)
if (nrow(df_incons_tbl) != nrow(df_ext))
  stop(sprintf(
    "compute_inconsistencies returned %d rows; expected %d (= nrow(df_ext)).",
    nrow(df_incons_tbl), nrow(df_ext)
  ))
inc_mat_tbl   <- as.matrix(df_incons_tbl[, c("Y_age_1", "Y_age_2", "Y_age_3",
                                               "Y_edu_1", "Y_edu_2", "Y_edu_3")])

incons_implied <- lapply(incons_labels, function(lbl) {
  fit <- incons_fits[[lbl]]
  if (is.null(fit)) return(NULL)
  tryCatch(implied_inconsistency(fit$params, inc_mat_tbl), error = function(e) NULL)
})
names(incons_implied) <- incons_labels

.incons_pct_row <- function(qty_name, row_label) {
  cells <- lapply(incons_labels, function(lbl) {
    imp <- incons_implied[[lbl]]
    v   <- if (!is.null(imp)) imp[[qty_name]] else NA_real_
    if (is.null(v)) v <- NA_real_
    se  <- .get_se(incons_se[[lbl]], qty_name)
    .fmt_pct(v, se)
  })
  list(
    est = c(row_label, sapply(cells, `[[`, 1L)),
    se  = c("",        sapply(cells, `[[`, 2L))
  )
}

.incons_param_row <- function(idx, row_label) {
  cells <- lapply(incons_labels, function(lbl) {
    fit <- incons_fits[[lbl]]
    if (is.null(fit) || is.null(fit$params$delta)) return(c("---", "(---)"))
    v  <- fit$params$delta[idx]
    se <- NA_real_  # delta SEs not yet bootstrapped
    .fmt_param(v, se)
  })
  list(
    est = c(row_label, sapply(cells, `[[`, 1L)),
    se  = c("",        sapply(cells, `[[`, 2L))
  )
}

.write_latex_table(
  col_headers = incons_col_headers,
  row_data    = list(
    list(header = "Transition parameters", rows = list(
      .incons_pct_row("entry_rate",      "Entry rate $\\theta_0$ (\\%)"),
      .incons_pct_row("exit_rate",       "Exit rate $1-\\theta_1$ (\\%)"),
      .incons_pct_row("employment_rate", "Employment rate $\\alpha^*$ (\\%)")
    )),
    list(header = "Misclassification (logistic link)", rows = list(
      .incons_param_row(1L, "$\\delta_0$ (intercept)"),
      .incons_param_row(2L, "$\\delta_1$ (age incons.)"),
      .incons_param_row(3L, "$\\delta_2$ (educ. incons.)"),
      .incons_pct_row("pi_base", "$\\pi^{\\text{base}}$ (\\%)"),
      .incons_pct_row("mean_pi", "Mean $\\bar{\\pi}$ (\\%)")
    )),
    list(header = NULL, rows = list(
      list(
        est = c("Log-likelihood", vapply(incons_labels, function(lbl) {
          fit <- incons_fits[[lbl]]
          if (is.null(fit)) "---" else .fmt_ll(fit$loglik)
        }, character(1L))),
        se = rep("", 3L)
      ),
      list(
        est = c("$N$", rep(.fmt_n(N_ext), 2L)),
        se  = rep("", 3L)
      )
    ))
  ),
  caption     = "Inconsistency-augmented EM model: implied probabilities",
  label       = "tab:inconsistency",
  path        = file.path(tables_ext_dir, "table_inconsistency.tex"),
  sub_headers = incons_sub_headers,
  note        = sprintf("Implied rates in \\%%. $\\delta$ SEs pending bootstrap. Bootstrap SE for other quantities ($B=%d$). Significance: $^{**}$ $p<0.01$, $^{*}$ $p<0.05$, $^{.}$ $p<0.10$.", B)
)

cat("\n\nAll tables written.\n")
