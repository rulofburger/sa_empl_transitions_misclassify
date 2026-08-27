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
#   Paper Table 1: Baseline free-initial specifications   (table_baseline_implied.tex)
#   Appendix Table A1: All baseline specifications        (table_baseline_params.tex)
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
table3_only <- identical(Sys.getenv("BUILD_TABLE3_ONLY"), "1")

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
cv_set4 <- prepare_covariate_matrix(df_ext, covariate_set = 4L)
X_list  <- list(X1 = cv_set1$X, X2 = cv_set2$X,
                X3 = cv_set3$X_transition, X4 = cv_set4$X_transition)
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
  log_tenure     = "$\\log(1+\\mathrm{tenure})$",
  log_time_since_work = "$\\log(1+\\mathrm{time\ since\ work})$",
  never_worked   = "Never worked",
  tenure_missing = "Tenure missing",
  timegap_missing = "Time since work missing",
  contracttype_1 = "Permanent contract",
  informal_sector = "Informal sector"
)

# ------------------------------------------------------------------------------
# 2. Table 1 and Appendix Table A1: baseline models
# ------------------------------------------------------------------------------

if (!table3_only) {
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

# Helper: parameter row on the raw scale.
.bl_param_row <- function(param_name, row_label, labels = baseline_labels) {
  cells <- lapply(labels, function(lbl) {
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

# Implied probability row (as percent).
.bl_implied_row <- function(qty_name, row_label, labels = baseline_labels) {
  cells <- lapply(labels, function(lbl) {
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

# Direction-specific error rates without collapsing the asymmetric model.
.bl_error_row <- function(direction = c("positive", "negative"), row_label,
                          labels = baseline_labels) {
  direction <- match.arg(direction)
  cells <- lapply(labels, function(lbl) {
    type <- model_types_bl[[lbl]]
    if (type == "none") return(c("0.00", ""))
    quantity <- if (type == "symmetric") "pi" else
      if (direction == "positive") "pi0" else "pi1"
    .fmt_pct_plain(bl_implied[[lbl]][[quantity]],
                   .get_se(bl_se[[lbl]], quantity))
  })
  list(est = c(row_label, sapply(cells, `[[`, 1L)),
       se = c("", sapply(cells, `[[`, 2L)))
}

.fmt_lr <- function(test) formatC(test$statistic, format = "f", digits = 2L,
                                  big.mark = ",")
.fmt_p_value <- function(test) {
  if (test$p_value < 0.001) "$<0.001$" else
    formatC(test$p_value, format = "f", digits = 3L)
}

lr_misclassification <- baseline_lr_test(
  bl_fits$none_free, bl_fits$sym_free, boundary = TRUE)
lr_asymmetry <- baseline_lr_test(bl_fits$sym_free, bl_fits$asym_free)
lr_stationarity <- list(
  none = baseline_lr_test(bl_fits$none_stat, bl_fits$none_free),
  symmetric = baseline_lr_test(bl_fits$sym_stat, bl_fits$sym_free),
  asymmetric = baseline_lr_test(bl_fits$asym_stat, bl_fits$asym_free))

# Main Table 1: free-initial-condition specifications only.
baseline_main_labels <- c("none_free", "sym_free", "asym_free")
baseline_main_headers <- c("", "(1)", "(2)", "(3)")
baseline_main_sub_headers <- c(
  paste(c("", "No miscl.", "Symmetric", "Asymmetric"),
        collapse = " & ") %+% " \\\\")
implied_rows_bl <- list(
  list(header = NULL, rows = list(
    .bl_implied_row("entry_rate", "Entry rate (\\%)", baseline_main_labels),
    .bl_implied_row("exit_rate", "Exit rate (\\%)", baseline_main_labels),
    .bl_error_row("positive", "False-positive rate (\\%)", baseline_main_labels),
    .bl_error_row("negative", "False-negative rate (\\%)", baseline_main_labels),
    .bl_implied_row("employment_rate", "Employment rate (\\%)",
                    baseline_main_labels)
  )),
  list(header = "Likelihood-ratio tests", rows = list(
    list(est = c("$LR$: no misclassification", "---",
                 .fmt_lr(lr_misclassification), "---"), se = rep("", 4L)),
    list(est = c("Boundary $p$-value", "---",
                 .fmt_p_value(lr_misclassification), "---"), se = rep("", 4L)),
    list(est = c("$LR$: symmetric errors", "---", "---",
                 .fmt_lr(lr_asymmetry)), se = rep("", 4L)),
    list(est = c("$p$-value", "---", "---", .fmt_p_value(lr_asymmetry)),
         se = rep("", 4L))
  )),
  list(header = NULL, rows = list(
    list(est = c("Log-likelihood", vapply(baseline_main_labels, function(lbl) {
      fit <- bl_fits[[lbl]]
      if (is.null(fit)) "---" else .fmt_ll(fit$loglik, digits = 3L)
    }, character(1L))), se = rep("", 4L)),
    list(est = c("$N$", rep(.fmt_n(N_baseline), 3L)), se = rep("", 4L))
  ))
)

.write_latex_table(
  col_headers = baseline_main_headers,
  row_data = implied_rows_bl,
  caption = "Baseline latent-state model",
  label = "tab:baseline_implied",
  path = file.path(tables_baseline_dir, "table_baseline_implied.tex"),
  sub_headers = baseline_main_sub_headers,
  note = paste0("All specifications estimate the initial employment probability ",
                "freely. Rates are percentages. Symmetric error imposes equal ",
                "false-positive and false-negative probabilities; the asymmetric ",
                "model reports them separately. Analytical SE in parentheses use an ",
                "individual-level survey-weighted sandwich covariance matrix and ",
                "the delta method. Employment is the model-implied steady-state ",
                "share. LR statistics use weights normalized to sum to $N$. The ",
                "test of no misclassification uses the 50:50 boundary mixture; ",
                "the asymmetry test uses $\\chi^2_1$. The calculation does not ",
                "incorporate strata or PSU clustering.")
)

# Appendix Table A1: all stationary and free specifications, raw parameters,
# implied rates, and tests of the stationarity restriction.
param_rows_bl <- list(
  list(header = "Transition parameters", rows = list(
    .bl_param_row("theta0", "$\\theta_0$ (entry)"),
    .bl_param_row("theta1", "$\\theta_1$ (persistence)"),
    .bl_param_row("alpha",  "$\\alpha$ (initial empl.)")
  )),
  list(header = "Misclassification parameters", rows = list(
    .bl_param_row("pi",  "$\\pi$ (symmetric)"),
    .bl_param_row("pi0", "$\\pi_0$ (false positive)"),
    .bl_param_row("pi1", "$\\pi_1$ (false negative)")
  )),
  list(header = "Implied probabilities", rows = list(
    .bl_implied_row("entry_rate", "Entry rate (\\%)"),
    .bl_implied_row("exit_rate", "Exit rate (\\%)"),
    .bl_error_row("positive", "False-positive rate (\\%)"),
    .bl_error_row("negative", "False-negative rate (\\%)"),
    .bl_implied_row("employment_rate", "Employment rate (\\%)")
  )),
  list(header = "Likelihood-ratio tests of stationarity", rows = list(
    list(est = c("$LR$: stationarity", "---", .fmt_lr(lr_stationarity$none),
      "---", .fmt_lr(lr_stationarity$symmetric), "---",
      .fmt_lr(lr_stationarity$asymmetric)), se = rep("", 7L)),
    list(est = c("$p$-value", "---", .fmt_p_value(lr_stationarity$none),
      "---", .fmt_p_value(lr_stationarity$symmetric), "---",
      .fmt_p_value(lr_stationarity$asymmetric)), se = rep("", 7L))
  )),
  list(header = NULL, rows = list(
    list(est = c("Log-likelihood", vapply(baseline_labels, function(lbl)
      .fmt_ll(bl_fits[[lbl]]$loglik, digits = 3L), character(1L))),
      se = rep("", 7L)),
    list(est = c("$N$", rep(.fmt_n(N_baseline), 6L)), se = rep("", 7L))
  ))
)

.write_latex_table(
  col_headers = baseline_col_headers,
  row_data = param_rows_bl,
  caption = "Baseline latent-state model: all specifications",
  label = "tab:baseline_params",
  path = file.path(tables_baseline_dir, "table_baseline_params.tex"),
  sub_headers = baseline_sub_headers,
  note = paste0("Rates are percentages. Analytical SE in parentheses use an ",
                "individual-level survey-weighted sandwich covariance matrix ",
                "and the delta method. $\\alpha$ is initial employment; in ",
                "stationary columns it is implied by the transition parameters. ",
                "Symmetric error imposes equal false-positive and false-negative ",
                "rates. Stationarity LR statistics appear in the corresponding ",
                "free-initial-condition column, use weights normalized to sum to ",
                "$N$, and have a $\\chi^2_1$ reference. The calculation does not ",
                "incorporate strata or PSU clustering.")
)

# ------------------------------------------------------------------------------
# 3. Covariate models: implied probabilities (split into stat/free panels)
# ------------------------------------------------------------------------------

cat("\n--- Building covariate tables ---\n")

# Sets 3--4 have transition-varying duration/job measures, so stationarity is
# not imposed on it. The main table consistently compares free-alpha models.
cov_labels_stat <- c("cov_s1_non_stat", "cov_s2_non_stat",
                     "cov_s1_sym_stat", "cov_s2_sym_stat")
cov_labels_free <- c("cov_s1_non_free", "cov_s2_non_free", "cov_s3_non_free", "cov_s4_non_free",
                     "cov_s1_sym_free", "cov_s2_sym_free", "cov_s3_sym_free", "cov_s4_sym_free")
cov_labels_all  <- unique(c(cov_labels_stat, cov_labels_free))

cov_col_headers_panel <- c("", paste0("(", 1:8, ")"))
cov_sub_headers_panel <- c(
  paste(c("", "\\multicolumn{4}{c}{No miscl.}",
          "\\multicolumn{4}{c}{Symmetric}"), collapse = " & ") %+% " \\\\",
  "\\cmidrule(lr){2-5} \\cmidrule(lr){6-9}",
  paste(c("", rep(paste("Set", 1:4), 2)),
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
  if (grepl("_s[34]_", lbl)) return(.se_map(cov_analytic[[lbl]]))
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
  cov_s3_non_stat = "none",      cov_s3_non_free = "none",
  cov_s4_sym_free = "symmetric", cov_s4_non_free = "none"
)
cov_xmat <- c(
  cov_s1_sym_stat = "X1", cov_s1_sym_free = "X1",
  cov_s1_non_stat = "X1", cov_s1_non_free = "X1",
  cov_s2_sym_stat = "X2", cov_s2_sym_free = "X2",
  cov_s2_non_stat = "X2", cov_s2_non_free = "X2",
  cov_s3_sym_stat = "X3", cov_s3_sym_free = "X3",
  cov_s3_non_stat = "X3", cov_s3_non_free = "X3",
  cov_s4_sym_free = "X4", cov_s4_non_free = "X4"
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
    if (grepl("_s3_", lbl)) return(3L)
    return(4L)
  }
  sets <- vapply(panel_labels, set_of, integer(1L))

  # Indicator: Yes if set includes covariate group
  # Set 1: Age, Age^2, Educ
  # Set 2: Set 1 + Race, Gender
  # Set 3 adds state-specific duration; Set 4 adds job characteristics.
  make_row <- function(label, min_set) {
    cells <- ifelse(sets >= min_set, "Yes", "No")
    list(est = c(label, cells), se = rep("", length(panel_labels) + 1L))
  }
  list(
    make_row("Age, Age$^2$, Educ.", 1L),
    make_row("+ Race, Gender", 2L),
    make_row("+ Tenure (exit), time since work (entry)", 3L),
    make_row("+ Contract type, informal sector (exit only)", 4L)
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

covrel_fit_path <- here::here("EM-baseline-ext", "output", "results",
                              "fit_cov_s4_reliability_free.rds")
covrel_inf_path <- here::here("EM-baseline-ext", "output", "results",
                              "analytical_se_cov_s4_reliability_free.rds")
if (!file.exists(covrel_fit_path) || !file.exists(covrel_inf_path)) {
  stop("Run EM-baseline-ext/estimate_table4_reliability.R before building Table 4.")
}
covrel_fit <- readRDS(covrel_fit_path)
covrel_inf <- readRDS(covrel_inf_path)
covrel_summary <- covrel_inf$summary
covrel_est <- setNames(covrel_summary$estimate, covrel_summary$quantity)
covrel_se <- setNames(covrel_summary$std_error, covrel_summary$quantity)

.covrel_cell <- function(quantity, formatter = .fmt_pct_plain) {
  formatter(unname(covrel_est[[quantity]]), unname(covrel_se[[quantity]]))
}

.table4_row <- function(old_quantity, new_quantity, row_label) {
  old <- .cov_plain_row(old_quantity, row_label, cov_labels_free, .fmt_pct_plain)
  new <- .covrel_cell(new_quantity)
  list(est = c(old$est, new[1L]), se = c(old$se, new[2L]))
}

table4_col_headers <- c("", paste0("(", 1:9, ")"))
table4_sub_headers <- c(
  paste(c("", "\\multicolumn{4}{c}{No miscl.}",
          "\\multicolumn{4}{c}{Constant symmetric}",
          "\\multicolumn{1}{c}{Reliability-dependent}"), collapse = " & ") %+% " \\\\",
  "\\cmidrule(lr){2-5} \\cmidrule(lr){6-9} \\cmidrule(lr){10-10}",
  paste(c("", rep(paste("Set", 1:4), 2), "Set 4"), collapse = " & ") %+% " \\\\"
)

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
  list(header = NULL, rows = list(
    .table4_row("mean_entry_rate", "mean_entry_rate", "Entry rate (\\%)"),
    .table4_row("mean_exit_rate", "mean_exit_rate", "Exit rate (\\%)"),
    .table4_row("pi", "mean_misclassification_rate", "Misclassification rate (\\%)"),
    .table4_row("mean_employment_rate", "mean_employment_rate", "Employment rate (\\%)")
  )),
  list(header = NULL, rows = list(
    list(est = c("Log-likelihood", vapply(cov_labels_free, function(lbl) {
      fit <- cov_fits[[lbl]]
      if (is.null(fit)) "---" else .fmt_ll(fit$loglik)
    }, character(1L)), .fmt_ll(covrel_fit$loglik)), se = rep("", 10L)),
    list(est = c("$N$", rep(.fmt_n(N_ext), 9L)), se = rep("", 10L))
  )),
  list(header = "Covariates included", rows = c(
    lapply(.cov_indicator_rows(cov_labels_free), function(x) {
      x$est <- c(x$est, "Yes")
      x$se <- c(x$se, "")
      x
    }),
    list(list(est = c("Reliability variables in error equation",
                       rep("No", 8L), "Yes"), se = rep("", 10L)))
  ))
)

.write_latex_table(
  col_headers = table4_col_headers,
  row_data = main_table4_rows,
  caption = "Risk-set and survey-weighted quarterly employment transitions",
  label = "tab:cov_risk_weighted",
  path = file.path(tables_ext_dir, "table_cov_risk_weighted.tex"),
  sub_headers = table4_sub_headers,
  note = paste0("All specifications estimate the initial employment probability freely. ",
                "Rates average predicted hazards using survey weights and posterior ",
                "probabilities of belonging to the relevant origin-state risk set. ",
                "Employment is the survey-weighted posterior employed share. ",
                "SE are bootstrap estimates when available and otherwise use ",
                "an individual-level survey-weighted sandwich/delta approximation. ",
                "The approximation does not incorporate strata or PSU clustering. ",
                "Column (9) retains the Set 4 transition equations and lets the ",
                "symmetric error probability depend on age and education inconsistencies, ",
                "their severity, and the two nested matching-quality indicators used in ",
                "Table 3, columns (5)--(6).")
)

# Secondary results for the reliability-dependent error equation. Transition
# coefficients are retained in the machine-readable summary CSV; this compact
# table reports the new error slopes and the probabilities needed to interpret
# them.
.covrel_secondary_row <- function(quantity, label, percent = FALSE) {
  cell <- if (percent) .covrel_cell(quantity) else
    .fmt_param_plain(unname(covrel_est[[quantity]]), unname(covrel_se[[quantity]]))
  list(est = c(label, cell[1L]), se = c("", cell[2L]))
}
covrel_secondary_rows <- list(
  list(header = "Symmetric-error equation: $0.5\\,\\Lambda(Z_{it}\\delta)$", rows = list(
    .covrel_secondary_row("error_intercept", "Intercept"),
    .covrel_secondary_row("age_inconsistency", "Age inconsistency"),
    .covrel_secondary_row("education_inconsistency", "Education inconsistency"),
    .covrel_secondary_row("large_age_inconsistency", "Large age inconsistency"),
    .covrel_secondary_row("large_education_inconsistency", "Large education inconsistency"),
    .covrel_secondary_row("panel_B_not_C", "Panel B but not C"),
    .covrel_secondary_row("panel_A_not_B", "Panel A but not B")
  )),
  list(header = "Implied misclassification probabilities", rows = list(
    .covrel_secondary_row("pi_base", "No inconsistency; panel C (\\%)", TRUE),
    .covrel_secondary_row("pi_age_mild", "Age inconsistency (\\%)", TRUE),
    .covrel_secondary_row("pi_age_severe", "Large age inconsistency (\\%)", TRUE),
    .covrel_secondary_row("pi_education_mild", "Education inconsistency (\\%)", TRUE),
    .covrel_secondary_row("pi_education_severe", "Large education inconsistency (\\%)", TRUE),
    .covrel_secondary_row("pi_both_mild", "Both inconsistencies (\\%)", TRUE),
    .covrel_secondary_row("pi_both_both_severe", "Both large inconsistencies (\\%)", TRUE),
    .covrel_secondary_row("pi_B_not_C", "Panel B but not C (\\%)", TRUE),
    .covrel_secondary_row("pi_A_not_B", "Panel A but not B (\\%)", TRUE)
  ))
)
.write_latex_table(
  col_headers = c("", "Set 4; reliability-dependent error"),
  row_data = covrel_secondary_rows,
  caption = "Reliability-dependent misclassification estimates for Table 4",
  label = "tab:cov_reliability_appendix",
  path = file.path(tables_ext_dir, "table_cov_reliability_appendix.tex"),
  note = paste0("Parentheses contain analytical individual-level ",
                "survey-weighted sandwich/delta-method standard errors. ",
                "Each implied rate changes the named indicator(s) from the ",
                "panel-C, internally consistent reference case while holding ",
                "the remaining indicators at zero.")
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
    active1 <- attr(.as_transition_design(X_list[[cov_xmat[[lbl]]]])$X12,
                    "persistence_active")
    if (is.null(active1)) active1 <- rep(TRUE, length(xnames))
    if ((block == "beta0" && !active[j]) || (block == "beta1" && !active1[j]))
      return(c("0.0000", "[fixed]"))
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
  list(header = "Other derived model quantities", rows = list(
    .cov_plain_row("alpha", "Initial employment rate (\\%)", cov_labels_free),
    .cov_plain_row("total_churn_flow", "Total employment churn (\\%)", cov_labels_free)
  )),
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
  list(header = "Set 4 discrete effects on the exit probability", rows = list(
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
  note = paste0("Raw probit coefficients are reported. Tenure enters persistence only; ",
                "time since work and never worked enter entry only. Contract type and sector are ",
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
cat("\nReliability-dependent Set 4 extension:\n")
print(data.frame(
  model = "cov_s4_reliability_free",
  entry = 100 * covrel_est[["mean_entry_rate"]],
  exit = 100 * covrel_est[["mean_exit_rate"]],
  employed_risk_share = 100 * covrel_est[["mean_employment_rate"]],
  pi = 100 * covrel_est[["mean_misclassification_rate"]],
  loglik = covrel_fit$loglik,
  N = covrel_fit$n_obs
), row.names = FALSE, digits = 4)

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
ame_labels <- c("cov_s1_sym_stat", "cov_s2_sym_stat")

ame_col_headers <- c("", "(1)", "(2)", "(3)", "(4)")
ame_sub_headers <- c(
  paste(c("", "\\multicolumn{2}{c}{Entry rate}",
          "\\multicolumn{2}{c}{Exit rate}"), collapse = " & ") %+% " \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  paste(c("", "Set 1", "Set 2", "Set 1", "Set 2"),
        collapse = " & ") %+% " \\\\"
)

# Build row-by-row
all_cov_nms <- if (!is.null(cov_implied[["cov_s2_sym_stat"]])) {
  names(cov_implied[["cov_s2_sym_stat"]]$ame_entry)
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
      est = c("$N$", rep(.fmt_n(N_ext), 4L)),
      se  = rep("", 5L)
    )
  ))
)

.write_latex_table(
  col_headers = ame_col_headers,
  row_data    = ame_rows_with_footer,
  caption     = "Average Marginal Effects on entry and exit rates (p.p.): symmetric, stationary models",
  label       = "tab:cov_ame",
  path        = file.path(tables_ext_dir, "table_cov_ame.tex"),
  col_align   = "lcccc",
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
}

# ------------------------------------------------------------------------------
# 6. Table: Inconsistency models
# ------------------------------------------------------------------------------

cat("\n--- Building inconsistency table ---\n")
audit_path <- here::here("EM-baseline-ext", "output", "results",
                         "fit_table6_inconsistency_audit.rds")
if (!file.exists(audit_path))
  stop("Missing Table 3 estimates. Run EM-baseline-ext/replicate_table6.R first.")
incons_audit <- readRDS(audit_path)

.audit_lookup <- function(tbl, quantity, estimate_col = "estimate",
                          se_col = "std_error") {
  i <- match(quantity, tbl$quantity)
  if (is.na(i)) return(c(NA_real_, NA_real_))
  c(tbl[[estimate_col]][i], tbl[[se_col]][i])
}
.audit_row <- function(tbl, quantity, label, probability = TRUE,
                       estimate_col = "estimate", se_col = "std_error") {
  x <- .audit_lookup(tbl, quantity, estimate_col, se_col)
  cell <- if (probability) .fmt_pct(x[1L], x[2L]) else .fmt_param(x[1L], x[2L])
  list(est = c(label, cell[1L]), se = c("", cell[2L]))
}

main_tbl <- incons_audit$table
increment_tbl <- incons_audit$extent_comparison_table
matching_tbl <- incons_audit$matching_comparison_table
transition_cov_path <- here::here("EM-baseline-ext", "output", "results",
                                  "fit_table3_transition_covariates.rds")
if (!file.exists(transition_cov_path))
  stop("Missing Table 3 transition-covariate estimates. Run EM-baseline-ext/estimate_table3_transition_covariates.R first.")
transition_cov <- readRDS(transition_cov_path)
cov4_tbl <- transition_cov$column4_inference$summary
cov5_tbl <- transition_cov$column5_inference$summary
if (is.null(transition_cov$column5_stationary_inference))
  stop("Missing the conditionally stationary Table 3 preferred model. Re-run EM-baseline-ext/estimate_table3_transition_covariates.R.")
cov5_stat_tbl <- transition_cov$column5_stationary_inference$summary
.model_cell <- function(tbl, quantity, specification, probability = TRUE) {
  x <- .audit_lookup(tbl, quantity,
    paste0("estimate_", specification), paste0("std_error_", specification))
  if (probability) .fmt_pct(x[1L], x[2L]) else .fmt_param(x[1L], x[2L])
}
.six_model_row <- function(basic_quantity, increment_quantity, label,
                           probability = TRUE,
                           basic_stationary_quantity = basic_quantity,
                           matching_quantity = increment_quantity) {
  cells <- list(
    .model_cell(main_tbl, basic_stationary_quantity, "stationary", probability),
    .model_cell(main_tbl, basic_quantity, "free", probability),
    .model_cell(increment_tbl, increment_quantity, "stationary", probability),
    .model_cell(increment_tbl, increment_quantity, "free", probability),
    .model_cell(matching_tbl, matching_quantity, "stationary", probability),
    .model_cell(matching_tbl, matching_quantity, "free", probability))
  list(est = c(label, vapply(cells, `[[`, character(1L), 1L)),
       se = c("", vapply(cells, `[[`, character(1L), 2L)))
}
.three_free_row <- function(basic_quantity, increment_quantity, label,
                            probability = TRUE,
                            matching_quantity = increment_quantity) {
  cells <- list(
    .model_cell(main_tbl, basic_quantity, "free", probability),
    .model_cell(increment_tbl, increment_quantity, "free", probability),
    .model_cell(matching_tbl, matching_quantity, "free", probability))
  list(est = c(label, vapply(cells, `[[`, character(1L), 1L)),
       se = c("", vapply(cells, `[[`, character(1L), 2L)))
}
.cov_cell <- function(tbl, quantity, probability = TRUE) {
  x <- .audit_lookup(tbl, quantity)
  if (probability) .fmt_pct(x[1L], x[2L]) else .fmt_param(x[1L], x[2L])
}
.two_cov_row <- function(quantity, label, probability = FALSE) {
  c4 <- .cov_cell(cov4_tbl, quantity, probability)
  c5 <- .cov_cell(cov5_tbl, quantity, probability)
  list(est = c(label, c4[1L], c5[1L]), se = c("", c4[2L], c5[2L]))
}
.cov5_only_row <- function(quantity, label, probability = FALSE) {
  c5 <- .cov_cell(cov5_tbl, quantity, probability)
  list(est = c(label, "---", c5[1L]), se = c("", "", c5[2L]))
}
.five_main_row <- function(basic_quantity, increment_quantity, cov_quantity,
                           label, matching_quantity = increment_quantity) {
  base <- .three_free_row(basic_quantity, increment_quantity, label,
                          matching_quantity = matching_quantity)
  c4 <- .cov_cell(cov4_tbl, cov_quantity)
  c5 <- .cov_cell(cov5_tbl, cov_quantity)
  list(est = c(base$est, c4[1L], c5[1L]),
       se = c(base$se, c4[2L], c5[2L]))
}
.three_selected_main_row <- function(basic_quantity, matching_quantity,
                                     final_quantity, label) {
  cells <- list(
    .model_cell(main_tbl, basic_quantity, "free"),
    .model_cell(matching_tbl, matching_quantity, "free"),
    .cov_cell(cov5_tbl, final_quantity))
  list(est = c(label, vapply(cells, `[[`, character(1L), 1L)),
       se = c("", vapply(cells, `[[`, character(1L), 2L)))
}
.eight_model_row <- function(basic_quantity, increment_quantity, label,
                             final_quantity, probability = TRUE,
                             basic_stationary_quantity = basic_quantity,
                             matching_quantity = increment_quantity) {
  first <- .six_model_row(
    basic_quantity, increment_quantity, label, probability,
    basic_stationary_quantity = basic_stationary_quantity,
    matching_quantity = matching_quantity)
  final <- list(
    .cov_cell(cov5_stat_tbl, final_quantity, probability),
    .cov_cell(cov5_tbl, final_quantity, probability))
  list(est = c(first$est, vapply(final, `[[`, character(1L), 1L)),
       se = c(first$se, vapply(final, `[[`, character(1L), 2L)))
}
.final_pair_row <- function(quantity, label, probability = FALSE,
                            prefix_est = rep("---", 6L),
                            prefix_se = rep("", 6L)) {
  final <- list(.cov_cell(cov5_stat_tbl, quantity, probability),
                .cov_cell(cov5_tbl, quantity, probability))
  list(est = c(label, prefix_est,
               vapply(final, `[[`, character(1L), 1L)),
       se = c("", prefix_se,
              vapply(final, `[[`, character(1L), 2L)))
}
.final_two_row <- function(quantity, label, probability = FALSE) {
  final <- list(.cov_cell(cov5_stat_tbl, quantity, probability),
                .cov_cell(cov5_tbl, quantity, probability))
  list(est = c(label, vapply(final, `[[`, character(1L), 1L)),
       se = c("", vapply(final, `[[`, character(1L), 2L)))
}
.increment_effect_row <- function(quantity, label, probability = TRUE) {
  stat <- .model_cell(increment_tbl, quantity, "stationary", probability)
  free <- .model_cell(increment_tbl, quantity, "free", probability)
  match_stat <- .model_cell(matching_tbl, quantity, "stationary", probability)
  match_free <- .model_cell(matching_tbl, quantity, "free", probability)
  list(est = c(label, "0.00", "0.00", stat[1L], free[1L],
               match_stat[1L], match_free[1L]),
       se = c("", "(fixed)", "(fixed)", stat[2L], free[2L],
              match_stat[2L], match_free[2L]))
}
.matching_only_row <- function(quantity, label, probability = TRUE,
                               omitted_value = "---",
                               omitted_se = "") {
  stat <- .model_cell(matching_tbl, quantity, "stationary", probability)
  free <- .model_cell(matching_tbl, quantity, "free", probability)
  list(est = c(label, rep(omitted_value, 4L), stat[1L], free[1L]),
       se = c("", rep(omitted_se, 4L), stat[2L], free[2L]))
}
.pad_six_row <- function(row, values = c("---", "---"),
                         ses = c("", "")) {
  list(est = c(row$est, values), se = c(row$se, ses))
}
.eight_matching_only_row <- function(quantity, label, probability = TRUE,
                                     omitted_value = "---",
                                     omitted_se = "") {
  .pad_six_row(.matching_only_row(
    quantity, label, probability, omitted_value, omitted_se))
}
.eight_increment_effect_row <- function(quantity, label,
                                        probability = TRUE) {
  .pad_six_row(.increment_effect_row(quantity, label, probability),
               values = c("0.00", "0.00"), ses = c("(fixed)", "(fixed)"))
}
.a3_select_row <- function(row) {
  keep <- c(1L, 2L, 3L, 6L, 7L, 8L, 9L)
  list(est = row$est[keep], se = row$se[keep])
}
.a3_model_row <- function(...) .a3_select_row(.eight_model_row(...))
.a3_matching_only_row <- function(...)
  .a3_select_row(.eight_matching_only_row(...))
.a3_increment_effect_row <- function(...)
  .a3_select_row(.eight_increment_effect_row(...))
.a3_rows <- function(rows) lapply(rows, .a3_select_row)
.stationarity_gap_row <- function() {
  basic_free <- .model_cell(main_tbl, "stationarity_gap", "free")
  increment_free <- .model_cell(increment_tbl, "stationarity_gap", "free")
  matching_free <- .model_cell(matching_tbl, "stationarity_gap", "free")
  list(est = c("Initial minus steady employment (p.p.)",
               "0.00", basic_free[1L], "0.00", increment_free[1L],
               "0.00", matching_free[1L]),
       se = c("", "(fixed)", basic_free[2L], "(fixed)", increment_free[2L],
              "(fixed)", matching_free[2L]))
}
table8_headers <- c("", paste0("(", 1:6, ")"))
table8_sub_headers <- c(
  paste(c("", "\\multicolumn{2}{c}{\\shortstack{Inconsistency\\\\ indicators}}",
          "\\multicolumn{2}{c}{\\shortstack{+ severity and\\\\ matching-rule indicators}}",
          "\\multicolumn{2}{c}{\\shortstack{+ transition covariates\\\\ and inconsistency effects}}"),
        collapse = " & ") %+% " \\\\ ",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
  paste(c("", rep(c("Stationary", "Free $\\alpha$"), 2L),
          "Conditional stat.", "Free $\\alpha$"),
        collapse = " & ") %+% " \\\\ "
)
table8_fit_rows <- lapply(list(
  list(est = c("Log-likelihood",
               .fmt_ll(incons_audit$stationary$loglik, 3L),
               .fmt_ll(incons_audit$free$loglik, 3L),
               .fmt_ll(incons_audit$extent_stationary$loglik, 3L),
               .fmt_ll(incons_audit$extent_free$loglik, 3L),
               .fmt_ll(incons_audit$matching_stationary$loglik, 3L),
               .fmt_ll(incons_audit$matching_free$loglik, 3L),
               .fmt_ll(transition_cov$column5_stationary$loglik, 3L),
               .fmt_ll(transition_cov$column5$loglik, 3L)),
       se = rep("", 9L)),
  list(est = c("$N$", rep(.fmt_n(N_ext), 8L)), se = rep("", 9L)),
  list(est = c("Stationarity imposed", rep(c("Yes", "No"), 3L),
               "Conditional", "No"), se = rep("", 9L)),
  list(est = c("Inconsistency-severity effects", "No", "No", rep("Yes", 4L),
               "No", "No"), se = rep("", 9L)),
  list(est = c("Matching-rule indicators", rep("No", 4L), "Yes", "Yes",
               "No", "No"), se = rep("", 9L)),
  list(est = c("Set 2 transition covariates", rep("No", 6L), "Yes", "Yes"),
       se = rep("", 9L)),
  list(est = c("Inconsistency effects in transitions", rep("No", 6L), "Yes", "Yes"),
       se = rep("", 9L))
), .a3_select_row)
table3_headers <- c("", paste0("(", 1:3, ")"))
table3_sub_headers <- paste(
  c("", "\\shortstack{Inconsistency\\\\ indicators}",
    "\\shortstack{+ severity and\\\\ matching-rule indicators}",
    "\\shortstack{+ transition covariates\\\\ and inconsistency effects}"),
  collapse = " & ") %+% " \\\\ "
table3_fit_rows <- list(
  list(est = c("Log-likelihood",
               .fmt_ll(incons_audit$free$loglik, 3L),
               .fmt_ll(incons_audit$matching_free$loglik, 3L),
               .fmt_ll(transition_cov$column5$loglik, 3L)),
       se = rep("", 4L)),
  list(est = c("$N$", rep(.fmt_n(N_ext), 3L)), se = rep("", 4L)),
  list(est = c("Inconsistency-severity effects", "No", "Yes", "No"),
       se = rep("", 4L)),
  list(est = c("Matching-rule indicators", "No", "Yes", "No"),
       se = rep("", 4L)),
  list(est = c("Set 2 transition covariates", "No", "No", "Yes"),
       se = rep("", 4L)),
  list(est = c("Inconsistency effects in transitions", "No", "No", "Yes"),
       se = rep("", 4L))
)

.write_latex_table(
  col_headers = table3_headers,
  sub_headers = table3_sub_headers,
  row_data = list(
    list(header = NULL, rows = list(
      .three_selected_main_row("entry_rate", "entry_rate", "mean_entry_rate",
                               "Entry rate (\\%)"),
      .three_selected_main_row("exit_rate", "exit_rate", "mean_exit_rate",
                               "Exit rate (\\%)"),
      .three_selected_main_row("mean_pi_survey_weighted",
        "mean_pi_survey_weighted", "mean_misclassification_rate",
        "Misclassification rate (\\%)"),
      .three_selected_main_row("initial_employment", "initial_employment",
                               "initial_employment",
                               "Initial employment rate (\\%)")
    )),
    list(header = NULL, rows = table3_fit_rows)
  ),
  caption = "Reported-characteristic inconsistencies",
  label = "tab:inconsistency",
  path = file.path(tables_ext_dir, "table_inconsistency.tex"),
  note = "All specifications estimate initial employment freely. Each misclassification equation includes age, education, race, and gender inconsistency indicators plus mutually exclusive indicators for exactly two, three, or four inconsistencies (zero or one is omitted); race and gender disagreements are attributed to waves using the same adjacent-wave logic as matching rule B in Table 2. Column (2) adds age and education severity effects and the matching-rule indicators. Column (3) instead combines this error equation with Set 2 transition covariates and the four origin-wave inconsistency effects in both entry and persistence. Its entry and exit rates are posterior risk-set, survey-weighted means. Mean misclassification averages individual-wave predictions using survey weights. Analytical standard errors use the survey-weighted sandwich covariance and delta method. An inconsistency may arise from response error, matching error, or both."
)

.write_latex_table(
  col_headers = table8_headers,
  sub_headers = table8_sub_headers,
  row_data = list(
    list(header = "Implied rates", rows = .a3_rows(list(
      .eight_model_row("entry_rate", "entry_rate", "Entry rate (\\%)",
                       "mean_entry_rate"),
      .eight_model_row("exit_rate", "exit_rate", "Exit rate (\\%)",
                       "mean_exit_rate"),
      .eight_model_row("initial_employment", "initial_employment",
                       "Initial employment rate (\\%)", "initial_employment",
                       basic_stationary_quantity = "steady_employment")
    ))),
    list(header = "Employment-state transformations", rows = .a3_rows(list(
      .pad_six_row(.six_model_row("steady_employment", "steady_employment",
                                  "Steady-state employment (\\%)")),
      .pad_six_row(.stationarity_gap_row())
    ))),
    list(header = "Distribution of predicted misclassification probabilities",
         rows = .a3_rows(list(
      .eight_model_row("min_pi_survey_weighted", "min_pi_survey_weighted",
                       "Lowest (\\%)", "min_misclassification_probability"),
      .eight_model_row("median_pi_survey_weighted", "median_pi_survey_weighted",
                       "Median (\\%)", "median_misclassification_probability"),
      .eight_model_row("mean_pi_survey_weighted", "mean_pi_survey_weighted",
                       "Mean (\\%)", "mean_misclassification_rate"),
      .eight_model_row("max_pi_survey_weighted", "max_pi_survey_weighted",
                       "Highest (\\%)", "max_misclassification_probability")
    ))),
    list(header = "Logistic-link coefficients", rows = .a3_rows(list(
      .eight_model_row("delta0", "delta0", "$\\delta_0$ (intercept)",
                       "error_intercept", FALSE),
      .eight_model_row("delta_age", "delta_age",
                       "$\\delta_1$ (age inconsistency)",
                       "age_inconsistency", FALSE),
      .eight_model_row("delta_education", "delta_education",
                       "$\\delta_2$ (education inconsistency)",
                       "education_inconsistency", FALSE),
      .eight_model_row("delta_race", "delta_race",
                       "$\\delta_3$ (race inconsistency)",
                       "race_inconsistency", FALSE),
      .eight_model_row("delta_gender", "delta_gender",
                       "$\\delta_4$ (gender inconsistency)",
                       "gender_inconsistency", FALSE),
      .eight_model_row("delta_two_inconsistencies", "delta_two_inconsistencies",
                       "Exactly two inconsistencies", "two_inconsistencies", FALSE),
      .eight_model_row("delta_three_inconsistencies", "delta_three_inconsistencies",
                       "Exactly three inconsistencies", "three_inconsistencies", FALSE),
      .eight_model_row("delta_four_inconsistencies", "delta_four_inconsistencies",
                       "Exactly four inconsistencies", "four_inconsistencies", FALSE),
      .eight_increment_effect_row("delta_age_severe",
                                  "$\\delta_5$ (large age increment)", FALSE),
      .eight_increment_effect_row("delta_education_severe",
                                  "$\\delta_6$ (large education increment)", FALSE),
      .eight_matching_only_row("delta_B_not_C",
                               "$\\delta_7$ (panel B but not C)",
                               FALSE, "0.000", "(fixed)"),
      .eight_matching_only_row("delta_A_not_B",
                               "$\\delta_8$ (panel A but not B)",
                               FALSE, "0.000", "(fixed)")
    ))),
    list(header = NULL, rows = table8_fit_rows)
  ),
  caption = "Reported-characteristic inconsistencies: implied rates and error equation",
  label = "tab:inconsistency_details",
  path = file.path(tables_ext_dir, "table_inconsistency_details.tex"),
  note = "The table retains stationary and free-initial-employment versions of the same three specifications shown in Table 3. Misclassification summaries describe the survey-weighted distribution of predicted person-wave probabilities over the estimation sample; lowest and highest are the observed-support extrema. Error-equation coefficients are on the half-logit index. Exact-count indicators are mutually exclusive and omit zero or one inconsistency. In the conditionally stationary preferred model, each person's wave-1 initial employment probability is implied by their wave-1 covariate-specific entry and exit hazards; later hazards may change with origin-wave inconsistencies. Parentheses contain analytical sandwich/delta standard errors."
)

.write_latex_table(
  col_headers = c("", "Conditional stat.", "Free $\\alpha$"),
  row_data = list(
    list(header = "Entry coefficients (probit index)", rows =
      lapply(transition_cov$design$column5, function(term)
        .final_two_row(paste0("entry_", term),
                       gsub("_", " ", term, fixed = TRUE)))),
    list(header = "Persistence coefficients (probit index)", rows =
      lapply(transition_cov$design$column5, function(term)
        .final_two_row(paste0("persistence_", term),
                       gsub("_", " ", term, fixed = TRUE))))
  ),
  caption = "Transition-covariate specification: transition coefficients",
  label = "tab:inconsistency_details_coefficients",
  path = file.path(tables_ext_dir,
                   "table_inconsistency_details_coefficients.tex"),
  note = "The conditional-stationarity column derives each person's wave-1 initial employment probability from their wave-1 covariate-specific entry and exit hazards; the free column estimates one common initial probability. Parentheses contain analytical sandwich standard errors."
)

.write_latex_table(
  col_headers = c("", "Conditional stat.", "Free $\\alpha$"),
  row_data = list(
    list(header = "Misclassification coefficients (half-logit index)", rows = list(
      .final_two_row("error_intercept", "Intercept"),
      .final_two_row("age_inconsistency", "Age inconsistency"),
      .final_two_row("education_inconsistency", "Education inconsistency"),
      .final_two_row("race_inconsistency", "Race inconsistency"),
      .final_two_row("gender_inconsistency", "Gender inconsistency"),
      .final_two_row("two_inconsistencies", "Exactly two inconsistencies"),
      .final_two_row("three_inconsistencies", "Exactly three inconsistencies"),
      .final_two_row("four_inconsistencies", "Exactly four inconsistencies")
    )),
    list(header = "Initial state and diagnostics", rows = list(
      .final_two_row("initial_employment", "Mean initial employment probability"),
      list(est = c("Log-likelihood",
        .fmt_ll(transition_cov$column5_stationary$loglik, 3L),
        .fmt_ll(transition_cov$column5$loglik, 3L)), se = rep("", 3L)),
      list(est = c("Information rank",
        paste0(transition_cov$column5_stationary_inference$diagnostics$information_rank,
               "/", transition_cov$column5_stationary_inference$diagnostics$K),
        paste0(transition_cov$column5_inference$diagnostics$information_rank,
               "/", transition_cov$column5_inference$diagnostics$K)),
        se = rep("", 3L)),
      list(est = c("Information condition number",
        formatC(transition_cov$column5_stationary_inference$diagnostics$information_condition,
                format = "fg", digits = 4),
        formatC(transition_cov$column5_inference$diagnostics$information_condition,
                format = "fg", digits = 4)), se = rep("", 3L)),
      list(est = c("$N$", rep(.fmt_n(N_ext), 2L)), se = rep("", 3L))
    ))
  ),
  caption = "Transition-covariate specification: misclassification coefficients and diagnostics",
  label = "tab:inconsistency_details_diagnostics",
  path = file.path(tables_ext_dir,
                   "table_inconsistency_details_diagnostics.tex"),
  note = "The conditional-stationarity column derives each person's wave-1 initial employment probability from their wave-1 covariate-specific entry and exit hazards; the free column estimates one common initial probability. The likelihoods are not nested because only the former allows observed initial-state heterogeneity. Parentheses contain analytical sandwich/delta standard errors."
)

.write_latex_table(
  col_headers = c("", "(A)", "(B)"),
  sub_headers = paste(c("", "Set 2 transition covariates",
                        "+ inconsistency effects, free $\\alpha$"),
                      collapse = " & ") %+% " \\\\ ",
  row_data = list(
    list(header = "Origin-wave inconsistency coefficients: entry", rows = list(
      .cov5_only_row("entry_origin_age_inconsistency", "Age inconsistency"),
      .cov5_only_row("entry_origin_education_inconsistency", "Education inconsistency"),
      .cov5_only_row("entry_origin_race_inconsistency", "Race inconsistency"),
      .cov5_only_row("entry_origin_gender_inconsistency", "Gender inconsistency")
    )),
    list(header = "Origin-wave inconsistency coefficients: persistence", rows = list(
      .cov5_only_row("persistence_origin_age_inconsistency", "Age inconsistency"),
      .cov5_only_row("persistence_origin_education_inconsistency", "Education inconsistency"),
      .cov5_only_row("persistence_origin_race_inconsistency", "Race inconsistency"),
      .cov5_only_row("persistence_origin_gender_inconsistency", "Gender inconsistency")
    )),
    list(header = "Misclassification-equation coefficients", rows = list(
      .two_cov_row("error_intercept", "Intercept"),
      .two_cov_row("age_inconsistency", "Age inconsistency"),
      .two_cov_row("education_inconsistency", "Education inconsistency"),
      .two_cov_row("race_inconsistency", "Race inconsistency"),
      .two_cov_row("gender_inconsistency", "Gender inconsistency")
    )),
    list(header = "Implied misclassification probabilities", rows = list(
      .two_cov_row("pi_base", "No inconsistency (\\%)", TRUE),
      .two_cov_row("pi_age_mild", "Age inconsistency (\\%)", TRUE),
      .two_cov_row("pi_education_mild", "Education inconsistency (\\%)", TRUE),
      .two_cov_row("pi_race", "Race inconsistency (\\%)", TRUE),
      .two_cov_row("pi_gender", "Gender inconsistency (\\%)", TRUE)
    )),
    list(header = NULL, rows = list(
      list(est = c("Log-likelihood",
        .fmt_ll(transition_cov$column4$loglik, 3L),
        .fmt_ll(transition_cov$column5$loglik, 3L)), se = rep("", 3L)),
      list(est = c("Information rank",
        paste0(transition_cov$column4_inference$diagnostics$information_rank, "/",
               transition_cov$column4_inference$diagnostics$K),
        paste0(transition_cov$column5_inference$diagnostics$information_rank, "/",
               transition_cov$column5_inference$diagnostics$K)), se = rep("", 3L)),
      list(est = c("Information condition number",
        formatC(transition_cov$column4_inference$diagnostics$information_condition,
                format = "fg", digits = 4),
        formatC(transition_cov$column5_inference$diagnostics$information_condition,
                format = "fg", digits = 4)), se = rep("", 3L)),
      list(est = c("$N$", rep(.fmt_n(N_ext), 2L)), se = rep("", 3L))
    ))
  ),
  caption = "Transition-covariate specifications",
  label = "tab:inconsistency_transition_covariates",
  path = file.path(tables_ext_dir, "table_inconsistency_transition_covariates.tex"),
  note = "Transition equations use probit links; the misclassification equation uses the half-logit link. Column (B) adds origin-wave inconsistency effects to both entry and persistence. The complete coefficients for this specification under free initial employment and conditional stationarity appear in Appendix Table A3, Panel B. Analytical standard errors use the survey-weighted sandwich covariance and delta method."
)

group_tbl <- incons_audit$reliability_group_table
.write_latex_table(
  col_headers = c("", "(1)"),
  row_data = list(
    list(header = "Reliable records", rows = list(
      .audit_row(group_tbl, "entry_reliable", "Entry rate (\\%)"),
      .audit_row(group_tbl, "exit_reliable", "Exit rate (\\%)"),
      .audit_row(group_tbl, "initial_reliable", "Initial employment (\\%)")
    )),
    list(header = "Records with any inconsistency", rows = list(
      .audit_row(group_tbl, "entry_unreliable", "Entry rate (\\%)"),
      .audit_row(group_tbl, "exit_unreliable", "Exit rate (\\%)"),
      .audit_row(group_tbl, "initial_unreliable", "Initial employment (\\%)")
    )),
    list(header = "Misclassification probabilities", rows = list(
      .audit_row(group_tbl, "pi_base", "No inconsistency (\\%)"),
      .audit_row(group_tbl, "pi_age", "Age inconsistency (\\%)"),
      .audit_row(group_tbl, "age_effect", "Age increment (p.p.)"),
      .audit_row(group_tbl, "pi_education", "Education inconsistency (\\%)"),
      .audit_row(group_tbl, "education_effect", "Education increment (p.p.)"),
      .audit_row(group_tbl, "mean_pi_survey_weighted", "Survey-weighted mean (\\%)")
    )),
    list(header = NULL, rows = list(list(est = c("$N$", .fmt_n(N_ext)), se = c("", ""))))
  ),
  caption = "Robustness: true transition rates by record-reliability group",
  label = "tab:inconsistency_reliability",
  path = file.path(tables_ext_dir, "table_inconsistency_reliability.tex"),
  note = "The true transition process and initial employment probability may differ between records with and without any age or education inconsistency. Other conventions follow Table~\\ref{tab:inconsistency}."
)

extent_tbl <- incons_audit$extent_table
.write_latex_table(
  col_headers = c("", "(1)"),
  row_data = list(
    list(header = "Age inconsistencies", rows = list(
      .audit_row(extent_tbl, "pi_age_mild", "Mild inconsistency (\\%)"),
      .audit_row(extent_tbl, "pi_age_severe", "Severe inconsistency (\\%)"),
      .audit_row(extent_tbl, "age_severity_effect", "Severe--mild difference (p.p.)"),
      .audit_row(extent_tbl, "delta_age_severe", "Severe increment, logit scale", FALSE)
    )),
    list(header = "Education inconsistencies", rows = list(
      .audit_row(extent_tbl, "pi_education_mild", "Mild inconsistency (\\%)"),
      .audit_row(extent_tbl, "pi_education_severe", "Severe inconsistency (\\%)"),
      .audit_row(extent_tbl, "education_severity_effect", "Severe--mild difference (p.p.)"),
      .audit_row(extent_tbl, "delta_education_severe", "Severe increment, logit scale", FALSE)
    )),
    list(header = NULL, rows = list(list(est = c("$N$", .fmt_n(N_ext)), se = c("", ""))))
  ),
  caption = "Robustness: misclassification by inconsistency severity",
  label = "tab:inconsistency_extent",
  path = file.path(tables_ext_dir, "table_inconsistency_extent.tex"),
  note = "A severe inconsistency lies at least two units beyond the admissible zero-to-one progression between waves. Analytical standard errors are in parentheses."
)

cat("\n\nAll tables written.\n")
