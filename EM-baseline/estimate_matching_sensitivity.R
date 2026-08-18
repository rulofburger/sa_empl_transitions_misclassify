# ==============================================================================
# Matching-rule sensitivity: unrestricted no-error and symmetric-error models
#
# Estimates the direct eight-cell MLE on df_qlfs_A/B/C.rds using the same
# ingestion and sample filters as the baseline Tables 2--3 pipeline. Initial
# employment is free, so stationarity is not imposed.
# ==============================================================================

library(here)
library(dplyr)
library(ggplot2) # ingestion constructs diagnostic ggplot objects

source(here::here("EM-baseline", "R", "source_all.R"))

RANDOM_SEED <- 1234L
N_STARTS <- 12L
set.seed(RANDOM_SEED)

panels <- data.frame(
  panel = c("A", "B", "C"),
  source_panel = paste0("df_qlfs_", c("A", "B", "C"), ".rds"),
  matching_rule = c("Baseline", "Stricter", "Strictest"),
  stringsAsFactors = FALSE
)

results_dir <- here::here("EM-baseline", "output", "results")
tables_dir <- here::here("EM-baseline", "output", "tables")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

.prepare_panel <- function(source_panel) {
  old_option <- getOption("sa_empl_transitions.qlfs_3wave_panel")
  on.exit(options(sa_empl_transitions.qlfs_3wave_panel = old_option), add = TRUE)
  options(sa_empl_transitions.qlfs_3wave_panel = source_panel)
  source(here::here("scripts", "ingest_data_3waves_SA.R"), local = environment())

  for (nm in c("y1", "y2", "y3")) {
    if (is.factor(df_qlfs[[nm]])) stop("Employment status is factor encoded: ", nm)
    values <- df_qlfs[[nm]][!is.na(df_qlfs[[nm]])]
    if (!all(values %in% c(0, 1))) stop("Non-binary employment status: ", nm)
  }

  df_qlfs |>
    filter(!is.na(y1), !is.na(y2), !is.na(y3),
           !is.na(weight), is.finite(weight), weight > 0) |>
    transmute(y1 = as.integer(y1), y2 = as.integer(y2), y3 = as.integer(y3),
              weight = as.numeric(weight)) |>
    as.data.frame()
}

.observed_start <- function(df, model_type) {
  w <- df$weight
  y1 <- df$y1; y2 <- df$y2; y3 <- df$y3
  denom0 <- sum(w * (y1 == 0)) + sum(w * (y2 == 0))
  denom1 <- sum(w * (y1 == 1)) + sum(w * (y2 == 1))
  p <- list(
    alpha = sum(w * y1) / sum(w),
    theta0 = (sum(w * (y1 == 0 & y2 == 1)) +
                sum(w * (y2 == 0 & y3 == 1))) / denom0,
    theta1 = (sum(w * (y1 == 1 & y2 == 1)) +
                sum(w * (y2 == 1 & y3 == 1))) / denom1
  )
  if (model_type == "symmetric") p$pi <- 0.025
  p
}

.make_starts <- function(df, model_type, nested = NULL) {
  seeds <- list()
  if (!is.null(nested)) seeds[[length(seeds) + 1L]] <- nested
  seeds[[length(seeds) + 1L]] <- .observed_start(df, model_type)
  seeds[[length(seeds) + 1L]] <- init_params(model_type, stationary = FALSE)
  eta <- pack_baseline_params(seeds[[1L]], model_type, stationary = FALSE)
  while (length(seeds) < N_STARTS) {
    seeds[[length(seeds) + 1L]] <- unpack_baseline_eta(
      eta + rnorm(length(eta), 0, 0.8), model_type, stationary = FALSE
    )
  }
  seeds
}

.lookup <- function(inference, quantity, field = "estimate") {
  row <- inference$summary[inference$summary$quantity == quantity, , drop = FALSE]
  if (nrow(row) != 1L) stop("Missing analytical quantity: ", quantity)
  row[[field]]
}

fits <- list()
analytical <- list()
summary_rows <- list()

for (i in seq_len(nrow(panels))) {
  meta <- panels[i, ]
  cat(sprintf("\n=== Panel %s: %s matching (%s) ===\n",
              meta$panel, tolower(meta$matching_rule), meta$source_panel))
  df_panel <- .prepare_panel(meta$source_panel)
  cells <- collapse_baseline_cells(df_panel)
  cat(sprintf("Estimation sample: N=%s; sum(weights)=%.3f\n",
              format(cells$n, big.mark = ","), cells$weight_sum))

  none_fit <- fit_baseline_mle(
    df_panel, model_type = "none", stationary = FALSE,
    starts = .make_starts(df_panel, "none"), compute_gamma = FALSE,
    verbose = 1L, source_panel = meta$source_panel
  )
  sym_nested <- none_fit$params
  sym_nested$pi <- 0.01
  sym_fit <- fit_baseline_mle(
    df_panel, model_type = "symmetric", stationary = FALSE,
    starts = .make_starts(df_panel, "symmetric", sym_nested),
    compute_gamma = FALSE, verbose = 1L, source_panel = meta$source_panel
  )

  if (!none_fit$converged || !sym_fit$converged)
    stop("Convergence diagnostics failed for panel ", meta$panel)
  if (sym_fit$loglik < none_fit$loglik - 1e-4)
    stop("Likelihood nesting check failed for panel ", meta$panel)

  panel_fits <- list(none_free = none_fit, sym_free = sym_fit)
  panel_inference <- lapply(panel_fits, analytical_se_baseline, df = df_panel)
  fits[[meta$panel]] <- panel_fits
  analytical[[meta$panel]] <- panel_inference

  for (model in names(panel_fits)) {
    fit <- panel_fits[[model]]
    inf <- panel_inference[[model]]
    imp <- implied_baseline(fit$params, fit$model_type)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      panel = meta$panel,
      matching_rule = meta$matching_rule,
      source_panel = meta$source_panel,
      model = if (fit$model_type == "none") "No error" else "Symmetric error",
      model_type = fit$model_type,
      stationary = FALSE,
      N = fit$sample$n,
      weight_sum = fit$sample$weight_sum,
      entry = imp$entry_rate,
      entry_se = .lookup(inf, "entry_rate", "se"),
      exit = imp$exit_rate,
      exit_se = .lookup(inf, "exit_rate", "se"),
      initial_employment = fit$params$alpha,
      initial_employment_se = .lookup(inf, "alpha", "se"),
      long_run_employment = imp$employment_rate,
      long_run_employment_se = .lookup(inf, "employment_rate", "se"),
      pi = fit$params$pi %||% NA_real_,
      pi_se = if (fit$model_type == "symmetric") .lookup(inf, "pi", "se") else NA_real_,
      loglik = fit$loglik,
      optimizer_code = fit$diagnostics$optimizer_code,
      max_abs_score = fit$diagnostics$max_abs_score,
      information_min_eigenvalue = fit$diagnostics$information_min_eigenvalue,
      sample_signature = fit$sample$signature,
      stringsAsFactors = FALSE
    )
  }
  rm(df_panel)
  gc(verbose = FALSE)
}

results <- do.call(rbind, summary_rows)
rownames(results) <- NULL

# Save reproducible intermediate objects; these paths are intentionally ignored.
for (panel in panels$panel) {
  for (model in names(fits[[panel]])) {
    stem <- paste0("matching_", tolower(panel), "_", model)
    saveRDS(fits[[panel]][[model]], file.path(results_dir, paste0("fit_", stem, ".rds")))
    saveRDS(analytical[[panel]][[model]],
            file.path(results_dir, paste0("analytical_se_", stem, ".rds")))
  }
}
write.csv(results, file.path(results_dir, "matching_sensitivity_summary.csv"),
          row.names = FALSE)

.cell <- function(est, se, scale = 100, digits = 2L) {
  if (is.na(est)) return(c("---", ""))
  c(sprintf(paste0("%.", digits, "f"), scale * est),
    sprintf(paste0("(%.", digits, "f)") , scale * se))
}
.ncell <- function(x) formatC(x, format = "d", big.mark = ",")

ordered <- results[order(match(results$panel, panels$panel),
                         match(results$model_type, c("none", "symmetric"))), ]
headers <- paste0(ordered$panel, ": ", ordered$model)
quantities <- list(
  c("entry", "entry_se", "Entry rate (\\%)"),
  c("exit", "exit_se", "Exit rate (\\%)"),
  c("initial_employment", "initial_employment_se", "Initial employment $\\alpha$ (\\%)"),
  c("long_run_employment", "long_run_employment_se", "Long-run employment (\\%)"),
  c("pi", "pi_se", "Misclassification $\\pi$ (\\%)")
)

lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Sensitivity to alternative panel-matching rules}",
  "\\label{table_matching_implied}",
  "\\begin{tabular}{lcccccc}",
  "\\toprule",
  " & \\multicolumn{2}{c}{Baseline (A)} & \\multicolumn{2}{c}{Stricter (B)} & \\multicolumn{2}{c}{Strictest (C)} \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}",
  paste(c("", rep(c("No error", "Symmetric error"), 3L)), collapse = " & ") |>
    paste0(" \\\\"),
  "\\midrule"
)
for (q in quantities) {
  cells <- lapply(seq_len(nrow(ordered)), function(j) .cell(ordered[[q[1L]]][j], ordered[[q[2L]]][j]))
  lines <- c(lines,
             paste(c(q[3L], vapply(cells, `[[`, character(1L), 1L)), collapse = " & ") |> paste0(" \\\\"))
  if (any(vapply(cells, `[[`, character(1L), 2L) != "")) {
    lines <- c(lines,
               paste(c("", vapply(cells, `[[`, character(1L), 2L)), collapse = " & ") |> paste0(" \\\\"))
  }
}
lines <- c(
  lines,
  paste(c("Log-likelihood (millions)", sprintf("%.3f", ordered$loglik / 1e6)), collapse = " & ") |> paste0(" \\\\"),
  paste(c("$N$", vapply(ordered$N, .ncell, character(1L))), collapse = " & ") |> paste0(" \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{0.98\\linewidth}",
  paste0("\\footnotesize \\textit{Note:} All models estimate the initial employment probability freely; ",
         "stationarity is not imposed. Long-run employment is $\\theta_0/(\\theta_0+1-\\theta_1)$. ",
         "Parentheses contain individual-level survey-weighted sandwich/delta standard errors. ",
         "Panels A, B, and C apply successively stricter matching rules and use otherwise identical sample filters."),
  "\\end{minipage}",
  "\\end{table}"
)
table_path <- file.path(tables_dir, "table_matching_implied.tex")
writeLines(lines, table_path)

screen <- ordered[, c("panel", "matching_rule", "model", "N", "entry", "entry_se",
                      "exit", "exit_se", "initial_employment", "initial_employment_se",
                      "long_run_employment", "long_run_employment_se", "pi", "pi_se",
                      "loglik", "max_abs_score", "information_min_eigenvalue")]
for (nm in c("entry", "entry_se", "exit", "exit_se", "initial_employment",
             "initial_employment_se", "long_run_employment", "long_run_employment_se",
             "pi", "pi_se")) screen[[nm]] <- 100 * screen[[nm]]

cat("\n=== Matching-rule sensitivity estimates (rates and SEs in percent) ===\n")
print(screen, row.names = FALSE, digits = 6)
cat("\nWritten table: ", table_path, "\n", sep = "")
cat("Written machine-readable summary: ",
    file.path(results_dir, "matching_sensitivity_summary.csv"), "\n", sep = "")
