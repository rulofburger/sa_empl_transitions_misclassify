# ==============================================================================
# Table 5 finite-mixture replication and identification audit
#
# The paper's Table 5 used the extension sample (complete age/education), while
# the FMM itself does not contain those covariates. We estimate both that legacy
# sample and the full baseline sample. The displayed replication uses the legacy
# sample so differences from the paper are not driven by sample composition.
# ==============================================================================

library(here)
library(dplyr)
library(ggplot2)

source(here::here("EM-baseline", "R", "source_all.R"))
source(here::here("EM-baseline-ext", "R", "source_all.R"))
source(here::here("scripts", "ingest_data_3waves_SA.R"))

RANDOM_SEED <- 1234L
N_STARTS <- 16L
RUN_BASELINE_SAMPLE <- identical(Sys.getenv("FMM_RUN_BASELINE_SAMPLE"), "1")
set.seed(RANDOM_SEED)

.as_fmm_panel <- function(df) {
  df |>
    filter(!is.na(y1), !is.na(y2), !is.na(y3),
           !is.na(weight), is.finite(weight), weight > 0) |>
    transmute(y1 = as.integer(y1), y2 = as.integer(y2), y3 = as.integer(y3),
              weight = as.numeric(weight)) |>
    as.data.frame()
}

baseline_keep <- !is.na(df_qlfs$y1) & !is.na(df_qlfs$y2) & !is.na(df_qlfs$y3) &
  !is.na(df_qlfs$weight) & is.finite(df_qlfs$weight) & df_qlfs$weight > 0
legacy_keep <- baseline_keep &
  !is.na(df_qlfs$age1) & !is.na(df_qlfs$age2) & !is.na(df_qlfs$age3) &
  !is.na(df_qlfs$educ1) & !is.na(df_qlfs$educ2) & !is.na(df_qlfs$educ3)

samples <- list(
  legacy_table5 = .as_fmm_panel(df_qlfs[legacy_keep, , drop = FALSE])
)
if (RUN_BASELINE_SAMPLE)
  samples$baseline <- .as_fmm_panel(df_qlfs[baseline_keep, , drop = FALSE])
rm(df_qlfs)

.random_start <- function(model_type, stationary) {
  p <- list(
    theta0_A = runif(1, 0.005, 0.45),
    theta1_A = runif(1, 0.75, 0.995),
    theta0_B = runif(1, 0.005, 0.55),
    theta1_B = runif(1, 0.35, 0.95),
    phi = runif(1, 0.03, 0.97)
  )
  if (stationary) {
    p$alpha_A <- stationary_alpha(p$theta0_A, p$theta1_A)
    p$alpha_B <- stationary_alpha(p$theta0_B, p$theta1_B)
  } else {
    p$alpha_A <- runif(1, 0.05, 0.95)
    p$alpha_B <- runif(1, 0.05, 0.95)
  }
  if (model_type == "symmetric") p$pi <- runif(1, 0.001, 0.12)
  .resolve_fmm_labels(p)
}

.make_starts <- function(model_type, stationary) {
  seeds <- list(init_params_fmm(model_type, stationary))
  if (model_type == "symmetric" && stationary) {
    # Values described in the current paper, included as a replication start.
    paper <- list(theta0_A = 0.0605, theta1_A = 1 - 0.0195,
                  theta0_B = 0.0213, theta1_B = 1 - 0.0544,
                  phi = 0.4326, pi = 0.0292)
    paper$alpha_A <- stationary_alpha(paper$theta0_A, paper$theta1_A)
    paper$alpha_B <- stationary_alpha(paper$theta0_B, paper$theta1_B)
    seeds[[length(seeds) + 1L]] <- paper
  }
  while (length(seeds) < N_STARTS)
    seeds[[length(seeds) + 1L]] <- .random_start(model_type, stationary)
  seeds
}

configs <- data.frame(
  label = c("fmm_non_stat", "fmm_non_free", "fmm_sym_stat", "fmm_sym_free"),
  model_type = c("none", "none", "symmetric", "symmetric"),
  stationary = c(TRUE, FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE
)

all_fits <- list()
for (sample_name in names(samples)) {
  cat(sprintf("\n========== %s sample: N=%s ==========\n", sample_name,
              format(nrow(samples[[sample_name]]), big.mark = ",")))
  sample_fits <- list()
  for (i in seq_len(nrow(configs))) {
    cfg <- configs[i, ]
    fit <- fit_fmm_mle(
      samples[[sample_name]], cfg$model_type, cfg$stationary,
      starts = .make_starts(cfg$model_type, cfg$stationary), verbose = 1L
    )
    sample_fits[[cfg$label]] <- fit
  }
  all_fits[[sample_name]] <- sample_fits
}

.fit_row <- function(sample_name, label) {
  fit <- all_fits[[sample_name]][[label]]
  p <- fit$params
  data.frame(
    sample = sample_name, label = label,
    model_type = fit$model_type, stationary = fit$stationary,
    N = fit$sample$n, loglik = fit$loglik,
    entry_A = p$theta0_A, exit_A = 1 - p$theta1_A, alpha_A = p$alpha_A,
    entry_B = p$theta0_B, exit_B = 1 - p$theta1_B, alpha_B = p$alpha_B,
    phi = p$phi, pi = p$pi %||% NA_real_,
    weighted_entry = p$phi * p$theta0_A + (1 - p$phi) * p$theta0_B,
    weighted_exit = p$phi * (1 - p$theta1_A) +
      (1 - p$phi) * (1 - p$theta1_B),
    jacobian_rank = fit$diagnostics$probability_jacobian_rank,
    parameter_count = fit$diagnostics$parameter_count,
    identified = fit$identified,
    optimizer_validated = fit$converged,
    max_abs_score = fit$diagnostics$max_abs_score,
    min_eigenvalue = fit$diagnostics$information_min_eigenvalue,
    stringsAsFactors = FALSE
  )
}

summary_rows <- do.call(rbind, lapply(names(samples), function(sample_name) {
  do.call(rbind, lapply(configs$label, function(label) .fit_row(sample_name, label)))
}))
rownames(summary_rows) <- NULL

results_dir <- here::here("EM-baseline-ext", "output", "results")
tables_dir <- here::here("EM-baseline-ext", "output", "tables")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(all_fits, file.path(results_dir, "fit_fmm_table5_audit.rds"))
write.csv(summary_rows, file.path(results_dir, "fmm_table5_audit.csv"), row.names = FALSE)

display <- summary_rows
for (nm in c("entry_A", "exit_A", "alpha_A", "entry_B", "exit_B", "alpha_B",
             "phi", "pi", "weighted_entry", "weighted_exit")) {
  display[[nm]] <- 100 * display[[nm]]
}
cat("\n========== Direct-MLE FMM results (percent) ==========\n")
print(display, row.names = FALSE, digits = 6)

# Report how much near-optimal solutions vary. Large variation at essentially
# the same likelihood is a practical manifestation of a likelihood ridge.
cat("\n========== Top-solution dispersion ==========\n")
for (label in configs$label) {
  fit <- all_fits$legacy_table5[[label]]
  top <- head(fit$candidates, 10L)
  cat(sprintf("%s: top-10 LL range %.6f; pi range %s; phi range %.4f--%.4f\n",
              label, max(top$loglik) - min(top$loglik),
              if (all(is.na(top$pi))) "NA" else sprintf("%.4f--%.4f", min(top$pi), max(top$pi)),
              min(top$phi), max(top$phi)))
}

# A diagnostic version of Table 5. Values in unidentified columns are retained
# solely to show one numerical point on the ridge; they are not estimates with
# a unique population interpretation.
tab <- summary_rows[summary_rows$sample == "legacy_table5", ]
tab <- tab[match(configs$label, tab$label), ]
.pct <- function(x) ifelse(is.na(x), "---", sprintf("%.2f", 100 * x))
lines <- c(
  "\\begin{table}[htbp]", "\\centering",
  "\\caption{Finite mixture model: direct-likelihood replication and identification audit}",
  "\\label{tab:fmm}", "\\begin{tabular}{lcccc}", "\\toprule",
  " & No error, stat. & No error, free & Symmetric, stat. & Symmetric, free \\\\",
  "\\midrule",
  paste(c("Type A entry (\\%)", .pct(tab$entry_A)), collapse = " & ") |> paste0(" \\\\"),
  paste(c("Type A exit (\\%)", .pct(tab$exit_A)), collapse = " & ") |> paste0(" \\\\"),
  paste(c("Type A initial employment (\\%)", .pct(tab$alpha_A)), collapse = " & ") |> paste0(" \\\\"),
  paste(c("Type B entry (\\%)", .pct(tab$entry_B)), collapse = " & ") |> paste0(" \\\\"),
  paste(c("Type B exit (\\%)", .pct(tab$exit_B)), collapse = " & ") |> paste0(" \\\\"),
  paste(c("Type B initial employment (\\%)", .pct(tab$alpha_B)), collapse = " & ") |> paste0(" \\\\"),
  paste(c("Type A share $\\phi$ (\\%)", .pct(tab$phi)), collapse = " & ") |> paste0(" \\\\"),
  paste(c("Misclassification $\\pi$ (\\%)", .pct(tab$pi)), collapse = " & ") |> paste0(" \\\\"),
  paste(c("Jacobian rank / parameters", paste0(tab$jacobian_rank, "/", tab$parameter_count)), collapse = " & ") |> paste0(" \\\\"),
  paste(c("Locally identified", ifelse(tab$identified, "Yes", "No")), collapse = " & ") |> paste0(" \\\\"),
  paste(c("Optimizer validated", ifelse(tab$optimizer_validated, "Yes", "No")), collapse = " & ") |> paste0(" \\\\"),
  paste(c("Log-likelihood (millions)", sprintf("%.3f", tab$loglik / 1e6)), collapse = " & ") |> paste0(" \\\\"),
  paste(c("$N$", formatC(tab$N, format = "d", big.mark = ",")), collapse = " & ") |> paste0(" \\\\"),
  "\\bottomrule", "\\end{tabular}",
  "\\begin{minipage}{0.98\\linewidth}",
  paste0("\\footnotesize \\textit{Note:} Direct maximisation of the eight-cell observed-data likelihood. ",
         "The legacy Table 5 complete-age/education sample is used. A rank below the parameter count ",
         "means the column is locally underidentified; its displayed values are one point on a likelihood ridge ",
         "and do not have unique standard errors or a unique structural interpretation."),
  "\\end{minipage}", "\\end{table}"
)
writeLines(lines, file.path(tables_dir, "table_fmm.tex"))
cat("\nWritten diagnostic Table 5 and audit outputs.\n")
