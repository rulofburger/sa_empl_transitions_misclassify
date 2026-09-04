# Four-wave AR(1) baseline for the corrected discrete-month duration model.
# The four blocks are cached independently because a full numerical-gradient
# iteration is expensive on all sixteen latent histories.

if (!file.exists("EM-tenure/R/source_all.R")) stop("Run from project root")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("EM-tenure/R/source_all.R")
source("scripts/ingest_data_4waves_SA.R")

outdir <- paste0("EM-tenure/output/results/job_change_monthly/",
  "four_wave_ar1")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
workers <- max(1L, min(8L, as.integer(Sys.getenv(
  "FOUR_WAVE_AR1_WORKERS", "2"))))
resume <- identical(tolower(Sys.getenv("FOUR_WAVE_AR1_RESUME", "true")),
  "true")
report_only <- identical(tolower(Sys.getenv("FOUR_WAVE_AR1_REPORT_ONLY",
  "false")), "true")
cycles <- as.integer(Sys.getenv("FOUR_WAVE_AR1_CYCLES", "2"))
full_maxit <- as.integer(Sys.getenv("FOUR_WAVE_AR1_FULL_MAXIT", "0"))

prepared_file <- file.path(outdir, "prepared_cells_latest.rds")
if (resume && file.exists(prepared_file)) {
  cached_data <- readRDS(prepared_file)
  df4 <- cached_data$df4
  df3_same <- cached_data$df3_same
  sample_flow <- cached_data$sample_flow
} else {
  n_status_complete <- nrow(df_qlfs)
  prepared <- prepare_eps_estimation_data_4w(
    add_nominal_interview_months_4w(df_qlfs), drop_incomplete = TRUE)
  n_duration_complete <- nrow(prepared)
  interview4 <- paste0("interview_month", 1:4)
  df4 <- collapse_eps_cells_4w(prepared, extra_cols = interview4)
  keep3 <- c(paste0("y", 1:3), paste0("tenure", 1:3),
    paste0("timegap_cat", 1:3), paste0("interview_month", 1:3),
    "weight")
  df3_same <- collapse_eps_cells(prepared[, keep3],
    allow_zero_tenure = TRUE,
    extra_cols = paste0("interview_month", 1:3))
  sample_flow <- data.frame(
    stage = c("Complete four-wave employment status",
      "Complete state-appropriate duration reports",
      "Collapsed four-wave likelihood cells",
      "Collapsed waves 1-3 cells on same-person sample"),
    observations = c(n_status_complete, n_duration_complete,
      nrow(df4), nrow(df3_same)))
  saveRDS(list(df4 = df4, df3_same = df3_same,
    sample_flow = sample_flow), prepared_file)
}

three_wave_file <- paste0("EM-tenure/output/results/job_change_monthly/",
  "separate_duration_reliability/fits_latest.rds")
if (!file.exists(three_wave_file))
  stop("Estimate the preferred three-wave model first")
three_wave <- readRDS(three_wave_file)$extension
start <- three_wave$params
start$tenure_measurement_model <- "monthly"

evaluate4 <- function(params) e_step_eps_4w(df4, params,
  check_df = FALSE)
start_cache <- file.path(outdir, "evaluation_start_latest.rds")
start_evaluation <- if (resume && file.exists(start_cache)) {
  cached_start <- readRDS(start_cache)
  if (identical(cached_start$params, start)) cached_start$evaluation else
    evaluate4(start)
} else evaluate4(start)
saveRDS(list(params = start, evaluation = start_evaluation), start_cache)
start_fit <- list(params = start, loglik = start_evaluation$loglik,
  gamma = start_evaluation$gamma,
  job_change_posterior = start_evaluation$job_change_posterior,
  duration_reliability_posterior =
    start_evaluation$duration_reliability_posterior,
  convergence = NA_integer_, iterations = 0L,
  par_unconstrained = .pack_four_wave_preferred(start))

transition_names <- c("alpha", paste0("log_hg", 1:5),
  paste0("log_hd", 1:5))
reporting_names <- c("pi", "eps", "eps_d", "job_change",
  "tenure_reliability_shift", "timegap_reliability_shift")
calendar_names <- c("tenure_heaping", "tenure_year_revision",
  "tenure_clean_anchor_revision", paste0("start_month_logit", 1:11),
  "tenure_exact_anchor_retention", "tenure_local_revision")

fit_stage <- function(label, incumbent, free_names, max_parameter_move) {
  path <- file.path(outdir, paste0("fit_", label, "_latest.rds"))
  if (resume && file.exists(path)) {
    fit <- readRDS(path)
    # Early checkpoints incorrectly used code zero for an accepted step.
    # A bounded single step never establishes optimizer convergence.
    if (!is.null(fit$backtrack)) {
      if (identical(fit$convergence, 0L))
        fit$original_step_acceptance_code <- fit$convergence
      fit$convergence <- 1L
      fit$converged <- FALSE
      fit$step_accepted <- any(fit$backtrack$selected & fit$backtrack$scale > 0)
      saveRDS(fit, path)
    }
    return(fit)
  }
  if (report_only) return(NULL)
  message("Four-wave AR(1) stage: ", label)
  fit <- four_wave_gradient_step(df4, incumbent$params,
    workers = workers, free_names = free_names,
    max_parameter_move = max_parameter_move)
  fit$stage <- label
  saveRDS(fit, path)
  fit
}

accept <- function(incumbent, candidate) {
  if (is.finite(candidate$loglik) &&
      candidate$loglik >= incumbent$loglik - 1e-5) candidate else incumbent
}

stages <- list(start = start_fit)
incumbent <- start_fit
blocks <- list(transition = transition_names, reporting = reporting_names,
  calendar = calendar_names)
for (cycle in seq_len(cycles)) for (block in names(blocks)) {
  label <- paste0(block, if (cycle == 1L) "" else cycle)
  candidate <- fit_stage(label, incumbent, blocks[[block]], .25 / cycle)
  if (is.null(candidate)) next
  stages[[label]] <- candidate
  incumbent <- accept(incumbent, candidate)
}
if (!report_only && full_maxit > 0L) {
  candidate <- fit_eps_piecewise_calendar_revision_monthly_4w(df4,
    incumbent$params, maxit = full_maxit, workers = workers,
    gradient_scheme = "central", verbose = 1L)
  candidate$objective_function <- NULL
  saveRDS(candidate, file.path(outdir, "fit_full_latest.rds"))
  stages$full <- candidate
  incumbent <- accept(incumbent, candidate)
}
best <- incumbent

summarise_fit <- function(label, fit) {
  rates <- duration_weighted_transition_rates_4w(df4, fit)[1L, ]
  params <- fit$params
  shifts <- .duration_reliability_shifts(params)
  posterior_job_change <- sum(df4$weight *
    fit$job_change_posterior$expected_changes) /
    sum(df4$weight * fit$job_change_posterior$opportunities)
  data.frame(stage = label, loglik = fit$loglik,
    gain_over_start = fit$loglik - start_fit$loglik,
    convergence = fit$convergence, iterations = fit$iterations,
    status = if (identical(fit$convergence, 0L))
      "Optimizer stopped; independent score check required" else
      "Provisional; convergence not established",
    entry_rate = rates$entry_rate, exit_rate = rates$exit_rate,
    initial_employment = params$alpha,
    status_misclassification = params$pi,
    tenure_contamination_midpoint = params$eps,
    timegap_contamination_midpoint = params$eps_d,
    posterior_job_change = posterior_job_change,
    tenure_reliability_shift = unname(shifts["tenure"]),
    timegap_reliability_shift = unname(shifts["timegap"]))
}
stage_summary <- do.call(rbind, Map(summarise_fit, names(stages), stages))
write.csv(sample_flow, file.path(outdir, "sample_flow_latest.csv"),
  row.names = FALSE)
write.csv(stage_summary, file.path(outdir, "stage_summary_latest.csv"),
  row.names = FALSE)
cat("\nFour-wave AR(1) checkpoints (provisional)\n")
print(stage_summary[, c("stage", "loglik", "entry_rate", "exit_rate",
  "status_misclassification", "initial_employment", "convergence")],
  row.names = FALSE, digits = 8)

# Hold the model fixed to distinguish the four-wave sample restriction from
# the extra-wave likelihood contribution. Re-estimation of the same-person
# three-wave comparator can be added after the four-wave baseline is stable.
same3_evaluation <- e_step_eps(df3_same, start, check_df = FALSE,
  suff_stats = FALSE)
same3_fit <- list(params = start, gamma = same3_evaluation$gamma)
same3_rates <- duration_weighted_transition_rates(df3_same, same3_fit)[1L, ]
options(sa_empl_transitions.qlfs_3wave_panel = "df_qlfs_A.rds",
  sa_empl_transitions.preserve_zero_tenure = TRUE)
source("scripts/ingest_data_3waves_SA.R")
full3_prepared <- prepare_eps_estimation_data(
  add_nominal_interview_months(df_qlfs), allow_zero_tenure = TRUE)
n_full3 <- nrow(full3_prepared)
full3_cells <- collapse_eps_cells(full3_prepared,
  allow_zero_tenure = TRUE,
  extra_cols = paste0("interview_month", 1:3))
full3_evaluation <- e_step_eps(full3_cells, start, check_df = FALSE,
  suff_stats = FALSE)
full3_rates <- duration_weighted_transition_rates(full3_cells,
  list(params = start, gamma = full3_evaluation$gamma))[1L, ]
sample_comparison <- data.frame(
  sample = c("Three-wave estimation sample",
    "Four-wave-linked people, waves 1-3",
    "Four-wave-linked people, waves 1-4"),
  parameters = c("Three-wave optimum", "Three-wave optimum",
    "Four-wave AR(1), provisional bounded refinement"),
  observations = c(n_full3,
    attr(df3_same, "n_original"), attr(df4, "n_original")),
  entry_rate = c(full3_rates$entry_rate, same3_rates$entry_rate,
    duration_weighted_transition_rates_4w(df4, best)$entry_rate[1L]),
  exit_rate = c(full3_rates$exit_rate, same3_rates$exit_rate,
    duration_weighted_transition_rates_4w(df4, best)$exit_rate[1L]))

cat("\nFour-wave sample construction\n")
print(sample_flow, row.names = FALSE)
cat("\nFour-wave AR(1) optimization\n")
print(stage_summary, row.names = FALSE, digits = 9)
cat("\nSample and wave comparison\n")
print(sample_comparison, row.names = FALSE, digits = 9)

write.csv(sample_flow, file.path(outdir, "sample_flow_latest.csv"),
  row.names = FALSE)
write.csv(stage_summary, file.path(outdir, "stage_summary_latest.csv"),
  row.names = FALSE)
write.csv(sample_comparison,
  file.path(outdir, "sample_comparison_latest.csv"), row.names = FALSE)
saveRDS(list(best = best, stages = stages, sample_flow = sample_flow,
  stage_summary = stage_summary, sample_comparison = sample_comparison,
  same_person_three_wave_reestimated = FALSE,
  independent_convergence_check_passed = FALSE,
  created_at = Sys.time(), session_info = sessionInfo()),
  file.path(outdir, "fits_latest.rds"))
