# Keep repaired empirical estimates separate from the historical likelihood.
# Run from the repository root after refine_four_wave_ar1_fast.R completes.
source("EM-tenure/R/source_all.R")
source("EM-tenure/R/four_wave_fast.R")
base <- "EM-tenure/output/results/job_change_monthly/four_wave_ar1"
outdir <- file.path(base, "continuous_clock_empirical")
fingerprint <- four_wave_fast_source_fingerprint()
gates <- readRDS(file.path(base, "recovery_continuous_clock/recovery_status_latest.rds"))
stopifnot(identical(gates$source_md5, fingerprint),
  isTRUE(gates$ready_for_recovery_fits), isTRUE(gates$recovery_fits_all_converged))
df <- readRDS(file.path(base, "prepared_cells_latest.rds"))$df4
legacy <- readRDS(file.path(base, "converged_comparison_latest.rds"))$best
paths <- list.files(base, "^fit_fast_.*_latest[.]rds$", full.names = TRUE)
fits <- lapply(paths, readRDS)
keep <- vapply(fits, function(f)
  identical(f$params$timegap_clock_model, "continuous_joint") &&
    identical(f$source_md5, fingerprint), logical(1))
fits <- fits[keep]; paths <- paths[keep]
if (!length(fits)) stop("No current-source repaired empirical fits")
names(fits) <- sub("_latest[.]rds$", "", sub("^fit_fast_", "", basename(paths)))
eligible <- vapply(fits, function(f) isTRUE(f$converged) &&
  isTRUE(f$independent_score_pass) && f$convergence == 0L &&
  length(f$projected_score) == 33L &&
  all(is.finite(f$projected_score)) &&
  max(abs(f$projected_score)) <= 1e-5, logical(1))
if (!any(eligible)) stop("No repaired fit passes independent convergence checks")
best_index <- which(eligible)[which.max(vapply(fits[eligible], `[[`, numeric(1), "loglik"))]
best <- fits[[best_index]]
stopifnot(isTRUE(legacy$converged), nrow(best$gamma) == nrow(df),
  nrow(legacy$gamma) == nrow(df), max(abs(rowSums(best$gamma)-1)) < 1e-10)
all_fits <- c(list(legacy = legacy), fits)
rate_rows <- lapply(names(all_fits), function(label) {
  f <- all_fits[[label]]
  r <- duration_weighted_transition_rates_4w(df, f)
  data.frame(model = label, r, row.names = NULL)
})
intervals <- do.call(rbind, rate_rows)
summary <- do.call(rbind, lapply(names(all_fits), function(label) {
  f <- all_fits[[label]]
  r <- intervals[intervals$model == label & intervals$wave == 0L, ]
  data.frame(model = label, n = attr(df, "n_original"), cells = nrow(df),
    entry_rate = r$entry_rate, exit_rate = r$exit_rate,
    misclassification_rate = f$params$pi, initial_employment = f$params$alpha,
    loglik = f$loglik, projected_score = max(abs(f$projected_score)),
    optimizer_code = f$convergence, converged = isTRUE(f$converged))
}))
selected <- summary[match(c("legacy", names(fits)[best_index]), summary$model), ]
measures <- c("entry_rate", "exit_rate", "misclassification_rate", "initial_employment")
comparison <- data.frame(measure = measures,
  legacy = as.numeric(selected[1, measures]),
  repaired = as.numeric(selected[2, measures]))
comparison$change_percentage_points <- 100 * (comparison$repaired - comparison$legacy)
hazards <- data.frame(hazard = rep(c("Entry", "Exit"), each = 5),
  duration = rep(c("0--3 months", "3--12 months", "1--3 years", "3--5 years", "5+ years"), 2),
  legacy = c(legacy$params$lambda_d, legacy$params$lambda_g),
  repaired = c(best$params$lambda_d, best$params$lambda_g))
coefficients <- data.frame(parameter = names(best$par_unconstrained),
  legacy = legacy$par_unconstrained[names(best$par_unconstrained)],
  repaired = best$par_unconstrained, row.names = NULL)

# Verify saved empirical posteriors and the repaired objective at the legacy
# parameter vector. Neither cross-formula likelihood difference is an LR test.
data <- prepare_four_wave_kernel_data(df)
reevaluated <- evaluate_four_wave_monthly_fast(data, best$params, threads = 8L)
stopifnot(abs(reevaluated$loglik - best$loglik) < 1e-6,
  max(abs(reevaluated$gamma - best$gamma)) < 1e-10)
old_at_new_law <- evaluate_four_wave_monthly_fast(data, legacy$params,
  timegap_clock = "continuous_joint", posterior = FALSE, threads = 8L)$loglik
stopifnot(best$loglik >= old_at_new_law - 1e-6)
# An independent full-sample R evaluation at the NEW fitted parameter vector.
reference <- e_step_eps_4w(df, best$params, check_df = FALSE, exact_risk = TRUE)
checks <- c(loglik_difference = best$loglik - reference$loglik,
  max_posterior_difference = max(abs(best$gamma - reference$gamma)))
stopifnot(abs(checks[1]) < 1e-6, checks[2] < 1e-10)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(summary, file.path(outdir, "fit_summary_latest.csv"), row.names = FALSE)
write.csv(comparison, file.path(outdir, "rate_comparison_latest.csv"), row.names = FALSE)
write.csv(intervals, file.path(outdir, "interval_rates_latest.csv"), row.names = FALSE)
write.csv(hazards, file.path(outdir, "hazards_latest.csv"), row.names = FALSE)
write.csv(coefficients, file.path(outdir, "transformed_coefficients_latest.csv"), row.names = FALSE)
saveRDS(list(best = best, best_label = names(fits)[best_index],
  summary = summary, comparison = comparison, hazards = hazards, intervals = intervals,
  legacy_loglik_under_repaired_law = old_at_new_law,
  repaired_objective_gain = best$loglik - old_at_new_law,
  reference_checks = checks, input_paths = paths,
  input_md5 = tools::md5sum(c(paths, file.path(base, "prepared_cells_latest.rds"))),
  source_md5 = fingerprint, n_parameters = length(best$par_unconstrained),
  empirical_fit_changed = TRUE, standard_errors_estimated = FALSE,
  same_person_three_wave_reestimated = FALSE, ar2_estimated = FALSE,
  session_info = sessionInfo()), file.path(outdir, "comparison_latest.rds"))
print(summary, row.names = FALSE, digits = 9)
print(comparison, row.names = FALSE, digits = 9)
print(hazards, row.names = FALSE, digits = 7)
print(checks)
cat("Gain on repaired objective from legacy parameters:", best$loglik - old_at_new_law, "\n")
